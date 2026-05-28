"""
Bulk RNA-seq Explorer
Version: bulk_rnaseq_explorer_v1_4

Scope for v1.4:
- Clean Streamlit product UI for count-matrix upload and sample grouping.
- Detect whether the uploaded gene IDs are Ensembl IDs, gene symbols, mixed, or unclear.
- Convert mouse Ensembl IDs to gene symbols when a local mapping can be parsed.
- Merge duplicated processed gene symbols by summing raw counts.
- Produce a clean processed count matrix for future QC.
- Add QC Overview data summaries without plots or normalization.
- Show concise internal readiness warnings without exposing a full validation/debug table.

To reduce Streamlit toolbar/menu visibility, users may create `.streamlit/config.toml` with:

[client]
toolbarMode = "minimal"

Explicitly out of scope:
- DESeq2, Rscript, DEG analysis, QC plots, PCA, heatmap, volcano plot,
  GSEA, ORA, pathway enrichment, cloud storage, login, Duke DCC, SLURM.
"""

from __future__ import annotations

import csv
import hashlib
import json
import re
from io import StringIO
from pathlib import Path
from typing import Any

import pandas as pd
import streamlit as st


APP_VERSION = "bulk_rnaseq_explorer_v1_4"

TEMPLATE_TEXT = """EnsemblID / Gene_symbol\tSample_1\tSample_2\tSample_3\tSample_4
ENSMUSG00000000001\t120\t98\t115\t130
ENSMUSG00000000028\t0\t4\t1\t3
Cxcl1\t50\t80\t320\t400
Actb\t10000\t9800\t10300\t9900
"""


def init_session_state() -> None:
    """Initialize session-scoped project state."""
    defaults: dict[str, Any] = {
        "app_version": APP_VERSION,
        "raw_counts_df": None,
        "processed_counts_df": None,
        "counts_file_name": None,
        "counts_file_signature": None,
        "gene_id_column": None,
        "sample_columns": [],
        "sample_group_assignments": {},
        "group_1_name": "Control",
        "group_2_name": "Treatment",
        "gene_map_df": None,
        "gene_map_status": {
            "detected": False,
            "parsed": False,
            "path": None,
            "n_mappings": 0,
            "message": "Not loaded",
        },
        "gene_id_mode": "unknown",
        "conversion_summary": default_conversion_summary(),
        "duplicate_summary_df": pd.DataFrame(columns=["Gene", "Original IDs", "Number of duplicated rows"]),
        "qc_results": None,
        "deg_results": None,
        "pathway_results": None,
        "plots": None,
    }
    for key, value in defaults.items():
        if key not in st.session_state:
            st.session_state[key] = value


def default_conversion_summary() -> dict[str, Any]:
    """Return default gene conversion and duplicate-merge summary."""
    return {
        "gene_id_mode": "unknown",
        "mapping_source": "Not available",
        "total_genes": 0,
        "converted_to_gene_symbols": 0,
        "unconverted_genes": 0,
        "duplicated_gene_symbols": 0,
        "processed_genes": 0,
    }


def clear_count_matrix_state() -> None:
    """Clear all count-dependent state after uploader file removal."""
    st.session_state["raw_counts_df"] = None
    st.session_state["processed_counts_df"] = None
    st.session_state["counts_file_name"] = None
    st.session_state["counts_file_signature"] = None
    st.session_state["gene_id_column"] = None
    st.session_state["sample_columns"] = []
    st.session_state["sample_group_assignments"] = {}
    st.session_state["gene_id_mode"] = "unknown"
    st.session_state["conversion_summary"] = default_conversion_summary()
    st.session_state["duplicate_summary_df"] = pd.DataFrame(columns=["Gene", "Original IDs", "Number of duplicated rows"])
    reset_analysis_state()


def reset_analysis_state() -> None:
    """Clear future analysis outputs when upstream data or assignments change."""
    st.session_state["qc_results"] = None
    st.session_state["deg_results"] = None
    st.session_state["pathway_results"] = None
    st.session_state["plots"] = None


def read_count_matrix_file(uploaded_file) -> pd.DataFrame:
    """Read an uploaded count matrix; tab-delimited input is preferred."""
    if uploaded_file is None:
        raise ValueError("No uploaded file was provided.")

    suffix = Path(uploaded_file.name).suffix.lower()
    uploaded_file.seek(0)

    if suffix == ".csv":
        return pd.read_csv(uploaded_file, sep=",")
    if suffix in {".tsv", ".txt"}:
        return pd.read_csv(uploaded_file, sep="\t")

    uploaded_file.seek(0)
    sample_bytes = uploaded_file.read(4096)
    uploaded_file.seek(0)
    try:
        sample_text = sample_bytes.decode("utf-8", errors="replace")
        dialect = csv.Sniffer().sniff(sample_text, delimiters=",\t;")
        return pd.read_csv(uploaded_file, sep=dialect.delimiter)
    except Exception:
        uploaded_file.seek(0)
        return pd.read_csv(uploaded_file, sep=None, engine="python")


def get_sample_columns(counts_df: pd.DataFrame | None, gene_id_column: str | None) -> list[str]:
    """Return count matrix columns other than the selected gene ID column."""
    if counts_df is None or gene_id_column is None:
        return []
    columns = [str(column) for column in counts_df.columns]
    if str(gene_id_column) not in columns:
        return []
    return [column for column in columns if column != str(gene_id_column)]


def detect_local_resources() -> pd.DataFrame:
    """Detect local files silently for internal loading."""
    cwd = Path.cwd()
    resource_specs = [
        ("source_mapping/mouse_ensembl_to_symbol", cwd / "source_mapping" / "mouse_ensembl_to_symbol"),
        ("source mapping/mouse_ensembl_to_symbol", cwd / "source mapping" / "mouse_ensembl_to_symbol"),
        ("mouse_ensembl_to_symbol.js", cwd / "mouse_ensembl_to_symbol.js"),
        ("source_mapping/", cwd / "source_mapping"),
        ("source mapping/", cwd / "source mapping"),
        ("database_raw/", cwd / "database_raw"),
    ]
    return pd.DataFrame(
        [
            {"Resource": name, "Path": str(path), "Detected": path.exists()}
            for name, path in resource_specs
        ]
    )


def _mapping_dataframe_from_object(data: Any) -> pd.DataFrame | None:
    """Convert common mapping object shapes to ensembl_id/gene_symbol columns."""
    if isinstance(data, dict):
        rows = []
        for ensembl_id, value in data.items():
            if isinstance(value, str):
                rows.append({"ensembl_id": ensembl_id, "gene_symbol": value})
            elif isinstance(value, dict):
                symbol = (
                    value.get("gene_symbol")
                    or value.get("symbol")
                    or value.get("Gene")
                    or value.get("gene_name")
                    or value.get("name")
                )
                if symbol:
                    rows.append({"ensembl_id": ensembl_id, "gene_symbol": symbol})
        return pd.DataFrame(rows) if rows else None

    if isinstance(data, list):
        ensembl_keys = ["ensembl_id", "EnsemblID", "ensembl", "gene_id"]
        symbol_keys = ["gene_symbol", "symbol", "Gene", "gene_name"]
        rows = []
        for item in data:
            if not isinstance(item, dict):
                continue
            ensembl_id = next((item.get(key) for key in ensembl_keys if item.get(key)), None)
            gene_symbol = next((item.get(key) for key in symbol_keys if item.get(key)), None)
            if ensembl_id and gene_symbol:
                rows.append({"ensembl_id": ensembl_id, "gene_symbol": gene_symbol})
        return pd.DataFrame(rows) if rows else None

    return None


def _extract_json_like_object_from_js(text: str) -> Any:
    """Extract the first JSON-like object from a JavaScript mapping file."""
    assignment_match = re.search(r"=\s*(\{.*\})\s*;?\s*$", text, flags=re.DOTALL)
    if assignment_match:
        object_text = assignment_match.group(1)
    else:
        first_brace = text.find("{")
        last_brace = text.rfind("}")
        if first_brace < 0 or last_brace <= first_brace:
            raise ValueError("No JSON-like mapping object was found.")
        object_text = text[first_brace : last_brace + 1]
    return json.loads(object_text.rstrip(";").strip())


def parse_gene_map_file(path: Path) -> tuple[pd.DataFrame | None, dict[str, Any]]:
    """Parse a local Ensembl-to-symbol mapping file if possible."""
    status = {
        "detected": path.exists(),
        "parsed": False,
        "path": str(path),
        "n_mappings": 0,
        "message": "File not detected",
    }
    if not path.exists():
        return None, status

    try:
        text = path.read_text(encoding="utf-8", errors="replace")
        if path.suffix.lower() == ".js":
            data = _extract_json_like_object_from_js(text)
        else:
            try:
                data = json.loads(text)
            except json.JSONDecodeError:
                data = _extract_json_like_object_from_js(text)

        mapping_df = _mapping_dataframe_from_object(data)
        if mapping_df is None or mapping_df.empty:
            status["message"] = "Detected but not parsed."
            return None, status

        mapping_df = mapping_df.astype("string").dropna().drop_duplicates()
        mapping_df["ensembl_base"] = mapping_df["ensembl_id"].str.replace(r"\.\d+$", "", regex=True)
        status.update(
            {
                "parsed": True,
                "n_mappings": int(mapping_df.shape[0]),
                "message": "Parsed successfully.",
            }
        )
        return mapping_df, status
    except Exception as exc:
        status["message"] = f"Detected but not parsed: {exc}"
        return None, status


def load_local_gene_map() -> tuple[pd.DataFrame | None, dict[str, Any]]:
    """Load local gene map using source mapping first and JS as fallback."""
    cwd = Path.cwd()
    candidates = [
        cwd / "source_mapping" / "mouse_ensembl_to_symbol",
        cwd / "source mapping" / "mouse_ensembl_to_symbol",
        cwd / "source_mapping" / "mouse_ensembl_to_symbol.json",
        cwd / "source mapping" / "mouse_ensembl_to_symbol.json",
        cwd / "mouse_ensembl_to_symbol.js",
    ]

    first_detected_status = None
    for path in candidates:
        if not path.exists():
            continue
        mapping_df, status = parse_gene_map_file(path)
        if first_detected_status is None:
            first_detected_status = status
        if mapping_df is not None:
            return mapping_df, status

    if first_detected_status is not None:
        return None, first_detected_status

    return None, {
        "detected": False,
        "parsed": False,
        "path": None,
        "n_mappings": 0,
        "message": "No local mapping available.",
    }


def detect_gene_id_mode(gene_ids: pd.Series) -> str:
    """Classify gene IDs as ensembl, symbol, mixed, or unknown."""
    values = gene_ids.dropna().astype("string").str.strip()
    values = values[values.ne("")]
    if values.empty:
        return "unknown"

    ensembl_mask = values.str.match(r"^ENSMUSG\d+(?:\.\d+)?$", case=False, na=False)
    symbol_mask = values.str.match(r"^[A-Za-z][A-Za-z0-9_.-]*$", na=False) & ~ensembl_mask
    ensembl_fraction = float(ensembl_mask.mean())
    symbol_fraction = float(symbol_mask.mean())

    if ensembl_fraction >= 0.6:
        return "ensembl"
    if symbol_fraction >= 0.6:
        return "symbol"
    if ensembl_fraction >= 0.2 and symbol_fraction >= 0.2:
        return "mixed"
    return "unknown"


def merge_duplicate_genes_by_sum(
    processed_df: pd.DataFrame,
    gene_column: str,
    sample_columns: list[str],
) -> pd.DataFrame:
    """Merge duplicated genes by summing raw counts across sample columns."""
    merge_df = processed_df[[gene_column] + sample_columns].copy()
    for sample in sample_columns:
        merge_df[sample] = pd.to_numeric(merge_df[sample], errors="coerce").fillna(0)
    merged = merge_df.groupby(gene_column, sort=False, as_index=False)[sample_columns].sum()
    return merged.rename(columns={gene_column: "Gene"})


def process_gene_symbols_and_merge_duplicates(
    counts_df: pd.DataFrame,
    gene_id_column: str,
    gene_map_df: pd.DataFrame | None,
    gene_map_status: dict[str, Any],
) -> tuple[pd.DataFrame, dict[str, Any], pd.DataFrame]:
    """Convert Ensembl IDs where possible, then merge duplicated processed genes."""
    sample_columns = get_sample_columns(counts_df, gene_id_column)
    gene_ids = counts_df[gene_id_column].astype("string").str.strip()
    gene_id_mode = detect_gene_id_mode(gene_ids)
    total_genes = int(counts_df.shape[0])
    converted_count = 0
    unconverted_count = 0

    processed_gene_ids = gene_ids.copy()
    mapping_source = "Not available"

    if gene_id_mode == "ensembl":
        if gene_map_df is not None and gene_map_status.get("parsed"):
            mapping_source = gene_map_status.get("path") or "Local gene map"
            mapping_by_exact = dict(zip(gene_map_df["ensembl_id"].astype(str), gene_map_df["gene_symbol"].astype(str)))
            mapping_by_base = dict(zip(gene_map_df["ensembl_base"].astype(str), gene_map_df["gene_symbol"].astype(str)))
            mapped_symbols = []
            for gene_id in gene_ids:
                gene_id_str = str(gene_id)
                gene_base = re.sub(r"\.\d+$", "", gene_id_str)
                symbol = mapping_by_exact.get(gene_id_str) or mapping_by_base.get(gene_base)
                mapped_symbols.append(symbol)
            mapped_series = pd.Series(mapped_symbols, index=gene_ids.index, dtype="string")
            converted_mask = mapped_series.notna() & mapped_series.ne("")
            processed_gene_ids = gene_ids.mask(converted_mask, mapped_series)
            converted_count = int(converted_mask.sum())
            unconverted_count = int(total_genes - converted_count)
        else:
            mapping_source = gene_map_status.get("message", "Local mapping was not parsed.")
            unconverted_count = total_genes
    elif gene_id_mode == "symbol":
        mapping_source = "Input gene symbols - conversion skipped"
    else:
        mapping_source = "Mixed or unclear gene IDs - conversion skipped"

    working_df = counts_df.copy()
    working_df["Gene"] = processed_gene_ids.fillna(gene_ids).astype(str)
    working_df["_original_gene_id"] = gene_ids.astype(str)

    duplicate_mask = working_df["Gene"].duplicated(keep=False)
    if duplicate_mask.any():
        duplicate_summary_df = (
            working_df.loc[duplicate_mask]
            .groupby("Gene", sort=False)
            .agg(
                **{
                    "Original IDs": ("_original_gene_id", lambda values: ", ".join(map(str, values))),
                    "Number of duplicated rows": ("_original_gene_id", "size"),
                }
            )
            .reset_index()
        )
    else:
        duplicate_summary_df = pd.DataFrame(columns=["Gene", "Original IDs", "Number of duplicated rows"])

    processed_counts_df = merge_duplicate_genes_by_sum(working_df, "Gene", sample_columns)
    conversion_summary = {
        "gene_id_mode": gene_id_mode,
        "mapping_source": mapping_source,
        "total_genes": total_genes,
        "converted_to_gene_symbols": converted_count,
        "unconverted_genes": unconverted_count,
        "duplicated_gene_symbols": int(duplicate_summary_df.shape[0]),
        "processed_genes": int(processed_counts_df.shape[0]),
    }
    return processed_counts_df, conversion_summary, duplicate_summary_df


def compute_basic_readiness() -> dict[str, Any]:
    """Compute concise readiness state for QC and future DEG."""
    raw_counts_df = st.session_state["raw_counts_df"]
    sample_columns = st.session_state["sample_columns"]
    assignments = st.session_state["sample_group_assignments"]
    group_1 = st.session_state["group_1_name"]
    group_2 = st.session_state["group_2_name"]

    errors: list[str] = []
    warnings: list[str] = []

    if raw_counts_df is None:
        return {"qc": "Not ready", "deg": "Not ready", "messages": ["Upload a count matrix first."]}

    if len(sample_columns) < 2:
        errors.append("At least 2 sample columns are required.")

    numeric_counts = raw_counts_df[sample_columns].apply(pd.to_numeric, errors="coerce") if sample_columns else pd.DataFrame()
    for sample in sample_columns:
        original_text = raw_counts_df[sample].astype("string").str.strip()
        non_missing_original = raw_counts_df[sample].notna() & original_text.ne("")
        if (non_missing_original & numeric_counts[sample].isna()).any():
            errors.append(f"Sample column `{sample}` contains non-numeric values.")

    if not numeric_counts.empty:
        if (numeric_counts < 0).sum().sum() > 0:
            errors.append("Negative count values were detected.")
        valid_mask = numeric_counts.notna()
        non_integer_entries = int((valid_mask & ((numeric_counts % 1) != 0)).sum().sum())
        if non_integer_entries > 0:
            warnings.append("Non-integer counts detected; DESeq2 expects raw integer counts.")
        library_sizes = numeric_counts.fillna(0).sum(axis=0)
        zero_library_samples = [sample for sample, total in library_sizes.items() if total == 0]
        if zero_library_samples:
            errors.append(f"Zero-library samples detected: {', '.join(map(str, zero_library_samples[:8]))}.")

    group_1_count = sum(1 for sample in sample_columns if assignments.get(sample) == group_1)
    group_2_count = sum(1 for sample in sample_columns if assignments.get(sample) == group_2)
    unassigned_count = sum(1 for sample in sample_columns if assignments.get(sample, "Unassigned") not in {group_1, group_2})

    if group_1 == group_2:
        errors.append("Group 1 and Group 2 names must be different.")
    if group_1_count == 0 or group_2_count == 0:
        errors.append("Both groups must contain at least one sample.")
    if group_1_count == 1 or group_2_count == 1:
        warnings.append("DESeq2 requires biological replicates for reliable inference.")
    if unassigned_count > 0:
        warnings.append(f"{unassigned_count} samples are unassigned.")

    qc_status = "Ready" if not any(msg for msg in errors if "Group" not in msg and "groups" not in msg) else "Not ready"
    if errors:
        deg_status = "Not ready"
    elif warnings:
        deg_status = "Warning"
    else:
        deg_status = "Ready"

    return {"qc": qc_status, "deg": deg_status, "messages": errors + warnings}


def compute_qc_summary(processed_counts_df: pd.DataFrame, sample_columns: list[str]) -> dict[str, Any]:
    """Compute sample-level and gene-level QC summaries from processed counts."""
    numeric_counts = processed_counts_df[sample_columns].apply(pd.to_numeric, errors="coerce").fillna(0)

    sample_qc_df = pd.DataFrame(
        {
            "Sample": sample_columns,
            "Library size": [float(numeric_counts[sample].sum()) for sample in sample_columns],
            "Detected genes": [int((numeric_counts[sample] > 0).sum()) for sample in sample_columns],
            "Zero-count genes": [int((numeric_counts[sample] == 0).sum()) for sample in sample_columns],
            "Zero fraction": [float((numeric_counts[sample] == 0).mean()) for sample in sample_columns],
            "Median count": [float(numeric_counts[sample].median()) for sample in sample_columns],
            "Mean count": [float(numeric_counts[sample].mean()) for sample in sample_columns],
        }
    )

    gene_totals = numeric_counts.sum(axis=1)
    gene_qc_summary = {
        "Total processed genes": int(processed_counts_df.shape[0]),
        "Expressed genes": int((gene_totals > 0).sum()),
        "All-zero genes": int((gene_totals == 0).sum()),
        "Low-count genes": int((gene_totals < 10).sum()),
        "Constant genes": int(numeric_counts.nunique(axis=1, dropna=False).eq(1).sum()),
    }
    return {"sample_qc_df": sample_qc_df, "gene_qc_summary": gene_qc_summary}


def format_path_for_display(path: str | None, max_length: int = 60) -> str:
    """Shorten long paths for stable page layout."""
    if not path:
        return "Not available"
    path_text = str(path)
    if len(path_text) <= max_length:
        return path_text
    parts = re.split(r"([\\/])", path_text)
    tail = "".join(parts[-5:]) if len(parts) >= 5 else path_text[-max_length + 3 :]
    return f"...{tail}"


def render_inline_badges(items: list[str], max_items: int = 30, separator: str = "/") -> None:
    """Render list items as separate inline code badges."""
    visible_items = items[:max_items]
    pieces = []
    for index, item in enumerate(visible_items):
        if index > 0:
            pieces.append(f" {separator} ")
        pieces.append(f"`{item}`")
    if len(items) > max_items:
        pieces.append(f" ... and {len(items) - max_items} more")
    st.write("".join(pieces) if pieces else "`None`")


def render_info_line(label: str, value: Any) -> None:
    """Render one compact label/value line with code-style value."""
    st.write(f"{label}: `{value}`")


def render_sidebar() -> None:
    """Render workflow-only sidebar."""
    st.sidebar.title("Workflow")
    st.sidebar.success("1. Upload Count Matrix")
    st.sidebar.success("2. Assign Sample Groups")
    st.sidebar.success("3. QC Overview")
    st.sidebar.caption("4. DEG Analysis - Coming soon / Locked")
    st.sidebar.caption("5. Visualization - Coming soon / Locked")
    st.sidebar.caption("6. Pathway Analysis - Coming soon / Locked")
    st.sidebar.caption("7. Export - Coming soon / Locked")


def _uploaded_file_signature(uploaded_file) -> str:
    uploaded_file.seek(0)
    content = uploaded_file.getvalue()
    uploaded_file.seek(0)
    return hashlib.sha256(content).hexdigest()


def _reprocess_uploaded_counts() -> None:
    raw_counts_df = st.session_state["raw_counts_df"]
    gene_id_column = st.session_state["gene_id_column"]
    if raw_counts_df is None or gene_id_column is None:
        return
    processed_counts_df, conversion_summary, duplicate_summary_df = process_gene_symbols_and_merge_duplicates(
        raw_counts_df,
        gene_id_column,
        st.session_state["gene_map_df"],
        st.session_state["gene_map_status"],
    )
    st.session_state["processed_counts_df"] = processed_counts_df
    st.session_state["conversion_summary"] = conversion_summary
    st.session_state["duplicate_summary_df"] = duplicate_summary_df
    st.session_state["gene_id_mode"] = conversion_summary["gene_id_mode"]
    reset_analysis_state()


def render_gene_symbol_conversion_section() -> None:
    """Render concise gene ID conversion status using compact info lines."""
    st.markdown("### Gene Symbols Conversion")
    summary = st.session_state["conversion_summary"]
    mode = summary.get("gene_id_mode", "unknown")

    if mode in {"mixed", "unknown"} and st.session_state["raw_counts_df"] is not None:
        st.warning(
            "Gene ID format appears mixed or unclear. Conversion was skipped. "
            "Please confirm whether the first column contains Ensembl IDs or gene symbols."
        )

    render_info_line("Mapping source", format_path_for_display(summary.get("mapping_source")))
    render_info_line("Total genes", f"{summary.get('total_genes', 0):,}")
    render_info_line("Converted to gene symbols", f"{summary.get('converted_to_gene_symbols', 0):,}")
    render_info_line("Unconverted genes", f"{summary.get('unconverted_genes', 0):,}")


def render_duplicate_gene_section() -> None:
    """Render duplicate gene-symbol merge status."""
    st.markdown("### Duplicated gene symbols")
    st.write("Duplicated gene symbols were merged by summing raw counts.")

    summary = st.session_state["conversion_summary"]
    duplicated_count = int(summary.get("duplicated_gene_symbols", 0))
    processed_genes = int(summary.get("processed_genes", 0))

    render_info_line("Duplicated gene symbols", f"{duplicated_count:,}")
    render_info_line("Processed genes after merge", f"{processed_genes:,}")

    duplicate_summary_df = st.session_state["duplicate_summary_df"]
    if duplicated_count > 0:
        with st.expander("View duplicated genes"):
            st.caption("Original IDs are the gene identifiers before conversion or merging.")
            if duplicate_summary_df.shape[0] > 500:
                st.info("Showing first 500 duplicated gene symbols.")
            st.dataframe(duplicate_summary_df.head(500), use_container_width=True, hide_index=True)
    elif st.session_state["raw_counts_df"] is not None:
        st.success("No duplicated processed gene symbols detected.")


def render_upload_count_matrix_tab() -> None:
    """Render simplified upload page."""
    st.subheader("Upload Count Matrix")
    st.write(
        "Upload a raw count matrix as a tab-delimited file (.tsv or .txt). "
        "CSV is also accepted. The first row should contain sample names."
    )

    st.markdown("### Sample template")
    st.dataframe(pd.read_csv(StringIO(TEMPLATE_TEXT), sep="\t"), use_container_width=True, hide_index=True)

    uploaded_file = st.file_uploader(
        "Upload raw count matrix",
        type=["tsv", "txt", "csv"],
        key="count_matrix_uploader_v1_4",
        help="Tab-delimited .tsv or .txt is recommended.",
    )

    if uploaded_file is None and st.session_state.get("counts_file_signature") is not None:
        clear_count_matrix_state()
        return

    if uploaded_file is not None:
        uploaded_signature = _uploaded_file_signature(uploaded_file)
        if st.session_state.get("counts_file_signature") != uploaded_signature:
            try:
                raw_counts_df = read_count_matrix_file(uploaded_file)
                st.session_state["raw_counts_df"] = raw_counts_df
                st.session_state["counts_file_name"] = uploaded_file.name
                st.session_state["counts_file_signature"] = uploaded_signature
                st.session_state["gene_id_column"] = str(raw_counts_df.columns[0]) if len(raw_counts_df.columns) else None
                st.session_state["sample_columns"] = get_sample_columns(raw_counts_df, st.session_state["gene_id_column"])
                st.session_state["sample_group_assignments"] = {
                    sample: "Unassigned" for sample in st.session_state["sample_columns"]
                }
                _reprocess_uploaded_counts()
            except Exception as exc:
                st.error(f"Could not read the uploaded count matrix: {exc}")

    raw_counts_df = st.session_state["raw_counts_df"]
    if raw_counts_df is None:
        return

    columns = [str(column) for column in raw_counts_df.columns]
    current_gene_column = st.session_state["gene_id_column"] if st.session_state["gene_id_column"] in columns else columns[0]
    selected_gene_column = st.selectbox(
        "Gene ID column",
        options=columns,
        index=columns.index(current_gene_column),
        key="gene_id_column_selector_v1_4",
    )
    if selected_gene_column != st.session_state["gene_id_column"]:
        st.session_state["gene_id_column"] = selected_gene_column
        st.session_state["sample_columns"] = get_sample_columns(raw_counts_df, selected_gene_column)
        st.session_state["sample_group_assignments"] = {
            sample: "Unassigned" for sample in st.session_state["sample_columns"]
        }
        _reprocess_uploaded_counts()

    sample_columns = st.session_state["sample_columns"]
    render_info_line("Uploaded file", st.session_state["counts_file_name"])
    st.write(
        "Detected shape: "
        f"`{raw_counts_df.shape[0]:,}` genes × `{len(sample_columns):,}` samples"
    )
    render_info_line("Detected gene ID column", st.session_state["gene_id_column"])
    st.write("Detected samples:")
    render_inline_badges(sample_columns)

    render_gene_symbol_conversion_section()
    render_duplicate_gene_section()


def update_group_names_safely(old_group_1: str, new_group_1: str, old_group_2: str, new_group_2: str) -> None:
    """Update assignments when group labels change."""
    updated = {}
    for sample, group in st.session_state["sample_group_assignments"].items():
        if group == old_group_1:
            updated[sample] = new_group_1
        elif group == old_group_2:
            updated[sample] = new_group_2
        else:
            updated[sample] = group
    st.session_state["sample_group_assignments"] = updated
    for sample, group in updated.items():
        st.session_state[f"group_assignment_{sample}"] = group


def render_assign_sample_groups_tab() -> None:
    """Render sample group assignment controls and readiness summary."""
    st.subheader("Assign Sample Groups")

    if st.session_state["raw_counts_df"] is None:
        st.info("Please upload a count matrix first.")
        return

    sample_columns = st.session_state["sample_columns"]
    if not sample_columns:
        st.warning("No sample columns were detected. Check the selected gene ID column.")
        return

    col1, col2 = st.columns(2)
    with col1:
        group_1_input = st.text_input("Group 1 name", value=st.session_state["group_1_name"])
    with col2:
        group_2_input = st.text_input("Group 2 name", value=st.session_state["group_2_name"])

    old_group_1 = st.session_state["group_1_name"]
    old_group_2 = st.session_state["group_2_name"]
    if group_1_input != old_group_1 or group_2_input != old_group_2:
        update_group_names_safely(old_group_1, group_1_input, old_group_2, group_2_input)
        st.session_state["group_1_name"] = group_1_input
        st.session_state["group_2_name"] = group_2_input
        reset_analysis_state()

    group_1 = st.session_state["group_1_name"]
    group_2 = st.session_state["group_2_name"]
    if group_1.strip() == "" or group_2.strip() == "":
        st.error("Both group names must be non-empty.")
        return
    if group_1 == group_2:
        st.error("Group 1 and Group 2 names must be different.")
        return

    action_col1, action_col2 = st.columns(2)
    with action_col1:
        if st.button("Assign first half to Group 1, second half to Group 2"):
            midpoint = len(sample_columns) // 2
            st.session_state["sample_group_assignments"] = {
                sample: group_1 if index < midpoint else group_2
                for index, sample in enumerate(sample_columns)
            }
            for sample, group in st.session_state["sample_group_assignments"].items():
                st.session_state[f"group_assignment_{sample}"] = group
            reset_analysis_state()
            st.rerun()
    with action_col2:
        if st.button("Clear all assignments"):
            st.session_state["sample_group_assignments"] = {sample: "Unassigned" for sample in sample_columns}
            for sample in sample_columns:
                st.session_state[f"group_assignment_{sample}"] = "Unassigned"
            reset_analysis_state()
            st.rerun()

    st.markdown("### Sample assignments")
    options = ["Unassigned", group_1, group_2]
    for sample in sample_columns:
        current_value = st.session_state["sample_group_assignments"].get(sample, "Unassigned")
        if current_value not in options:
            current_value = "Unassigned"
        widget_key = f"group_assignment_{sample}"
        if st.session_state.get(widget_key) not in options:
            st.session_state[widget_key] = current_value
        selected_group = st.selectbox(sample, options=options, index=options.index(current_value), key=widget_key)
        if selected_group != st.session_state["sample_group_assignments"].get(sample, "Unassigned"):
            st.session_state["sample_group_assignments"][sample] = selected_group
            reset_analysis_state()

    assignments = st.session_state["sample_group_assignments"]
    summary_df = pd.DataFrame(
        [
            {"Sample": sample, "Assigned group": assignments.get(sample, "Unassigned")}
            for sample in sample_columns
        ]
    )
    group_1_count = int((summary_df["Assigned group"] == group_1).sum())
    group_2_count = int((summary_df["Assigned group"] == group_2).sum())
    unassigned_count = int((summary_df["Assigned group"] == "Unassigned").sum())

    st.markdown("### Group assignment summary")
    st.dataframe(summary_df, use_container_width=True, hide_index=True)

    metric_col1, metric_col2, metric_col3 = st.columns(3)
    metric_col1.metric(group_1, group_1_count)
    metric_col2.metric(group_2, group_2_count)
    metric_col3.metric("Unassigned", unassigned_count)

    readiness = compute_basic_readiness()
    st.markdown("### Readiness")
    col_qc, col_deg = st.columns(2)
    col_qc.metric("Readiness for QC", readiness["qc"])
    col_deg.metric("Readiness for DEG", readiness["deg"])
    messages = readiness["messages"]
    if messages:
        st.warning("Review before continuing:")
        for message in messages:
            st.write(f"- {message}")
    else:
        st.success("Inputs are ready for QC and future DEG setup.")


def render_qc_overview_tab() -> None:
    """Render QC Overview data summaries without plots."""
    st.subheader("QC Overview")

    processed_counts_df = st.session_state["processed_counts_df"]
    sample_columns = st.session_state["sample_columns"]
    if processed_counts_df is None:
        st.info("Please upload and process a count matrix first.")
        return
    if not sample_columns:
        st.warning("No sample columns are available for QC.")
        return

    qc_summary = compute_qc_summary(processed_counts_df, sample_columns)
    sample_qc_df = qc_summary["sample_qc_df"]
    gene_qc_summary = qc_summary["gene_qc_summary"]

    st.markdown("### Dataset summary")
    render_info_line("Total processed genes", f"{gene_qc_summary['Total processed genes']:,}")
    render_info_line("Number of samples", f"{len(sample_columns):,}")
    render_info_line("Expressed genes", f"{gene_qc_summary['Expressed genes']:,}")
    render_info_line("All-zero genes", f"{gene_qc_summary['All-zero genes']:,}")
    render_info_line("Low-count genes", f"{gene_qc_summary['Low-count genes']:,}")
    render_info_line("Constant genes", f"{gene_qc_summary['Constant genes']:,}")

    st.markdown("### Sample-level QC")
    display_sample_qc = sample_qc_df.copy()
    display_sample_qc["Library size"] = display_sample_qc["Library size"].round(0).astype("int64")
    display_sample_qc["Zero fraction"] = display_sample_qc["Zero fraction"].round(4)
    display_sample_qc["Median count"] = display_sample_qc["Median count"].round(3)
    display_sample_qc["Mean count"] = display_sample_qc["Mean count"].round(3)
    st.dataframe(display_sample_qc, use_container_width=True, hide_index=True)

    st.markdown("### QC readiness")
    readiness = compute_basic_readiness()
    render_info_line("Readiness for QC", readiness["qc"])
    render_info_line("Readiness for DEG", readiness["deg"])
    if readiness["messages"]:
        st.warning("Review before continuing:")
        for message in readiness["messages"]:
            st.write(f"- {message}")
    else:
        st.success("Inputs are ready for QC and future DEG setup.")


def main() -> None:
    st.set_page_config(
        page_title="Bulk RNA-seq Explorer",
        page_icon="",
        layout="wide",
        initial_sidebar_state="expanded",
    )
    init_session_state()

    if st.session_state["gene_map_status"].get("message") == "Not loaded":
        gene_map_df, gene_map_status = load_local_gene_map()
        st.session_state["gene_map_df"] = gene_map_df
        st.session_state["gene_map_status"] = gene_map_status

    st.title("Bulk RNA-seq Explorer")
    st.caption("Python/Streamlit MVP migrated from validated browser prototype")

    render_sidebar()

    upload_tab, groups_tab, qc_tab = st.tabs(["Upload Count Matrix", "Assign Sample Groups", "QC Overview"])
    with upload_tab:
        render_upload_count_matrix_tab()
    with groups_tab:
        render_assign_sample_groups_tab()
    with qc_tab:
        render_qc_overview_tab()


if __name__ == "__main__":
    main()
