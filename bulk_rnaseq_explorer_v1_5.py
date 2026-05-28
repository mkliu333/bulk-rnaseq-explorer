"""
Bulk RNA-seq Explorer
Version: bulk_rnaseq_explorer_v1_5

Scope for v1.5:
- Clean Streamlit product UI for count-matrix upload and sample grouping.
- Detect whether the uploaded gene IDs are Ensembl IDs, gene symbols, mixed, or unclear.
- Convert mouse Ensembl IDs to gene symbols when a local mapping can be parsed.
- Merge duplicated processed gene symbols by summing raw counts.
- Produce a clean processed count matrix for future QC.
- Add Quality Control dataset summary and Plotly bar plots.
- Add multiselect-based QC grouping for stable sample/group QC views.

To reduce Streamlit toolbar/menu visibility, users may create `.streamlit/config.toml` with:

[client]
toolbarMode = "minimal"

Required for QC plots:
pip install plotly kaleido

Optional for faster gene-map cache:
pip install pyarrow

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
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st


APP_VERSION = "bulk_rnaseq_explorer_v1_5"

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
        "qc_grouping_sets": {},
        "active_qc_grouping_set": None,
        "current_qc_group_editor": {
            "grouping_set_name": "QC grouping 1",
            "groups": {
                "Group 1": [],
                "Group 2": [],
            },
        },
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
    st.session_state["qc_grouping_sets"] = {}
    st.session_state["active_qc_grouping_set"] = None
    st.session_state["current_qc_group_editor"] = {
        "grouping_set_name": "QC grouping 1",
        "groups": {"Group 1": [], "Group 2": []},
    }
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


def load_cached_gene_map() -> tuple[pd.DataFrame | None, dict[str, Any]]:
    """Load Python-friendly cached mapping when available."""
    cache_path = Path.cwd() / "assets" / "mouse_ensembl_to_symbol.parquet"
    status = {
        "detected": cache_path.exists(),
        "parsed": False,
        "path": str(cache_path),
        "n_mappings": 0,
        "message": "Cache not detected.",
    }
    if not cache_path.exists():
        return None, status
    try:
        mapping_df = pd.read_parquet(cache_path)
        if mapping_df.empty or not {"ensembl_id", "gene_symbol"}.issubset(mapping_df.columns):
            status["message"] = "Cache exists but does not contain expected mapping columns."
            return None, status
        if "ensembl_base" not in mapping_df.columns:
            mapping_df["ensembl_base"] = mapping_df["ensembl_id"].astype("string").str.replace(r"\.\d+$", "", regex=True)
        status.update({"parsed": True, "n_mappings": int(mapping_df.shape[0]), "message": "Loaded cached mapping."})
        return mapping_df, status
    except Exception as exc:
        status["message"] = f"Cache read failed: {exc}"
        return None, status


def save_gene_map_cache(gene_map_df: pd.DataFrame) -> dict[str, Any]:
    """Save mapping cache if parquet support is available."""
    cache_path = Path.cwd() / "assets" / "mouse_ensembl_to_symbol.parquet"
    status = {"path": str(cache_path), "saved": False, "message": ""}
    try:
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        gene_map_df.to_parquet(cache_path, index=False)
        status.update({"saved": True, "message": "Gene map cache saved."})
    except Exception as exc:
        status["message"] = f"Gene map cache not saved: {exc}"
    return status


def load_local_gene_map() -> tuple[pd.DataFrame | None, dict[str, Any]]:
    """Load local gene map, preferring parquet cache and falling back to source files."""
    cwd = Path.cwd()
    cached_df, cached_status = load_cached_gene_map()
    if cached_df is not None:
        return cached_df, cached_status

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
            cache_status = save_gene_map_cache(mapping_df)
            if not cache_status["saved"]:
                status["message"] = f"{status['message']} {cache_status['message']}"
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


def aggregate_sample_qc_by_group(
    sample_qc_df: pd.DataFrame,
    grouping_dict: dict[str, list[str]],
    metric_col: str,
    aggregation: str,
) -> tuple[pd.DataFrame, list[str]]:
    """Aggregate one sample-level QC metric by QC group."""
    rows = []
    missing_samples = []
    sample_set = set(sample_qc_df["Sample"].astype(str))
    for group_name, samples in grouping_dict.items():
        valid_samples = [sample for sample in samples if sample in sample_set]
        missing_samples.extend([sample for sample in samples if sample not in sample_set])
        values = sample_qc_df.loc[sample_qc_df["Sample"].isin(valid_samples), metric_col]
        if values.empty:
            continue
        if aggregation == "Median":
            value = float(values.median())
        elif aggregation == "Sum":
            value = float(values.sum())
        else:
            value = float(values.mean())
        rows.append({"Group": group_name, "Value": value, "N samples": len(valid_samples)})
    return pd.DataFrame(rows), sorted(set(missing_samples))


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
    st.sidebar.success("2. Quality Control")
    st.sidebar.caption("3. DEG Analysis - Coming soon / Locked")
    st.sidebar.caption("4. Visualization - Coming soon / Locked")
    st.sidebar.caption("5. Pathway Analysis - Coming soon / Locked")
    st.sidebar.caption("6. Export - Coming soon / Locked")


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
                st.session_state["current_qc_group_editor"] = {
                    "grouping_set_name": "QC grouping 1",
                    "groups": {"Group 1": [], "Group 2": []},
                }
                st.session_state["qc_grouping_sets"] = {}
                st.session_state["active_qc_grouping_set"] = None
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
        st.session_state["current_qc_group_editor"] = {
            "grouping_set_name": "QC grouping 1",
            "groups": {"Group 1": [], "Group 2": []},
        }
        st.session_state["qc_grouping_sets"] = {}
        st.session_state["active_qc_grouping_set"] = None
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


def make_qc_bar_plot(
    data_df: pd.DataFrame,
    x_col: str,
    y_col: str,
    title: str,
    y_axis_title: str,
    group_col: str | None = None,
    y_tick_format: str | None = None,
) -> go.Figure:
    """Create a clean Plotly QC bar plot."""
    marker_color = "#4C78A8" if group_col is None else "#59A14F"
    fig = go.Figure(
        data=[
            go.Bar(
                x=data_df[x_col],
                y=data_df[y_col],
                marker_color=marker_color,
                customdata=data_df[[group_col]].to_numpy() if group_col and group_col in data_df else None,
                hovertemplate=f"{x_col}: %{{x}}<br>{y_axis_title}: %{{y:,.4g}}<extra></extra>",
            )
        ]
    )
    fig.update_layout(
        title=title,
        height=520,
        plot_bgcolor="white",
        paper_bgcolor="white",
        margin=dict(l=70, r=30, t=70, b=120),
        xaxis=dict(tickangle=45, title="", showgrid=False),
        yaxis=dict(title=y_axis_title, showgrid=True, gridcolor="#E5E7EB", tickformat=y_tick_format),
        bargap=0.25,
        font=dict(size=13),
    )
    return fig


def render_qc_grouping_section(sample_columns: list[str]) -> None:
    """Render stable multiselect-based QC grouping editor.

    Drag-and-drop grouping can be added later via a custom Streamlit component;
    v1.5 uses multiselect-based grouping for stability and maintainability.
    """
    st.markdown("### Assign QC grouping")
    editor = st.session_state["current_qc_group_editor"]
    grouping_name = st.text_input("Grouping set name", value=editor["grouping_set_name"], key="qc_grouping_name_input")
    editor["grouping_set_name"] = grouping_name

    group_names = list(editor["groups"].keys())
    if st.button("Add group"):
        next_index = len(group_names) + 1
        editor["groups"][f"Group {next_index}"] = []
        st.rerun()

    updated_groups: dict[str, list[str]] = {}
    for index, old_group_name in enumerate(group_names):
        current_samples = editor["groups"].get(old_group_name, [])
        assigned_elsewhere = {
            sample
            for other_group, samples in editor["groups"].items()
            if other_group != old_group_name
            for sample in samples
        }
        options = [sample for sample in sample_columns if sample not in assigned_elsewhere or sample in current_samples]
        cols = st.columns([2, 5, 1])
        with cols[0]:
            new_group_name = st.text_input("Group name", value=old_group_name, key=f"qc_group_name_{index}")
        with cols[1]:
            selected_samples = st.multiselect(
                "Samples",
                options=options,
                default=[sample for sample in current_samples if sample in options],
                key=f"qc_group_samples_{index}",
            )
        with cols[2]:
            remove = st.button("Remove", key=f"qc_group_remove_{index}", disabled=len(group_names) <= 2)
        if remove and len(group_names) > 2:
            editor["groups"].pop(old_group_name, None)
            st.session_state["current_qc_group_editor"] = editor
            st.rerun()
        if not remove:
            clean_name = new_group_name.strip() or old_group_name
            updated_groups[clean_name] = selected_samples

    if len(updated_groups) < 2:
        updated_groups = editor["groups"]
        st.warning("At least two groups are required.")
    editor["groups"] = updated_groups
    st.session_state["current_qc_group_editor"] = editor

    assigned_samples = {sample for samples in editor["groups"].values() for sample in samples}
    unassigned_samples = [sample for sample in sample_columns if sample not in assigned_samples]
    st.write("Unassigned samples:")
    render_inline_badges(unassigned_samples)

    if st.button("Save QC grouping"):
        clean_name = editor["grouping_set_name"].strip() or "QC grouping 1"
        st.session_state["qc_grouping_sets"][clean_name] = {
            group: samples for group, samples in editor["groups"].items() if group.strip()
        }
        st.session_state["active_qc_grouping_set"] = clean_name
        st.success(f"Saved QC grouping: {clean_name}")

    grouping_names = list(st.session_state["qc_grouping_sets"].keys())
    if grouping_names:
        current_active = st.session_state["active_qc_grouping_set"]
        active_index = grouping_names.index(current_active) if current_active in grouping_names else 0
        active = st.selectbox("Active QC grouping set", grouping_names, index=active_index)
        st.session_state["active_qc_grouping_set"] = active
        if st.button("Delete active grouping set"):
            st.session_state["qc_grouping_sets"].pop(active, None)
            remaining = list(st.session_state["qc_grouping_sets"].keys())
            st.session_state["active_qc_grouping_set"] = remaining[0] if remaining else None
            st.rerun()


def render_qc_barplot_section(
    title: str,
    description: str,
    sample_qc_df: pd.DataFrame,
    metric_col: str,
    y_axis_title: str,
    default_aggregation: str,
    y_tick_format: str | None = None,
) -> None:
    """Render controls and Plotly bar plot for one QC metric."""
    st.markdown(f"### {title}")
    st.write(description)

    grouping_sets = st.session_state["qc_grouping_sets"]
    active_grouping = st.session_state["active_qc_grouping_set"]
    can_group = bool(active_grouping and active_grouping in grouping_sets)

    col1, col2, col3 = st.columns(3)
    with col1:
        plot_by_options = ["Sample name"] + (["QC group"] if can_group else [])
        plot_by = st.selectbox("Plot by", plot_by_options, key=f"{metric_col}_plot_by")
    with col2:
        if can_group:
            grouping_name = st.selectbox(
                "QC grouping set",
                list(grouping_sets.keys()),
                index=list(grouping_sets.keys()).index(active_grouping),
                key=f"{metric_col}_grouping_set",
            )
            st.session_state["active_qc_grouping_set"] = grouping_name
        else:
            st.selectbox("QC grouping set", ["No QC grouping selected"], disabled=True, key=f"{metric_col}_grouping_set_disabled")
            grouping_name = None
    with col3:
        aggregations = ["Mean", "Median", "Sum"]
        aggregation = st.selectbox(
            "Aggregation",
            aggregations,
            index=aggregations.index(default_aggregation),
            key=f"{metric_col}_aggregation",
        )

    if plot_by == "QC group" and grouping_name:
        plot_df, missing_samples = aggregate_sample_qc_by_group(
            sample_qc_df,
            grouping_sets[grouping_name],
            metric_col,
            aggregation,
        )
        if missing_samples:
            st.warning(f"Ignored samples not found in QC table: {', '.join(missing_samples[:8])}")
        if plot_df.empty:
            st.warning("No valid samples are assigned in the selected QC grouping set.")
            return
        fig = make_qc_bar_plot(plot_df, "Group", "Value", title, y_axis_title, group_col="Group", y_tick_format=y_tick_format)
    else:
        plot_df = sample_qc_df[["Sample", metric_col]].rename(columns={metric_col: "Value"})
        fig = make_qc_bar_plot(plot_df, "Sample", "Value", title, y_axis_title, y_tick_format=y_tick_format)

    st.plotly_chart(fig, use_container_width=True)


def render_quality_control_tab() -> None:
    """Render Quality Control summaries, grouping, and bar plots."""
    st.subheader("Quality Control")

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

    st.markdown("### Dataset summary")
    display_sample_qc = sample_qc_df.copy()
    display_sample_qc["Library size"] = display_sample_qc["Library size"].round(0).astype("int64")
    display_sample_qc["Zero fraction"] = display_sample_qc["Zero fraction"].round(4)
    display_sample_qc["Median count"] = display_sample_qc["Median count"].round(3)
    display_sample_qc["Mean count"] = display_sample_qc["Mean count"].round(3)
    st.dataframe(display_sample_qc, use_container_width=True, hide_index=True)

    render_qc_grouping_section(sample_columns)
    render_qc_barplot_section(
        "Library Size",
        "Total raw counts per sample or QC group.",
        sample_qc_df,
        "Library size",
        "Total raw counts",
        "Mean",
    )
    render_qc_barplot_section(
        "Detected Genes",
        "Genes with non-zero raw counts per sample or QC group.",
        sample_qc_df,
        "Detected genes",
        "Detected genes",
        "Mean",
    )
    render_qc_barplot_section(
        "Zero-count Fraction",
        "Fraction of genes with zero raw counts.",
        sample_qc_df,
        "Zero fraction",
        "Zero-count fraction",
        "Mean",
        y_tick_format=".0%",
    )


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

    upload_tab, qc_tab = st.tabs(["Upload Count Matrix", "Quality Control"])
    with upload_tab:
        render_upload_count_matrix_tab()
    with qc_tab:
        render_quality_control_tab()


if __name__ == "__main__":
    main()
