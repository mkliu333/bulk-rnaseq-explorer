"""
Bulk RNA-seq Explorer
Version: bulk_rnaseq_explorer_v1_9

Scope for v1.9:
- Clean Streamlit product UI for count-matrix upload and sample grouping.
- Detect whether the uploaded gene IDs are Ensembl IDs, gene symbols, mixed, or unclear.
- Convert mouse Ensembl IDs to gene symbols when a local mapping can be parsed.
- Merge duplicated processed gene symbols by summing raw counts.
- Produce a clean processed count matrix for future QC.
- Fix Quality Control barplot reset, color settings, advanced controls, and exports.
- Keep QC barplot state deterministic across Streamlit reruns.

To reduce Streamlit toolbar/menu visibility, users may create `.streamlit/config.toml` with:

[client]
toolbarMode = "minimal"

Required for QC plots:
pip install plotly kaleido

Optional for faster gene-map cache:
pip install pyarrow

Explicitly out of scope:
- DESeq2, Rscript, DEG analysis, normalization, PCA, sample correlation,
- heatmap, volcano plot,
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
import plotly.graph_objects as go
import streamlit as st


APP_VERSION = "bulk_rnaseq_explorer_v1_9"

DEFAULT_QC_COLORS = [
    "#355070", "#6d597a", "#b56576", "#e56b6f",
    "#4d908e", "#577590", "#bc6c25", "#8d99ae",
    "#2a9d8f", "#7f5539", "#4361ee", "#588157",
]

QC_PLOT_DEFINITIONS = {
    "library_size": {
        "title": "Library Size",
        "metric_col": "Library size",
        "description": "Total raw counts per sample or QC group.",
        "y_axis_title": "Total raw counts",
        "y_tick_format": None,
    },
    "detected_genes": {
        "title": "Detected Genes",
        "metric_col": "Detected genes",
        "description": "Genes with non-zero raw counts per sample or QC group.",
        "y_axis_title": "Detected genes",
        "y_tick_format": None,
    },
    "zero_fraction": {
        "title": "Zero-count Fraction",
        "metric_col": "Zero fraction",
        "description": "Fraction of genes with zero raw counts.",
        "y_axis_title": "Zero-count fraction",
        "y_tick_format": ".0%",
    },
}

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
            "groups": [
                {"id": "group_1", "name": "Group 1", "samples": []},
                {"id": "group_2", "name": "Group 2", "samples": []},
            ],
        },
        "qc_group_id_counter": 2,
        "qc_plot_settings": default_qc_plot_settings(),
        "qc_export_bytes": {},
        "qc_summary": None,
        "qc_summary_signature": None,
        "qc_plot_data_cache": {},
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
    st.session_state["app_version"] = APP_VERSION
    plot_settings = st.session_state.setdefault("qc_plot_settings", default_qc_plot_settings())
    for plot_id in QC_PLOT_DEFINITIONS:
        setting = plot_settings.setdefault(plot_id, get_default_qc_plot_setting(plot_id))
        for setting_key, setting_value in get_default_qc_plot_setting(plot_id).items():
            setting.setdefault(setting_key, setting_value)
    st.session_state.setdefault("qc_export_bytes", {})
    st.session_state.setdefault("qc_summary", None)
    st.session_state.setdefault("qc_summary_signature", None)
    st.session_state.setdefault("qc_plot_data_cache", {})
    st.session_state.setdefault("qc_group_id_counter", 2)


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
        "groups": [
            {"id": "group_1", "name": "Group 1", "samples": []},
            {"id": "group_2", "name": "Group 2", "samples": []},
        ],
    }
    st.session_state["qc_group_id_counter"] = 2
    st.session_state["qc_plot_settings"] = default_qc_plot_settings()
    st.session_state["qc_export_bytes"] = {}
    st.session_state["qc_summary"] = None
    st.session_state["qc_summary_signature"] = None
    st.session_state["qc_plot_data_cache"] = {}
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


def _qc_summary_signature() -> str:
    """Return a stable signature for the current processed count matrix and samples."""
    sample_columns = st.session_state.get("sample_columns", [])
    return json.dumps(
        {
            "counts": st.session_state.get("counts_file_signature"),
            "gene_id_column": st.session_state.get("gene_id_column"),
            "samples": sample_columns,
            "processed_shape": getattr(st.session_state.get("processed_counts_df"), "shape", None),
        },
        sort_keys=True,
        default=str,
    )


def compute_and_store_qc_summary() -> dict[str, Any] | None:
    """Compute QC summary once per processed matrix and store it in session state."""
    processed_counts_df = st.session_state.get("processed_counts_df")
    sample_columns = st.session_state.get("sample_columns", [])
    if processed_counts_df is None or not sample_columns:
        st.session_state["qc_summary"] = None
        st.session_state["qc_summary_signature"] = None
        return None
    signature = _qc_summary_signature()
    qc_summary = compute_qc_summary(processed_counts_df, sample_columns)
    st.session_state["qc_summary"] = qc_summary
    st.session_state["qc_summary_signature"] = signature
    return qc_summary


def get_cached_qc_summary() -> dict[str, Any] | None:
    """Return cached QC summary, recomputing only when the processed matrix changed."""
    signature = _qc_summary_signature()
    if st.session_state.get("qc_summary") is None or st.session_state.get("qc_summary_signature") != signature:
        return compute_and_store_qc_summary()
    return st.session_state["qc_summary"]


def aggregate_sample_qc_by_group(
    sample_qc_df: pd.DataFrame,
    grouping_dict: dict[str, list[str]],
    metric_col: str,
    aggregation: str,
) -> tuple[pd.DataFrame, pd.DataFrame, list[str]]:
    """Aggregate one sample-level QC metric by QC group and return overlay sample dots."""
    rows = []
    overlay_rows = []
    missing_samples = []
    sample_set = set(sample_qc_df["Sample"].astype(str))
    for group_name, samples in grouping_dict.items():
        valid_samples = [sample for sample in samples if sample in sample_set]
        missing_samples.extend([sample for sample in samples if sample not in sample_set])
        values = sample_qc_df.loc[sample_qc_df["Sample"].isin(valid_samples), metric_col]
        if values.empty:
            continue
        method = aggregation if aggregation in {"Mean", "Median", "Sum"} else "Mean"
        if aggregation == "Median":
            value = float(values.median())
        elif aggregation == "Sum":
            value = float(values.sum())
        else:
            value = float(values.mean())
        rows.append({"Label": group_name, "Value": value, "ColorKey": f"group:{group_name}", "N samples": len(valid_samples)})
        for sample in valid_samples:
            sample_value = sample_qc_df.loc[sample_qc_df["Sample"] == sample, metric_col].iloc[0]
            overlay_rows.append(
                {
                    "Group": group_name,
                    "Sample": sample,
                    "Value": float(sample_value),
                    "ColorKey": f"group:{group_name}",
                    "Aggregation": method,
                }
            )
    return pd.DataFrame(rows), pd.DataFrame(overlay_rows), sorted(set(missing_samples))


def default_qc_plot_settings() -> dict[str, dict[str, Any]]:
    """Return default settings for all QC bar plots."""
    return {plot_id: get_default_qc_plot_setting(plot_id) for plot_id in QC_PLOT_DEFINITIONS}


def get_default_qc_plot_setting(plot_id: str) -> dict[str, Any]:
    """Return default settings for one QC bar plot."""
    return {
        "plot_by": "Sample name",
        "grouping_set": None,
        "aggregation": "Mean",
        "width": 980,
        "height": 480,
        "x_axis_angle": 45,
        "axis_title_font_size": 14,
        "colors": {},
    }


def reset_qc_plot_setting(plot_id: str) -> None:
    """Reset one QC plot, including Streamlit widget state that can repopulate old values."""
    st.session_state["qc_plot_settings"][plot_id] = get_default_qc_plot_setting(plot_id)
    widget_suffixes = [
        "plot_by",
        "grouping_set",
        "grouping_set_disabled",
        "aggregation",
        "width",
        "height",
        "x_axis_angle",
        "axis_title_font_size",
    ]
    for suffix in widget_suffixes:
        st.session_state.pop(f"{plot_id}_{suffix}", None)
    for key in list(st.session_state.keys()):
        if str(key).startswith(f"{plot_id}_color_"):
            st.session_state.pop(key, None)
    clear_qc_export_cache(plot_id)
    clear_qc_plot_cache(plot_id)
    st.rerun()


def clear_qc_export_cache(plot_id: str | None = None) -> None:
    """Clear cached static Plotly export bytes."""
    if plot_id is None:
        st.session_state["qc_export_bytes"] = {}
        return
    for key in list(st.session_state.get("qc_export_bytes", {})):
        if plot_id in str(key):
            st.session_state["qc_export_bytes"].pop(key, None)


def clear_qc_plot_cache(plot_id: str | None = None) -> None:
    """Clear cached QC plot data."""
    if plot_id is None:
        st.session_state["qc_plot_data_cache"] = {}
        clear_qc_export_cache()
        return
    cache = st.session_state.setdefault("qc_plot_data_cache", {})
    for key in list(cache):
        if key.startswith(f"{plot_id}:"):
            cache.pop(key, None)


def default_qc_group_editor(sample_columns: list[str]) -> dict[str, Any]:
    """Build a default QC grouping draft."""
    return {
        "grouping_set_name": "QC grouping 1",
        "groups": [
            {"id": "group_1", "name": "Group 1", "samples": []},
            {"id": "group_2", "name": "Group 2", "samples": []},
        ],
    }


def normalize_qc_group_editor(editor: dict[str, Any] | None, sample_columns: list[str]) -> dict[str, Any]:
    """Keep the grouping draft compatible with current samples and stable group IDs."""
    if not isinstance(editor, dict):
        return default_qc_group_editor(sample_columns)
    sample_set = set(sample_columns)
    raw_groups = editor.get("groups")
    normalized_groups: list[dict[str, Any]] = []

    if isinstance(raw_groups, dict):
        iterable_groups = [
            {"id": f"group_{index}", "name": group_name, "samples": samples}
            for index, (group_name, samples) in enumerate(raw_groups.items(), start=1)
        ]
    elif isinstance(raw_groups, list):
        iterable_groups = raw_groups
    else:
        iterable_groups = []

    used_ids: set[str] = set()
    for index, group in enumerate(iterable_groups, start=1):
        if not isinstance(group, dict):
            continue
        group_id = str(group.get("id") or f"group_{index}").strip() or f"group_{index}"
        if group_id in used_ids:
            group_id = f"group_{index}"
        used_ids.add(group_id)
        clean_name = str(group.get("name") or f"Group {index}").strip() or f"Group {index}"
        samples = group.get("samples", [])
        unique_samples = []
        for sample in samples if isinstance(samples, list) else []:
            sample_text = str(sample)
            if sample_text in sample_set and sample_text not in unique_samples:
                unique_samples.append(sample_text)
        normalized_groups.append({"id": group_id, "name": clean_name, "samples": unique_samples})

    while len(normalized_groups) < 2:
        next_index = len(normalized_groups) + 1
        normalized_groups.append({"id": f"group_{next_index}", "name": f"Group {next_index}", "samples": []})

    max_numeric_id = 0
    for group in normalized_groups:
        match = re.match(r"^group_(\d+)$", str(group["id"]))
        if match:
            max_numeric_id = max(max_numeric_id, int(match.group(1)))
    st.session_state["qc_group_id_counter"] = max(int(st.session_state.get("qc_group_id_counter", 2)), max_numeric_id, 2)

    return {
        "grouping_set_name": str(editor.get("grouping_set_name") or "QC grouping 1"),
        "groups": normalized_groups,
    }


def get_editor_groups_as_dict(groups_or_editor: Any) -> dict[str, list[str]]:
    """Convert editor groups or saved grouping dict to {group_name: samples}."""
    if isinstance(groups_or_editor, dict) and isinstance(groups_or_editor.get("groups"), list):
        return {
            str(group.get("name", "")).strip(): list(group.get("samples", []))
            for group in groups_or_editor["groups"]
            if str(group.get("name", "")).strip()
        }
    if isinstance(groups_or_editor, list):
        return {
            str(group.get("name", "")).strip(): list(group.get("samples", []))
            for group in groups_or_editor
            if isinstance(group, dict) and str(group.get("name", "")).strip()
        }
    if isinstance(groups_or_editor, dict):
        return {
            str(group_name).strip(): list(samples)
            for group_name, samples in groups_or_editor.items()
            if str(group_name).strip() and isinstance(samples, list)
        }
    return {}


def get_unassigned_samples(groups_or_editor: Any, sample_columns: list[str]) -> list[str]:
    """Return samples not assigned to any QC group."""
    groups = get_editor_groups_as_dict(groups_or_editor)
    assigned = {sample for samples in groups.values() for sample in samples}
    return [sample for sample in sample_columns if sample not in assigned]


def get_assigned_samples_by_group(groups_or_editor: Any) -> dict[str, set[str]]:
    """Return selected samples for each QC group as sets."""
    groups = get_editor_groups_as_dict(groups_or_editor)
    return {group_name: set(samples) for group_name, samples in groups.items()}


def get_available_samples_for_group(
    group_id: str,
    editor: dict[str, Any],
    sample_columns: list[str],
) -> list[str]:
    """Return multiselect options with samples assigned elsewhere filtered out."""
    normalized_editor = normalize_qc_group_editor(editor, sample_columns)
    groups = normalized_editor["groups"]
    current_group = next((group for group in groups if group["id"] == group_id), None)
    current_group_samples = set(current_group.get("samples", []) if current_group else [])
    assigned_elsewhere = {
        sample
        for group in groups
        if group["id"] != group_id
        for sample in group.get("samples", [])
    }
    return [
        sample for sample in sample_columns
        if sample not in assigned_elsewhere or sample in current_group_samples
    ]


def get_duplicate_assigned_samples(groups_or_editor: Any) -> list[str]:
    """Return samples assigned to more than one QC group."""
    groups = get_editor_groups_as_dict(groups_or_editor)
    counts: dict[str, int] = {}
    for samples in groups.values():
        for sample in samples:
            counts[sample] = counts.get(sample, 0) + 1
    return sorted(sample for sample, count in counts.items() if count > 1)


def validate_qc_grouping(groups_or_editor: Any, sample_columns: list[str]) -> list[str]:
    """Validate a QC grouping before saving."""
    if isinstance(groups_or_editor, dict) and isinstance(groups_or_editor.get("groups"), list):
        group_names = [str(group.get("name", "")).strip() for group in groups_or_editor["groups"]]
    elif isinstance(groups_or_editor, list):
        group_names = [str(group.get("name", "")).strip() for group in groups_or_editor if isinstance(group, dict)]
    elif isinstance(groups_or_editor, dict):
        group_names = [str(group_name).strip() for group_name in groups_or_editor]
    else:
        group_names = []
    groups = get_editor_groups_as_dict(groups_or_editor)
    errors = []
    if len(group_names) < 2:
        errors.append("At least two groups are required.")
    blank_group_names = [name for name in group_names if not name.strip()]
    if blank_group_names:
        errors.append("Group names cannot be blank.")
    duplicate_group_names = len(group_names) != len({name for name in group_names if name})
    if duplicate_group_names:
        errors.append("Group names must be unique.")
    duplicates = get_duplicate_assigned_samples(groups)
    if duplicates:
        errors.append(f"Samples assigned to multiple groups: {', '.join(duplicates)}")
    sample_set = set(sample_columns)
    unknown = sorted({sample for samples in groups.values() for sample in samples if sample not in sample_set})
    if unknown:
        errors.append(f"Unknown samples: {', '.join(unknown)}")
    return errors


def save_qc_grouping_set(grouping_name: str, groups_or_editor: Any) -> str:
    """Persist a QC grouping set and make it active."""
    clean_name = grouping_name.strip() or "QC grouping 1"
    st.session_state["qc_grouping_sets"][clean_name] = get_editor_groups_as_dict(groups_or_editor)
    st.session_state["active_qc_grouping_set"] = clean_name
    clear_qc_plot_cache()
    reset_analysis_state()
    return clean_name


def _clear_qc_group_widget_keys(group_id: str | None = None) -> None:
    """Clear grouping editor widget keys so reset/remove actions do not leave stale UI state."""
    prefixes = ["qc_grouping_set_name_input"]
    if group_id is None:
        prefixes.extend(["qc_group_name_input_", "qc_group_samples_input_"])
    else:
        prefixes.extend([f"qc_group_name_input_{group_id}", f"qc_group_samples_input_{group_id}"])
    for key in list(st.session_state.keys()):
        if any(str(key).startswith(prefix) for prefix in prefixes):
            st.session_state.pop(key, None)


def update_qc_grouping_name_from_widget() -> None:
    """Sync the grouping set name widget into the editor draft."""
    editor = normalize_qc_group_editor(st.session_state.get("current_qc_group_editor"), st.session_state.get("sample_columns", []))
    editor["grouping_set_name"] = st.session_state.get("qc_grouping_set_name_input", "QC grouping 1")
    st.session_state["current_qc_group_editor"] = editor


def update_qc_group_name_from_widget(group_id: str) -> None:
    """Sync one group name widget into the editor draft."""
    editor = normalize_qc_group_editor(st.session_state.get("current_qc_group_editor"), st.session_state.get("sample_columns", []))
    widget_key = f"qc_group_name_input_{group_id}"
    for group in editor["groups"]:
        if group["id"] == group_id:
            group["name"] = str(st.session_state.get(widget_key, group["name"])).strip() or group["name"]
            break
    st.session_state["current_qc_group_editor"] = editor


def update_qc_group_samples_from_widget(group_id: str) -> None:
    """Sync one sample multiselect widget into the editor draft."""
    sample_columns = st.session_state.get("sample_columns", [])
    editor = normalize_qc_group_editor(st.session_state.get("current_qc_group_editor"), sample_columns)
    widget_key = f"qc_group_samples_input_{group_id}"
    selected_samples = [
        sample for sample in st.session_state.get(widget_key, [])
        if sample in sample_columns
    ]
    for group in editor["groups"]:
        if group["id"] == group_id:
            group["samples"] = selected_samples
            break
    st.session_state["current_qc_group_editor"] = editor


def add_qc_group_to_editor(sample_columns: list[str]) -> None:
    """Add one empty group to the QC grouping draft."""
    editor = normalize_qc_group_editor(st.session_state.get("current_qc_group_editor"), sample_columns)
    st.session_state["qc_group_id_counter"] = int(st.session_state.get("qc_group_id_counter", 2)) + 1
    group_id = f"group_{st.session_state['qc_group_id_counter']}"
    group_number = len(editor["groups"]) + 1
    editor["groups"].append({"id": group_id, "name": f"Group {group_number}", "samples": []})
    st.session_state["current_qc_group_editor"] = editor


def remove_qc_group_from_editor(group_id: str, sample_columns: list[str]) -> None:
    """Remove one group from the QC grouping draft."""
    editor = normalize_qc_group_editor(st.session_state.get("current_qc_group_editor"), sample_columns)
    if len(editor["groups"]) > 2:
        editor["groups"] = [group for group in editor["groups"] if group["id"] != group_id]
        _clear_qc_group_widget_keys(group_id)
    st.session_state["current_qc_group_editor"] = normalize_qc_group_editor(editor, sample_columns)


def clear_qc_group_editor(sample_columns: list[str]) -> None:
    """Reset only the unsaved QC grouping editor draft."""
    st.session_state["current_qc_group_editor"] = default_qc_group_editor(sample_columns)
    st.session_state["qc_group_id_counter"] = 2
    _clear_qc_group_widget_keys()


def get_qc_color(plot_id: str, entity_key: str, index: int) -> str:
    """Get a stable color for one sample or group in one plot."""
    settings = st.session_state["qc_plot_settings"].setdefault(plot_id, get_default_qc_plot_setting(plot_id))
    colors = settings.setdefault("colors", {})
    if entity_key not in colors:
        colors[entity_key] = DEFAULT_QC_COLORS[index % len(DEFAULT_QC_COLORS)]
    return colors[entity_key]


def set_qc_color(plot_id: str, entity_key: str, color: str) -> None:
    """Store a custom QC color."""
    settings = st.session_state["qc_plot_settings"].setdefault(plot_id, get_default_qc_plot_setting(plot_id))
    settings.setdefault("colors", {})[entity_key] = color


def render_qc_color_settings(plot_id: str, entity_keys: list[str], labels: list[str]) -> None:
    """Render color pickers for the current QC plot basis."""
    if not entity_keys:
        st.caption("Color settings are available after plot entities are available.")
        return
    st.markdown(
        """
        <style>
        div[data-testid="stColorPicker"] {
            margin-bottom: 0.25rem;
        }
        </style>
        """,
        unsafe_allow_html=True,
    )
    n_cols = min(8, max(1, len(entity_keys)))
    columns = st.columns(n_cols)
    for index, (entity_key, label) in enumerate(zip(entity_keys, labels)):
        with columns[index % n_cols]:
            color = st.color_picker(
                label,
                value=get_qc_color(plot_id, entity_key, index),
                key=f"{plot_id}_color_{entity_key}",
            )
            set_qc_color(plot_id, entity_key, color)


def prepare_qc_plot_data(
    sample_qc_df: pd.DataFrame,
    plot_id: str,
    metric_col: str,
    settings: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame | None, list[str]]:
    """Prepare sample-level or group-level QC plot data."""
    if settings.get("plot_by") != "QC assignment group":
        plot_df = sample_qc_df[["Sample", metric_col]].rename(columns={"Sample": "Label", metric_col: "Value"})
        plot_df["ColorKey"] = plot_df["Label"].map(lambda sample: f"sample:{sample}")
        return plot_df, None, []

    grouping_name = settings.get("grouping_set")
    grouping_sets = st.session_state.get("qc_grouping_sets", {})
    grouping_dict = grouping_sets.get(grouping_name)
    if not grouping_dict:
        return pd.DataFrame(), pd.DataFrame(), []
    aggregation = settings.get("aggregation", "Mean")
    if plot_id == "zero_fraction" and aggregation == "Sum":
        aggregation = "Mean"
        settings["aggregation"] = "Mean"
    return aggregate_sample_qc_by_group(sample_qc_df, grouping_dict, metric_col, aggregation)


def make_qc_plot_cache_key(
    plot_id: str,
    metric_col: str,
    settings: dict[str, Any],
) -> str:
    """Build a lightweight cache key for QC plot data, excluding visual-only settings."""
    grouping_set = settings.get("grouping_set") if settings.get("plot_by") == "QC assignment group" else None
    grouping_signature = ""
    if grouping_set:
        grouping_signature = json.dumps(
            st.session_state.get("qc_grouping_sets", {}).get(grouping_set, {}),
            sort_keys=True,
        )
    payload = {
        "plot_id": plot_id,
        "metric_col": metric_col,
        "plot_by": settings.get("plot_by", "Sample name"),
        "grouping_set": grouping_set,
        "aggregation": settings.get("aggregation", "Mean"),
        "qc_summary_signature": st.session_state.get("qc_summary_signature"),
        "grouping_signature": grouping_signature,
    }
    return f"{plot_id}:{hashlib.sha256(json.dumps(payload, sort_keys=True).encode('utf-8')).hexdigest()}"


def get_cached_qc_plot_data(
    sample_qc_df: pd.DataFrame,
    plot_id: str,
    metric_col: str,
    settings: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame | None, list[str], list[str], list[str]]:
    """Return cached QC plot data; visual settings and colors do not invalidate it."""
    cache_key = make_qc_plot_cache_key(plot_id, metric_col, settings)
    cache = st.session_state.setdefault("qc_plot_data_cache", {})
    if cache_key not in cache:
        plot_df, overlay_df, missing_samples = prepare_qc_plot_data(sample_qc_df, plot_id, metric_col, settings)
        cache[cache_key] = {
            "plot_df": plot_df,
            "overlay_df": overlay_df,
            "missing_samples": missing_samples,
            "entity_keys": list(plot_df["ColorKey"]) if "ColorKey" in plot_df else [],
            "entity_labels": list(plot_df["Label"]) if "Label" in plot_df else [],
        }
    cached = cache[cache_key]
    return (
        cached["plot_df"],
        cached["overlay_df"],
        cached["missing_samples"],
        cached["entity_keys"],
        cached["entity_labels"],
    )


def _hex_to_rgba(hex_color: str, alpha: float) -> str:
    """Convert #RRGGBB to an rgba color string."""
    color = hex_color.lstrip("#")
    if len(color) != 6:
        return f"rgba(53, 80, 112, {alpha})"
    red, green, blue = (int(color[i : i + 2], 16) for i in (0, 2, 4))
    return f"rgba({red}, {green}, {blue}, {alpha})"


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
    compute_and_store_qc_summary()
    clear_qc_plot_cache()
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
                    "groups": [
                        {"id": "group_1", "name": "Group 1", "samples": []},
                        {"id": "group_2", "name": "Group 2", "samples": []},
                    ],
                }
                st.session_state["qc_group_id_counter"] = 2
                st.session_state["qc_grouping_sets"] = {}
                st.session_state["active_qc_grouping_set"] = None
                st.session_state["qc_plot_settings"] = default_qc_plot_settings()
                st.session_state["qc_export_bytes"] = {}
                st.session_state["qc_plot_data_cache"] = {}
                st.session_state["qc_summary"] = None
                st.session_state["qc_summary_signature"] = None
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
            "groups": [
                {"id": "group_1", "name": "Group 1", "samples": []},
                {"id": "group_2", "name": "Group 2", "samples": []},
            ],
        }
        st.session_state["qc_group_id_counter"] = 2
        st.session_state["qc_grouping_sets"] = {}
        st.session_state["active_qc_grouping_set"] = None
        st.session_state["qc_plot_settings"] = default_qc_plot_settings()
        st.session_state["qc_export_bytes"] = {}
        st.session_state["qc_plot_data_cache"] = {}
        st.session_state["qc_summary"] = None
        st.session_state["qc_summary_signature"] = None
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
    plot_df: pd.DataFrame,
    plot_id: str,
    title: str,
    x_col: str,
    y_col: str,
    y_axis_title: str,
    width: int,
    height: int,
    x_axis_angle: int,
    colors: dict[str, str],
    axis_title_font_size: int = 14,
    y_tick_format: str | None = None,
    overlay_df: pd.DataFrame | None = None,
) -> go.Figure:
    """Create a clean Plotly QC bar plot."""
    bar_colors = [
        colors.get(color_key, get_qc_color(plot_id, str(color_key), index))
        for index, color_key in enumerate(plot_df["ColorKey"])
    ]
    customdata_cols = ["N samples"] if "N samples" in plot_df.columns else []
    customdata = plot_df[customdata_cols].to_numpy() if customdata_cols else None
    hovertemplate = (
        "Group: %{x}<br>Value: %{y:,.4g}<br>N samples: %{customdata[0]}<extra></extra>"
        if customdata_cols
        else "Sample: %{x}<br>Value: %{y:,.4g}<extra></extra>"
    )
    fig = go.Figure()
    fig.add_trace(
        go.Bar(
            x=plot_df[x_col],
            y=plot_df[y_col],
            marker_color=bar_colors,
            customdata=customdata,
            hovertemplate=hovertemplate,
            name=title,
        )
    )
    if overlay_df is not None and not overlay_df.empty:
        dot_colors = [
            _hex_to_rgba(colors.get(color_key, get_qc_color(plot_id, str(color_key), index)), 0.82)
            for index, color_key in enumerate(overlay_df["ColorKey"])
        ]
        fig.add_trace(
            go.Scatter(
                x=overlay_df["Group"],
                y=overlay_df["Value"],
                mode="markers",
                marker=dict(size=7, color=dot_colors, line=dict(width=1, color="rgba(17, 24, 39, 0.55)")),
                customdata=overlay_df[["Sample", "Group", "Aggregation"]].to_numpy(),
                hovertemplate=(
                    "Sample: %{customdata[0]}<br>"
                    "Group: %{customdata[1]}<br>"
                    "Value: %{y:,.4g}<br>"
                    "Aggregation: %{customdata[2]}<extra></extra>"
                ),
                name="Samples",
            )
        )
    fig.update_layout(
        title=title,
        width=width,
        height=height,
        plot_bgcolor="white",
        paper_bgcolor="white",
        margin=dict(l=70, r=30, t=70, b=110),
        xaxis=dict(tickangle=x_axis_angle, title="", showgrid=False, title_font=dict(size=axis_title_font_size)),
        yaxis=dict(
            title=y_axis_title,
            showgrid=True,
            gridcolor="#E5E7EB",
            tickformat=y_tick_format,
            title_font=dict(size=axis_title_font_size),
        ),
        bargap=0.25,
        font=dict(size=13),
        showlegend=False,
    )
    return fig


# Streamlit reruns the script on widget changes. v1.8 minimizes expensive
# recomputation and uses immediate session-state guards, but true instant
# drag/drop or fully local UI updates would require a custom Streamlit
# component or a React frontend.
def render_qc_group_editor_form(sample_columns: list[str]) -> None:
    """Render the immediate-update QC grouping draft editor."""
    editor = normalize_qc_group_editor(st.session_state.get("current_qc_group_editor"), sample_columns)
    if "qc_grouping_set_name_input" in st.session_state:
        editor["grouping_set_name"] = st.session_state["qc_grouping_set_name_input"]
    for group in editor["groups"]:
        name_key = f"qc_group_name_input_{group['id']}"
        samples_key = f"qc_group_samples_input_{group['id']}"
        if name_key in st.session_state:
            group["name"] = str(st.session_state[name_key]).strip() or group["name"]
        if samples_key in st.session_state:
            group["samples"] = [
                sample for sample in st.session_state[samples_key]
                if sample in sample_columns
            ]
    editor = normalize_qc_group_editor(editor, sample_columns)
    st.session_state["current_qc_group_editor"] = editor

    with st.container(border=True):
        if "qc_grouping_set_name_input" not in st.session_state:
            st.session_state["qc_grouping_set_name_input"] = editor["grouping_set_name"]
        st.text_input(
            "Grouping set name",
            key="qc_grouping_set_name_input",
            on_change=update_qc_grouping_name_from_widget,
        )

        for group in editor["groups"]:
            group_id = group["id"]
            name_key = f"qc_group_name_input_{group_id}"
            samples_key = f"qc_group_samples_input_{group_id}"
            options = get_available_samples_for_group(group_id, editor, sample_columns)
            st.session_state[name_key] = group["name"]
            st.session_state[samples_key] = [sample for sample in group.get("samples", []) if sample in options]

            cols = st.columns([2.0, 5.0, 1.35])
            with cols[0]:
                st.text_input(
                    "Group name",
                    key=name_key,
                    on_change=update_qc_group_name_from_widget,
                    args=(group_id,),
                )
            with cols[1]:
                st.multiselect(
                    "Samples",
                    options=options,
                    key=samples_key,
                    on_change=update_qc_group_samples_from_widget,
                    args=(group_id,),
                )
            with cols[2]:
                st.write("")
                st.button(
                    "Remove group",
                    key=f"qc_group_remove_{group_id}",
                    disabled=len(editor["groups"]) <= 2,
                    on_click=remove_qc_group_from_editor,
                    args=(group_id, sample_columns),
                )

        action_cols = st.columns([1.25, 1.35, 3.4, 1.05])
        with action_cols[0]:
            save_grouping = st.button("Save QC grouping", key="qc_group_save_button", type="primary")
        with action_cols[1]:
            st.button(
                "Clear grouping info",
                key="qc_group_clear_button",
                on_click=clear_qc_group_editor,
                args=(sample_columns,),
            )
        with action_cols[3]:
            st.button(
                "Add group",
                key="qc_group_add_button",
                on_click=add_qc_group_to_editor,
                args=(sample_columns,),
            )

    if save_grouping:
        update_qc_grouping_name_from_widget()
        current_editor = normalize_qc_group_editor(st.session_state["current_qc_group_editor"], sample_columns)
        for group in current_editor["groups"]:
            update_qc_group_name_from_widget(group["id"])
            update_qc_group_samples_from_widget(group["id"])
        current_editor = normalize_qc_group_editor(st.session_state["current_qc_group_editor"], sample_columns)
        errors = validate_qc_grouping(current_editor, sample_columns)
        if errors:
            for error in errors:
                st.error(error)
        else:
            saved_name = save_qc_grouping_set(current_editor["grouping_set_name"], current_editor)
            st.success(f"Saved QC grouping: {saved_name}")

    draft_groups = st.session_state["current_qc_group_editor"]["groups"]
    st.write("Unassigned samples:")
    render_inline_badges(get_unassigned_samples(draft_groups, sample_columns))
    duplicates = get_duplicate_assigned_samples(draft_groups)
    if duplicates:
        st.warning(f"Duplicate assignment in current draft: {', '.join(duplicates)}")


def render_qc_grouping_section(sample_columns: list[str]) -> None:
    """Render the immediate-update QC grouping editor and saved-set controls."""
    st.markdown("### Assign QC grouping")
    render_qc_group_editor_form(sample_columns)

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
            clear_qc_plot_cache()
            st.rerun()


def render_qc_barplot_section(
    plot_id: str,
    sample_qc_df: pd.DataFrame,
) -> None:
    """Render controls and Plotly bar plot for one QC metric."""
    definition = QC_PLOT_DEFINITIONS[plot_id]
    title = definition["title"]
    description = definition["description"]
    metric_col = definition["metric_col"]
    y_axis_title = definition["y_axis_title"]
    y_tick_format = definition["y_tick_format"]
    settings_store = st.session_state.setdefault("qc_plot_settings", default_qc_plot_settings())
    settings = settings_store.setdefault(plot_id, get_default_qc_plot_setting(plot_id))

    st.markdown(f"### {title}")
    st.write(description)

    grouping_sets = st.session_state["qc_grouping_sets"]
    grouping_names = list(grouping_sets.keys())
    if settings.get("grouping_set") not in grouping_sets:
        settings["grouping_set"] = st.session_state.get("active_qc_grouping_set") if st.session_state.get("active_qc_grouping_set") in grouping_sets else None
    if settings.get("plot_by") not in {"Sample name", "QC assignment group"}:
        settings["plot_by"] = "Sample name"
    settings.setdefault("axis_title_font_size", 14)
    if plot_id == "zero_fraction" and settings.get("aggregation") == "Sum":
        settings["aggregation"] = "Mean"
        st.session_state.pop(f"{plot_id}_aggregation", None)

    control_cols = st.columns([1.1, 1.45, 0.95, 0.55, 0.75, 0.75])
    with control_cols[0]:
        settings["plot_by"] = st.selectbox(
            "Plot by",
            ["Sample name", "QC assignment group"],
            index=["Sample name", "QC assignment group"].index(settings.get("plot_by", "Sample name")),
            key=f"{plot_id}_plot_by",
        )
    plot_by_group = settings["plot_by"] == "QC assignment group"
    effective_group_mode = plot_by_group and bool(grouping_names)
    with control_cols[1]:
        if effective_group_mode:
            selected_grouping = settings.get("grouping_set") if settings.get("grouping_set") in grouping_names else grouping_names[0]
            settings["grouping_set"] = st.selectbox(
                "QC assignment set",
                grouping_names,
                index=grouping_names.index(selected_grouping),
                key=f"{plot_id}_grouping_set",
            )
        else:
            st.selectbox(
                "QC assignment set",
                ["No QC assignment selected"],
                disabled=True,
                key=f"{plot_id}_grouping_set_disabled",
            )
            if not effective_group_mode:
                settings["grouping_set"] = None
    with control_cols[2]:
        aggregation_options = ["Mean", "Median"] if plot_id == "zero_fraction" else ["Mean", "Median", "Sum"]
        settings["aggregation"] = st.selectbox(
            "Aggregation",
            aggregation_options,
            index=aggregation_options.index(settings.get("aggregation", "Mean")),
            disabled=not effective_group_mode,
            key=f"{plot_id}_aggregation",
        )
    with control_cols[3]:
        st.write("")
        if st.button("Reset", key=f"{plot_id}_reset"):
            reset_qc_plot_setting(plot_id)

    with st.expander("Advanced settings"):
        adv_cols = st.columns(4)
        with adv_cols[0]:
            settings["width"] = st.slider("Plot width", 640, 1440, int(settings.get("width", 980)), 20, key=f"{plot_id}_width")
        with adv_cols[1]:
            settings["height"] = st.slider("Plot height", 360, 900, int(settings.get("height", 480)), 20, key=f"{plot_id}_height")
        with adv_cols[2]:
            angle_options = [0, 30, 45, 60, 90]
            settings["x_axis_angle"] = st.selectbox(
                "X-axis angle",
                angle_options,
                index=angle_options.index(int(settings.get("x_axis_angle", 45))),
                key=f"{plot_id}_x_axis_angle",
            )
        with adv_cols[3]:
            settings["axis_title_font_size"] = st.slider(
                "Axis title font size",
                10,
                28,
                int(settings.get("axis_title_font_size", 14)),
                1,
                key=f"{plot_id}_axis_title_font_size",
            )

    if plot_by_group and not grouping_names:
        st.warning("Please create and save a QC grouping set first.")
        settings["plot_by"] = "Sample name"

    plot_df, overlay_df, missing_samples, entity_keys, entity_labels = get_cached_qc_plot_data(
        sample_qc_df,
        plot_id,
        metric_col,
        settings,
    )
    if missing_samples:
        st.warning(f"Ignored samples not found in QC table: {', '.join(missing_samples[:8])}")
    if plot_df.empty:
        st.warning("No valid samples are assigned in the selected QC grouping set.")
        return

    with st.expander("Color settings"):
        render_qc_color_settings(plot_id, entity_keys, entity_labels)

    color_map = {
        entity_key: get_qc_color(plot_id, entity_key, index)
        for index, entity_key in enumerate(entity_keys)
    }
    fig = make_qc_bar_plot(
        plot_df,
        plot_id,
        title,
        "Label",
        "Value",
        y_axis_title,
        int(settings["width"]),
        int(settings["height"]),
        int(settings["x_axis_angle"]),
        color_map,
        axis_title_font_size=int(settings.get("axis_title_font_size", 14)),
        y_tick_format=y_tick_format,
        overlay_df=overlay_df,
    )

    export_signature = hashlib.sha256(
        json.dumps(
            {
                "plot_by": settings.get("plot_by"),
                "grouping_set": settings.get("grouping_set"),
                "aggregation": settings.get("aggregation"),
                "width": settings.get("width"),
                "height": settings.get("height"),
                "x_axis_angle": settings.get("x_axis_angle"),
                "axis_title_font_size": settings.get("axis_title_font_size"),
                "colors": settings.get("colors", {}),
                "qc_summary_signature": st.session_state.get("qc_summary_signature"),
            },
            sort_keys=True,
        ).encode("utf-8")
    ).hexdigest()[:8]
    filename_base = f"{plot_id}_{APP_VERSION}_{export_signature}"
    render_qc_export_buttons(fig, plot_id, filename_base, control_cols[4], control_cols[5])
    st.plotly_chart(fig, use_container_width=False, key=f"{plot_id}_plotly_chart")


def get_plotly_export_bytes(fig: go.Figure, format: str, cache_key: str) -> bytes | None:
    """Return cached Plotly static image bytes, generating them only for the current visual state."""
    export_cache = st.session_state.setdefault("qc_export_bytes", {})
    if cache_key in export_cache:
        return export_cache[cache_key]
    try:
        kwargs: dict[str, Any] = {"format": format}
        if format == "png":
            kwargs["scale"] = 3
        export_bytes = fig.to_image(**kwargs)
        if isinstance(export_bytes, str):
            export_bytes = export_bytes.encode("utf-8")
        export_cache[cache_key] = export_bytes
        return export_bytes
    except Exception:
        export_cache[cache_key] = None
        return None


def render_qc_export_buttons(
    fig: go.Figure,
    plot_id: str,
    filename_base: str,
    png_column,
    svg_column,
) -> None:
    """Render direct PNG/SVG download buttons for one QC plot."""
    png_key = f"{plot_id}:{filename_base}:png"
    svg_key = f"{plot_id}:{filename_base}:svg"
    png_bytes = get_plotly_export_bytes(fig, "png", png_key)
    svg_bytes = get_plotly_export_bytes(fig, "svg", svg_key)
    with png_column:
        st.write("")
        st.download_button(
            "Download PNG",
            data=png_bytes or b"",
            file_name=f"{filename_base}.png",
            mime="image/png",
            disabled=png_bytes is None,
            key=f"{plot_id}_download_png",
        )
    with svg_column:
        st.write("")
        st.download_button(
            "Download SVG",
            data=svg_bytes or b"",
            file_name=f"{filename_base}.svg",
            mime="image/svg+xml",
            disabled=svg_bytes is None,
            key=f"{plot_id}_download_svg",
        )
    if png_bytes is None or svg_bytes is None:
        st.warning("Static image export requires kaleido. Install with: pip install kaleido")


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

    qc_summary = get_cached_qc_summary()
    if qc_summary is None:
        st.warning("QC summary is unavailable for the current count matrix.")
        return
    sample_qc_df = qc_summary["sample_qc_df"]

    st.markdown("### Dataset summary")
    display_sample_qc = sample_qc_df.copy()
    display_sample_qc["Library size"] = display_sample_qc["Library size"].round(0).astype("int64")
    display_sample_qc["Zero fraction"] = display_sample_qc["Zero fraction"].round(4)
    display_sample_qc["Mean count"] = display_sample_qc["Mean count"].round(3)
    display_sample_qc = display_sample_qc[
        ["Sample", "Library size", "Detected genes", "Zero-count genes", "Zero fraction", "Mean count"]
    ]
    st.dataframe(display_sample_qc, use_container_width=True, hide_index=True)

    render_qc_grouping_section(sample_columns)
    render_qc_barplot_section("library_size", sample_qc_df)
    render_qc_barplot_section("detected_genes", sample_qc_df)
    render_qc_barplot_section("zero_fraction", sample_qc_df)


def main() -> None:
    st.set_page_config(
        page_title="Bulk RNA-seq Explorer",
        page_icon="",
        layout="wide",
        initial_sidebar_state="expanded",
    )
    init_session_state()
    st.markdown(
        """
        <style>
        div[data-testid="stButton"] button,
        div[data-testid="stDownloadButton"] button {
            white-space: nowrap;
        }
        </style>
        """,
        unsafe_allow_html=True,
    )

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
