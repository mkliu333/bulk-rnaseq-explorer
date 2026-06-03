"""
Bulk RNA-seq Explorer
Version: bulk_rnaseq_explorer_v2_3

Scope for v2.3:
- Clean Streamlit product UI for count-matrix upload and sample grouping.
- Detect whether the uploaded gene IDs are Ensembl IDs, gene symbols, mixed, or unclear.
- Convert mouse Ensembl IDs to gene symbols when a local mapping can be parsed.
- Merge duplicated processed gene symbols by summing raw counts.
- Produce a clean processed count matrix for future QC.
- Fix Quality Control grouping editor save/reset behavior and barplot reset state.
- Use placeholder-based QC grouping inputs and nonce-based QC plot widgets.
- Polish Quality Control button responsiveness and barplot axis label controls.
- Tighten Quality Control action-button rows on wide and narrow screens.
- Use adaptive QC button widths without truncating labels.
- Refactor Quality Control action rows into compact responsive clusters.
- Align Quality Control buttons to their related input grids with small gaps.
- Automatically run Normalization when processed count input changes.
- Keep QC barplot rendering fast by lazily preparing PNG/SVG export bytes.
- Refine Normalized matrix search, pagination, and CSV controls.
- Add PCA and Sample Correlation QC plots.
- Require normalized expression matrices for PCA and Sample Correlation.
- Refine Sample Correlation heatmap annotation layout.

To reduce Streamlit toolbar/menu visibility, users may create `.streamlit/config.toml` with:

[client]
toolbarMode = "minimal"

Required for QC plots:
pip install plotly kaleido

Optional for faster gene-map cache:
pip install pyarrow

Explicitly out of scope:
- DEG analysis, DEG heatmap, volcano plot,
  GSEA, ORA, pathway enrichment, cloud storage, login, Duke DCC, SLURM.
"""

from __future__ import annotations

import csv
import hashlib
import json
import re
import subprocess
import tempfile
from datetime import datetime
from io import StringIO
from pathlib import Path
from typing import Any

try:
    import numpy as np
except ImportError:  # pragma: no cover - surfaced in the Normalization UI.
    np = None

import pandas as pd
import plotly.graph_objects as go
import streamlit as st


APP_VERSION = "bulk_rnaseq_explorer_v2_3"

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

NORMALIZATION_OUTPUTS = {
    "Raw counts": "raw_counts.csv",
    "CPM": "cpm.csv",
    "log2(CPM + 1)": "log2_cpm_plus1.csv",
    "DESeq2 normalized counts": "deseq2_normalized_counts.csv",
    "DESeq2 VST": "deseq2_vst.csv",
    "edgeR TMM CPM": "edger_tmm_cpm.csv",
    "edgeR TMM logCPM": "edger_tmm_logcpm.csv",
}

NORMALIZATION_FACTOR_FILES = {
    "DESeq2 size factors": "deseq2_size_factors.csv",
    "edgeR TMM normalization factors": "edger_tmm_norm_factors.csv",
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
            "grouping_set_name": "",
            "groups": [
                {"id": "group_1", "name": "", "samples": []},
                {"id": "group_2", "name": "", "samples": []},
            ],
        },
        "qc_group_id_counter": 2,
        "qc_plot_settings": default_qc_plot_settings(),
        "qc_plot_reset_nonce": {},
        "qc_export_bytes": {},
        "qc_summary": None,
        "qc_summary_signature": None,
        "qc_plot_data_cache": {},
        "qc_active_view": "QC Summary & Grouping",
        "qc_pca_settings": default_qc_pca_settings(),
        "qc_pca_reset_nonce": 0,
        "qc_pca_cache": {},
        "qc_corr_settings": default_qc_corr_settings(),
        "qc_corr_reset_nonce": 0,
        "qc_corr_cache": {},
        "qc_expression_matrix_cache": {},
        "normalization_results": None,
        "normalization_output_dir": None,
        "normalization_report": None,
        "normalization_tables": {},
        "normalization_run_status": None,
        "normalization_input_signature": None,
        "normalization_selected_matrix": "Raw counts",
        "normalization_table_page": 1,
        "normalization_table_rows_per_page": 25,
        "normalization_table_search": "",
        "normalization_selected_gene": "",
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
        sync_qc_axis_label_font_size(setting)
    plot_reset_nonce = st.session_state.setdefault("qc_plot_reset_nonce", {})
    for plot_id in QC_PLOT_DEFINITIONS:
        plot_reset_nonce.setdefault(plot_id, 0)
    st.session_state.setdefault("qc_export_bytes", {})
    st.session_state.setdefault("qc_summary", None)
    st.session_state.setdefault("qc_summary_signature", None)
    st.session_state.setdefault("qc_plot_data_cache", {})
    st.session_state.setdefault("qc_active_view", "QC Summary & Grouping")
    if st.session_state.get("qc_active_view") == "QC Summary":
        st.session_state["qc_active_view"] = "QC Summary & Grouping"
    st.session_state.setdefault("qc_pca_settings", default_qc_pca_settings())
    st.session_state.setdefault("qc_pca_reset_nonce", 0)
    st.session_state.setdefault("qc_pca_cache", {})
    st.session_state.setdefault("qc_corr_settings", default_qc_corr_settings())
    st.session_state.setdefault("qc_corr_reset_nonce", 0)
    st.session_state.setdefault("qc_corr_cache", {})
    st.session_state.setdefault("qc_expression_matrix_cache", {})
    st.session_state.setdefault("qc_group_id_counter", 2)
    st.session_state.setdefault("normalization_results", None)
    st.session_state.setdefault("normalization_output_dir", None)
    st.session_state.setdefault("normalization_report", None)
    st.session_state.setdefault("normalization_tables", {})
    st.session_state.setdefault("normalization_run_status", None)
    st.session_state.setdefault("normalization_input_signature", None)
    st.session_state.setdefault("normalization_selected_matrix", "Raw counts")
    st.session_state.setdefault("normalization_table_page", 1)
    st.session_state.setdefault("normalization_table_rows_per_page", 25)
    st.session_state.setdefault("normalization_table_search", "")
    st.session_state.setdefault("normalization_selected_gene", "")


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
        "grouping_set_name": "",
        "groups": [
            {"id": "group_1", "name": "", "samples": []},
            {"id": "group_2", "name": "", "samples": []},
        ],
    }
    st.session_state["qc_group_id_counter"] = 2
    st.session_state["qc_plot_settings"] = default_qc_plot_settings()
    st.session_state["qc_plot_reset_nonce"] = {plot_id: 0 for plot_id in QC_PLOT_DEFINITIONS}
    st.session_state["qc_export_bytes"] = {}
    st.session_state["qc_summary"] = None
    st.session_state["qc_summary_signature"] = None
    st.session_state["qc_plot_data_cache"] = {}
    clear_qc_expression_analysis_cache()
    st.session_state["gene_id_mode"] = "unknown"
    st.session_state["conversion_summary"] = default_conversion_summary()
    st.session_state["duplicate_summary_df"] = pd.DataFrame(columns=["Gene", "Original IDs", "Number of duplicated rows"])
    clear_normalization_state()
    reset_analysis_state()


def reset_analysis_state() -> None:
    """Clear future analysis outputs when upstream data or assignments change."""
    st.session_state["qc_results"] = None
    st.session_state["deg_results"] = None
    st.session_state["pathway_results"] = None
    st.session_state["plots"] = None


def clear_normalization_state() -> None:
    """Clear normalization outputs when the processed count matrix changes."""
    st.session_state["normalization_results"] = None
    st.session_state["normalization_output_dir"] = None
    st.session_state["normalization_report"] = None
    st.session_state["normalization_tables"] = {}
    st.session_state["normalization_run_status"] = None
    st.session_state["normalization_input_signature"] = None
    st.session_state["normalization_selected_matrix"] = "Raw counts"
    st.session_state["normalization_table_page"] = 1
    st.session_state["normalization_table_rows_per_page"] = 25
    st.session_state["normalization_table_search"] = ""
    st.session_state["normalization_selected_gene"] = ""
    clear_qc_expression_analysis_cache()


def clear_qc_expression_analysis_cache() -> None:
    """Clear cached transformed matrices, PCA scores, and sample correlations."""
    st.session_state["qc_expression_matrix_cache"] = {}
    st.session_state["qc_pca_cache"] = {}
    st.session_state["qc_corr_cache"] = {}


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


def default_qc_pca_settings() -> dict[str, Any]:
    """Return default settings for the PCA QC plot."""
    return {
        "matrix": "DESeq2 VST",
        "plot_by": "Sample name",
        "grouping_set": None,
        "component_pair": "PC1 vs PC2",
        "width": 980,
        "height": 560,
        "axis_label_font_size": 14,
        "point_size": 10,
        "show_sample_labels": False,
        "colors": {},
    }


def default_qc_corr_settings() -> dict[str, Any]:
    """Return default settings for the sample-correlation QC plot."""
    return {
        "matrix": "DESeq2 VST",
        "plot_by": "Sample name",
        "grouping_set": None,
        "top_variable_genes": "1000",
        "plot_size": 860,
        "label_size": 12,
        "x_axis_angle": 45,
        "colorbar_thickness": 16,
        "colorbar_length": 0.72,
        "show_correlation_values": False,
        "show_group_annotation": True,
        "colors": {},
    }


def get_default_qc_plot_setting(plot_id: str) -> dict[str, Any]:
    """Return default settings for one QC bar plot."""
    return {
        "plot_by": "Sample name",
        "grouping_set": None,
        "aggregation": "Mean",
        "width": 980,
        "height": 480,
        "x_axis_angle": 45,
        "axis_label_font_size": 14,
        "axis_title_font_size": 14,
        "colors": {},
    }


def sync_qc_axis_label_font_size(settings: dict[str, Any]) -> int:
    """Keep axis-label font size populated for QC plot settings."""
    if "axis_label_font_size" not in settings:
        settings["axis_label_font_size"] = settings.get("axis_title_font_size", 14)
    return int(settings.get("axis_label_font_size", 14))


def reset_qc_plot_setting(plot_id: str) -> None:
    """Reset one QC plot, including Streamlit widget state that can repopulate old values."""
    st.session_state["qc_plot_settings"][plot_id] = get_default_qc_plot_setting(plot_id).copy()
    reset_nonce = st.session_state.setdefault("qc_plot_reset_nonce", {})
    reset_nonce[plot_id] = int(reset_nonce.get(plot_id, 0)) + 1
    widget_suffixes = [
        "plot_by",
        "grouping_set",
        "grouping_set_disabled",
        "aggregation",
        "width",
        "height",
        "x_axis_angle",
        "axis_label_font_size",
        "axis_title_font_size",
    ]
    legacy_prefixes = [f"{plot_id}_{suffix}" for suffix in widget_suffixes]
    nonce_prefixes = [f"{plot_id}*{suffix}*" for suffix in widget_suffixes]
    extra_prefixes = [f"{plot_id}*color*", f"{plot_id}*download*"]
    for key in list(st.session_state.keys()):
        key_text = str(key)
        if any(key_text.startswith(prefix) for prefix in [*legacy_prefixes, *nonce_prefixes, *extra_prefixes]):
            st.session_state.pop(key, None)
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
        return
    cache = st.session_state.setdefault("qc_plot_data_cache", {})
    for key in list(cache):
        if key.startswith(f"{plot_id}:"):
            cache.pop(key, None)


def default_qc_group_editor(sample_columns: list[str]) -> dict[str, Any]:
    """Build a default QC grouping draft."""
    return {
        "grouping_set_name": "",
        "groups": [
            {"id": "group_1", "name": "", "samples": []},
            {"id": "group_2", "name": "", "samples": []},
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
        clean_name = str(group.get("name", "")).strip()
        samples = group.get("samples", [])
        unique_samples = []
        for sample in samples if isinstance(samples, list) else []:
            sample_text = str(sample)
            if sample_text in sample_set and sample_text not in unique_samples:
                unique_samples.append(sample_text)
        normalized_groups.append({"id": group_id, "name": clean_name, "samples": unique_samples})

    while len(normalized_groups) < 2:
        next_index = len(normalized_groups) + 1
        normalized_groups.append({"id": f"group_{next_index}", "name": "", "samples": []})

    max_numeric_id = 0
    for group in normalized_groups:
        match = re.match(r"^group_(\d+)$", str(group["id"]))
        if match:
            max_numeric_id = max(max_numeric_id, int(match.group(1)))
    st.session_state["qc_group_id_counter"] = max(int(st.session_state.get("qc_group_id_counter", 2)), max_numeric_id, 2)

    return {
        "grouping_set_name": str(editor.get("grouping_set_name", "")).strip(),
        "groups": normalized_groups,
    }


def _effective_qc_group_name(group: dict[str, Any], index: int) -> str:
    """Return the saved group name, using placeholder text only at save/validation time."""
    return str(group.get("name", "")).strip() or f"Group {index}"


def get_editor_groups_as_dict(groups_or_editor: Any) -> dict[str, list[str]]:
    """Convert editor groups or saved grouping dict to {group_name: samples}."""
    if isinstance(groups_or_editor, dict) and isinstance(groups_or_editor.get("groups"), list):
        return {
            _effective_qc_group_name(group, index): list(group.get("samples", []))
            for index, group in enumerate(groups_or_editor["groups"], start=1)
        }
    if isinstance(groups_or_editor, list):
        return {
            _effective_qc_group_name(group, index): list(group.get("samples", []))
            for index, group in enumerate(groups_or_editor, start=1)
            if isinstance(group, dict)
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
        group_names = [
            _effective_qc_group_name(group, index)
            for index, group in enumerate(groups_or_editor["groups"], start=1)
        ]
    elif isinstance(groups_or_editor, list):
        group_names = [
            _effective_qc_group_name(group, index)
            for index, group in enumerate(groups_or_editor, start=1)
            if isinstance(group, dict)
        ]
    elif isinstance(groups_or_editor, dict):
        group_names = [str(group_name).strip() for group_name in groups_or_editor]
    else:
        group_names = []
    groups = get_editor_groups_as_dict(groups_or_editor)
    errors = []
    if len(group_names) < 2:
        errors.append("At least two groups are required.")
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
    clear_qc_export_cache()
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
    editor["grouping_set_name"] = st.session_state.get("qc_grouping_set_name_input", "")
    st.session_state["current_qc_group_editor"] = editor


def update_qc_group_name_from_widget(group_id: str) -> None:
    """Sync one group name widget into the editor draft."""
    editor = normalize_qc_group_editor(st.session_state.get("current_qc_group_editor"), st.session_state.get("sample_columns", []))
    widget_key = f"qc_group_name_input_{group_id}"
    for group in editor["groups"]:
        if group["id"] == group_id:
            group["name"] = str(st.session_state.get(widget_key, group["name"])).strip()
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
    editor["groups"].append({"id": group_id, "name": "", "samples": []})
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
    settings = get_qc_visual_settings(plot_id)
    colors = settings.setdefault("colors", {})
    if entity_key not in colors:
        colors[entity_key] = DEFAULT_QC_COLORS[index % len(DEFAULT_QC_COLORS)]
    return colors[entity_key]


def set_qc_color(plot_id: str, entity_key: str, color: str) -> None:
    """Store a custom QC color."""
    settings = get_qc_visual_settings(plot_id)
    settings.setdefault("colors", {})[entity_key] = color


def get_qc_visual_settings(plot_id: str) -> dict[str, Any]:
    """Return the settings dict that owns colors for a QC plot."""
    if plot_id == "pca_plot":
        return st.session_state.setdefault("qc_pca_settings", default_qc_pca_settings())
    if plot_id == "sample_correlation":
        return st.session_state.setdefault("qc_corr_settings", default_qc_corr_settings())
    return st.session_state["qc_plot_settings"].setdefault(plot_id, get_default_qc_plot_setting(plot_id))


def get_qc_reset_nonce(plot_id: str) -> int:
    """Return the widget nonce for a QC plot."""
    if plot_id == "pca_plot":
        return int(st.session_state.setdefault("qc_pca_reset_nonce", 0))
    if plot_id == "sample_correlation":
        return int(st.session_state.setdefault("qc_corr_reset_nonce", 0))
    return int(st.session_state.setdefault("qc_plot_reset_nonce", {}).setdefault(plot_id, 0))


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
    nonce = st.session_state.setdefault("qc_plot_reset_nonce", {}).setdefault(plot_id, 0)
    for index, (entity_key, label) in enumerate(zip(entity_keys, labels)):
        with columns[index % n_cols]:
            color = st.color_picker(
                label,
                value=get_qc_color(plot_id, entity_key, index),
                key=f"{plot_id}*color*{entity_key}_{nonce}",
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


def align_button_with_input() -> None:
    """Align a button with neighboring labeled input widgets."""
    st.markdown("<div style='height: 1.75rem;'></div>", unsafe_allow_html=True)


def render_sidebar() -> None:
    """Render workflow-only sidebar."""
    st.sidebar.title("Workflow")
    st.sidebar.success("1. Upload Count Matrix")
    st.sidebar.success("2. Quality Control")
    st.sidebar.success("3. Normalization")
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
    clear_normalization_state()
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
                    "grouping_set_name": "",
                    "groups": [
                        {"id": "group_1", "name": "", "samples": []},
                        {"id": "group_2", "name": "", "samples": []},
                    ],
                }
                st.session_state["qc_group_id_counter"] = 2
                st.session_state["qc_grouping_sets"] = {}
                st.session_state["active_qc_grouping_set"] = None
                st.session_state["qc_plot_settings"] = default_qc_plot_settings()
                st.session_state["qc_plot_reset_nonce"] = {plot_id: 0 for plot_id in QC_PLOT_DEFINITIONS}
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
            "grouping_set_name": "",
            "groups": [
                {"id": "group_1", "name": "", "samples": []},
                {"id": "group_2", "name": "", "samples": []},
            ],
        }
        st.session_state["qc_group_id_counter"] = 2
        st.session_state["qc_grouping_sets"] = {}
        st.session_state["active_qc_grouping_set"] = None
        st.session_state["qc_plot_settings"] = default_qc_plot_settings()
        st.session_state["qc_plot_reset_nonce"] = {plot_id: 0 for plot_id in QC_PLOT_DEFINITIONS}
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
    axis_label_font_size: int = 14,
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
        xaxis=dict(
            tickangle=x_axis_angle,
            title="",
            showgrid=False,
            tickfont=dict(size=axis_label_font_size),
        ),
        yaxis=dict(
            title=y_axis_title,
            showgrid=True,
            gridcolor="#E5E7EB",
            tickformat=y_tick_format,
            title_font=dict(size=axis_label_font_size),
            tickfont=dict(size=axis_label_font_size),
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
            group["name"] = str(st.session_state[name_key]).strip()
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
            placeholder="QC grouping 1",
            on_change=update_qc_grouping_name_from_widget,
        )

        for index, group in enumerate(editor["groups"], start=1):
            group_id = group["id"]
            name_key = f"qc_group_name_input_{group_id}"
            samples_key = f"qc_group_samples_input_{group_id}"
            options = get_available_samples_for_group(group_id, editor, sample_columns)
            st.session_state[name_key] = group["name"]
            st.session_state[samples_key] = [sample for sample in group.get("samples", []) if sample in options]

            row_cols = st.columns([2.0, 4.8, 1.55], gap="small")
            with row_cols[0]:
                st.text_input(
                    "Group name",
                    key=name_key,
                    placeholder=f"Group {index}",
                    on_change=update_qc_group_name_from_widget,
                    args=(group_id,),
                )
            with row_cols[1]:
                st.multiselect(
                    "Samples",
                    options=options,
                    key=samples_key,
                    on_change=update_qc_group_samples_from_widget,
                    args=(group_id,),
                )
            with row_cols[2]:
                align_button_with_input()
                st.button(
                    "Remove group",
                    key=f"qc_group_remove_{group_id}",
                    disabled=len(editor["groups"]) <= 2,
                    on_click=remove_qc_group_from_editor,
                    args=(group_id, sample_columns),
                    help="Remove this group",
                )

        action_cols = st.columns([2.0, 4.8, 1.55], gap="small")
        with action_cols[0]:
            save_grouping = st.button(
                "Save QC grouping",
                key="qc_group_save_button",
                type="primary",
            )
        with action_cols[1]:
            st.button(
                "Clear grouping info",
                key="qc_group_clear_button",
                on_click=clear_qc_group_editor,
                args=(sample_columns,),
            )
        with action_cols[2]:
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
            st.session_state["qc_group_save_message"] = f"Saved QC grouping: {saved_name}"
            clear_qc_group_editor(sample_columns)
            st.rerun()

    draft_groups = st.session_state["current_qc_group_editor"]["groups"]
    st.write("Unassigned samples:")
    render_inline_badges(get_unassigned_samples(draft_groups, sample_columns))
    duplicates = get_duplicate_assigned_samples(draft_groups)
    if duplicates:
        st.warning(f"Duplicate assignment in current draft: {', '.join(duplicates)}")


def render_qc_grouping_section(sample_columns: list[str]) -> None:
    """Render the immediate-update QC grouping editor and saved-set controls."""
    st.markdown("### Assign QC grouping")
    if st.session_state.get("qc_group_save_message"):
        st.success(st.session_state.pop("qc_group_save_message"))
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
            clear_qc_export_cache()
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
    axis_label_font_size = sync_qc_axis_label_font_size(settings)
    if plot_id == "zero_fraction" and settings.get("aggregation") == "Sum":
        settings["aggregation"] = "Mean"
        st.session_state.pop(f"{plot_id}_aggregation", None)
    nonce = get_qc_reset_nonce(plot_id)

    control_cols = st.columns([1.35, 1.75, 1.05, 0.72, 0.55, 0.55, 1.8], gap="small")
    with control_cols[0]:
        settings["plot_by"] = st.selectbox(
            "Plot by",
            ["Sample name", "QC assignment group"],
            index=["Sample name", "QC assignment group"].index(settings.get("plot_by", "Sample name")),
            key=f"{plot_id}*plot_by*{nonce}",
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
                key=f"{plot_id}*grouping_set*{nonce}",
            )
        else:
            st.selectbox(
                "QC assignment set",
                ["No QC assignment selected"],
                disabled=True,
                key=f"{plot_id}*grouping_set_disabled*{nonce}",
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
            key=f"{plot_id}*aggregation*{nonce}",
        )

    with control_cols[3]:
        align_button_with_input()
        if st.button("Reset", key=f"{plot_id}_reset"):
            reset_qc_plot_setting(plot_id)

    with st.expander("Advanced settings"):
        adv_cols = st.columns(4)
        with adv_cols[0]:
            settings["width"] = st.slider("Plot width", 640, 1440, int(settings.get("width", 980)), 20, key=f"{plot_id}*width*{nonce}")
        with adv_cols[1]:
            settings["height"] = st.slider("Plot height", 360, 900, int(settings.get("height", 480)), 20, key=f"{plot_id}*height*{nonce}")
        with adv_cols[2]:
            settings["axis_label_font_size"] = st.slider(
                "Axis label font size",
                10,
                28,
                axis_label_font_size,
                1,
                key=f"{plot_id}*axis_label_font_size*{nonce}",
            )
        with adv_cols[3]:
            angle_options = [0, 30, 45, 60, 90]
            settings["x_axis_angle"] = st.selectbox(
                "X-axis angle",
                angle_options,
                index=angle_options.index(int(settings.get("x_axis_angle", 45))),
                key=f"{plot_id}*x_axis_angle*{nonce}",
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
        axis_label_font_size=int(settings.get("axis_label_font_size", 14)),
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
                "axis_label_font_size": settings.get("axis_label_font_size"),
                "colors": settings.get("colors", {}),
                "qc_summary_signature": st.session_state.get("qc_summary_signature"),
            },
            sort_keys=True,
        ).encode("utf-8")
    ).hexdigest()[:8]
    filename_base = f"{plot_id}_{APP_VERSION}_{export_signature}"
    render_qc_export_buttons(fig, plot_id, filename_base, export_signature, control_cols[4], control_cols[5], nonce)
    st.plotly_chart(fig, use_container_width=False, key=f"{plot_id}_plotly_chart")


def get_cached_plot_export_bytes(plot_id: str, export_signature: str, fig: go.Figure, fmt: str) -> bytes | None:
    """Return cached Plotly static image bytes for one plot state."""
    cache_key = f"{plot_id}:{export_signature}:{fmt}"
    export_cache = st.session_state.setdefault("qc_export_bytes", {})
    if cache_key in export_cache:
        return export_cache[cache_key]
    try:
        kwargs: dict[str, Any] = {"format": fmt}
        if fmt == "png":
            kwargs["scale"] = 3
        export_bytes = fig.to_image(**kwargs)
        if isinstance(export_bytes, str):
            export_bytes = export_bytes.encode("utf-8")
        export_cache[cache_key] = export_bytes
        return export_bytes
    except Exception:
        return None


def render_qc_export_buttons(
    fig: go.Figure,
    plot_id: str,
    filename_base: str,
    export_signature: str,
    png_column,
    svg_column,
    nonce: int,
) -> None:
    """Render lazy PNG/SVG export controls for one QC plot."""
    export_cache = st.session_state.setdefault("qc_export_bytes", {})

    def render_one_export(column, fmt: str, mime: str) -> None:
        cache_key = f"{plot_id}:{export_signature}:{fmt}"
        export_bytes = export_cache.get(cache_key)
        with column:
            align_button_with_input()
            if export_bytes is None:
                if st.button(fmt.upper(), key=f"{plot_id}*prepare*{fmt}*{filename_base}*{nonce}", help=f"Prepare {fmt.upper()} export"):
                    export_bytes = get_cached_plot_export_bytes(plot_id, export_signature, fig, fmt)
                    if export_bytes is None:
                        st.warning("Static image export requires kaleido. Install with: pip install kaleido")
                    else:
                        st.rerun()
                return
            st.download_button(
                fmt.upper(),
                data=export_bytes,
                file_name=f"{filename_base}.{fmt}",
                mime=mime,
                key=f"{plot_id}*download*{fmt}*{filename_base}*{nonce}",
                help=f"Download {fmt.upper()}",
            )

    render_one_export(png_column, "png", "image/png")
    render_one_export(svg_column, "svg", "image/svg+xml")


def validate_normalization_input(counts_df: pd.DataFrame | None, sample_columns: list[str]) -> dict[str, Any]:
    """Validate the processed raw count matrix used by normalization."""
    result: dict[str, Any] = {
        "valid": False,
        "errors": [],
        "warnings": [],
        "summary": {},
        "numeric_counts": None,
    }
    if counts_df is None or counts_df.empty:
        result["errors"].append("Processed count matrix is empty.")
        return result
    if "Gene" not in counts_df.columns:
        result["errors"].append("Processed count matrix must contain a Gene column.")
        return result
    if len(sample_columns) < 2:
        result["errors"].append("At least two sample columns are required for normalization.")
        return result
    missing_samples = [sample for sample in sample_columns if sample not in counts_df.columns]
    if missing_samples:
        result["errors"].append(f"Sample columns missing from processed matrix: {', '.join(missing_samples[:8])}")
        return result

    gene_series = counts_df["Gene"].astype("string").fillna("").str.strip()
    empty_genes = int((gene_series == "").sum())
    duplicated_genes = int(gene_series[gene_series != ""].duplicated().sum())
    numeric_counts = counts_df[sample_columns].apply(pd.to_numeric, errors="coerce")
    if numeric_counts.isna().any().any():
        result["errors"].append("Non-numeric count values detected.")
    if (numeric_counts < 0).any().any():
        result["errors"].append("Negative count values detected.")
    if empty_genes:
        result["errors"].append(f"Empty gene IDs detected: {empty_genes}")

    filled_counts = numeric_counts.fillna(0)
    all_zero_genes = int((filled_counts.sum(axis=1) == 0).sum())
    total_counts = filled_counts.sum(axis=0)
    deviations = (filled_counts - filled_counts.round()).abs()
    max_deviation = float(deviations.max().max()) if not deviations.empty else 0.0
    non_integer_mask = deviations > 1e-6
    non_integer_fraction = float(non_integer_mask.to_numpy().sum() / filled_counts.size) if filled_counts.size else 0.0
    if max_deviation > 1e-6 or non_integer_fraction > 0:
        result["warnings"].append(
            "Non-integer count values detected. DESeq2 and edgeR will use rounded counts in the R workflow."
        )
    sample_totals_near_cpm = bool(((total_counts - 1_000_000).abs() < 50_000).all()) if len(total_counts) else False
    if sample_totals_near_cpm and non_integer_fraction > 0.05:
        result["warnings"].append("Input may be TPM/CPM-like because sample totals are near 1e6 and decimals are common.")

    result["valid"] = not result["errors"]
    result["summary"] = {
        "genes": int(counts_df.shape[0]),
        "samples": int(len(sample_columns)),
        "duplicated_gene_ids": duplicated_genes,
        "all_zero_genes": all_zero_genes,
        "max_deviation_from_integer": max_deviation,
        "non_integer_value_fraction": non_integer_fraction,
    }
    result["numeric_counts"] = numeric_counts
    return result


def get_normalization_input_signature() -> str:
    """Return a signature for the current processed matrix and selected sample columns."""
    processed_counts_df = st.session_state.get("processed_counts_df")
    sample_columns = st.session_state.get("sample_columns", [])
    processed_hash = None
    if processed_counts_df is not None and "Gene" in processed_counts_df.columns and sample_columns:
        signature_df = processed_counts_df[["Gene", *sample_columns]].copy()
        processed_hash = hashlib.sha256(
            pd.util.hash_pandas_object(signature_df, index=True).to_numpy().tobytes()
        ).hexdigest()
    payload = {
        "counts": st.session_state.get("counts_file_signature"),
        "gene_id_column": st.session_state.get("gene_id_column"),
        "samples": sample_columns,
        "processed_shape": getattr(processed_counts_df, "shape", None),
        "processed_hash": processed_hash,
    }
    return hashlib.sha256(json.dumps(payload, sort_keys=True, default=str).encode("utf-8")).hexdigest()


def create_normalization_output_dir() -> Path:
    """Create a unique output directory for one normalization run."""
    base_dir = Path.cwd() / "outputs" / "normalization"
    base_dir.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    signature = st.session_state.get("counts_file_signature")
    if signature:
        short_hash = str(signature)[:8]
    else:
        short_hash = hashlib.sha256(f"{timestamp}:{tempfile.gettempdir()}".encode("utf-8")).hexdigest()[:8]
    output_dir = base_dir / f"{APP_VERSION}_{timestamp}_{short_hash}"
    output_dir.mkdir(parents=True, exist_ok=False)
    return output_dir


def prepare_normalization_input(
    processed_counts_df: pd.DataFrame,
    sample_columns: list[str],
    output_dir: Path,
) -> Path:
    """Write the processed raw count matrix that will be consumed by the R workflow."""
    output_dir.mkdir(parents=True, exist_ok=True)
    input_df = processed_counts_df[["Gene", *sample_columns]].copy()
    input_df["Gene"] = input_df["Gene"].astype(str)
    for sample in sample_columns:
        input_df[sample] = pd.to_numeric(input_df[sample], errors="coerce")
    input_path = output_dir / "input_processed_raw_counts.csv"
    input_df.to_csv(input_path, index=False)
    return input_path


def run_r_normalization(
    counts_path: Path,
    output_dir: Path,
    gene_col: str = "Gene",
    prior_count: float = 1.0,
) -> dict[str, Any]:
    """Run the standalone R normalization script and capture its process output."""
    try:
        subprocess.run(["Rscript", "--version"], capture_output=True, text=True, timeout=30)
    except FileNotFoundError:
        return {
            "success": False,
            "stdout": "",
            "stderr": "Rscript was not found. Please install R and make sure Rscript is available in PATH.",
            "returncode": None,
            "command": ["Rscript", "--version"],
        }
    except Exception as exc:
        return {
            "success": False,
            "stdout": "",
            "stderr": f"Could not check Rscript availability: {exc}",
            "returncode": None,
            "command": ["Rscript", "--version"],
        }

    script_path = Path.cwd() / "r_scripts" / "normalize_counts.R"
    command = [
        "Rscript",
        str(script_path),
        "--counts",
        str(counts_path),
        "--outdir",
        str(output_dir),
        "--gene_col",
        gene_col,
        "--prior_count",
        str(prior_count),
    ]
    try:
        completed = subprocess.run(command, capture_output=True, text=True, timeout=600)
        return {
            "success": completed.returncode == 0,
            "stdout": completed.stdout,
            "stderr": completed.stderr,
            "returncode": completed.returncode,
            "command": command,
        }
    except subprocess.TimeoutExpired as exc:
        return {
            "success": False,
            "stdout": exc.stdout or "",
            "stderr": "R normalization timed out after 600 seconds.",
            "returncode": None,
            "command": command,
        }
    except Exception as exc:
        return {
            "success": False,
            "stdout": "",
            "stderr": f"R normalization failed before completion: {exc}",
            "returncode": None,
            "command": command,
        }


def load_normalization_report(output_dir: Path) -> dict[str, Any] | None:
    """Load normalization_report.json if it exists."""
    report_path = output_dir / "normalization_report.json"
    if not report_path.exists():
        return None
    return json.loads(report_path.read_text(encoding="utf-8"))


def load_normalization_table(output_dir: Path, label: str) -> pd.DataFrame:
    """Load and cache one normalization matrix table by label."""
    filename = NORMALIZATION_OUTPUTS[label]
    cache_key = f"{output_dir.resolve()}:{filename}"
    cache = st.session_state.setdefault("normalization_tables", {})
    if cache_key not in cache:
        cache[cache_key] = pd.read_csv(output_dir / filename)
    return cache[cache_key]


def _matrix_signature(df: pd.DataFrame, sample_columns: list[str]) -> str:
    """Return a stable hash for a Gene + sample expression matrix."""
    signature_df = df[["Gene", *sample_columns]].copy()
    return hashlib.sha256(pd.util.hash_pandas_object(signature_df, index=True).to_numpy().tobytes()).hexdigest()


def _align_expression_matrix(df: pd.DataFrame, sample_columns: list[str]) -> pd.DataFrame | None:
    """Return a Gene + sample matrix aligned to current sample columns."""
    if df is None or df.empty or "Gene" not in df.columns:
        return None
    if any(sample not in df.columns for sample in sample_columns):
        return None
    aligned = df[["Gene", *sample_columns]].copy()
    aligned["Gene"] = aligned["Gene"].astype(str)
    for sample in sample_columns:
        aligned[sample] = pd.to_numeric(aligned[sample], errors="coerce")
    aligned = aligned.replace([np.inf, -np.inf], np.nan).dropna(subset=sample_columns)
    return aligned


def get_available_expression_matrices() -> dict[str, dict[str, Any]]:
    """Return cached expression matrices available for PCA and sample correlation."""
    processed_counts_df = st.session_state.get("processed_counts_df")
    sample_columns = st.session_state.get("sample_columns", [])
    if processed_counts_df is None or not sample_columns or np is None:
        return {}
    input_signature = get_normalization_input_signature()
    output_dir = None
    results = st.session_state.get("normalization_results")
    if results and st.session_state.get("normalization_input_signature") == input_signature:
        output_dir_text = results.get("output_dir")
        output_dir = Path(output_dir_text) if output_dir_text else None
    cache_signature = {
        "input_signature": input_signature,
        "output_dir": str(output_dir) if output_dir else None,
    }
    cache_key = hashlib.sha256(json.dumps(cache_signature, sort_keys=True).encode("utf-8")).hexdigest()
    cache = st.session_state.setdefault("qc_expression_matrix_cache", {})
    if cache.get("cache_key") == cache_key:
        return cache.get("matrices", {})

    matrices: dict[str, dict[str, Any]] = {}
    if output_dir and output_dir.exists():
        for label in ["DESeq2 VST", "log2(CPM + 1)", "edgeR TMM logCPM", "DESeq2 normalized counts"]:
            try:
                matrix_df = load_normalization_table(output_dir, label)
            except Exception:
                continue
            aligned = _align_expression_matrix(matrix_df, sample_columns)
            if aligned is not None and not aligned.empty:
                matrices[label] = {
                    "df": aligned,
                    "signature": _matrix_signature(aligned, sample_columns),
                    "source": "normalization output",
                }

    ordered_labels = ["DESeq2 VST", "log2(CPM + 1)", "edgeR TMM logCPM", "DESeq2 normalized counts"]
    ordered_matrices = {label: matrices[label] for label in ordered_labels if label in matrices}
    cache["cache_key"] = cache_key
    cache["matrices"] = ordered_matrices
    return ordered_matrices


def select_effective_expression_values(
    matrix_df: pd.DataFrame,
    sample_columns: list[str],
    top_n: int | None = None,
) -> tuple[Any, int]:
    """Return finite non-constant gene rows, optionally limited by variance."""
    values = matrix_df[sample_columns].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    finite_mask = np.isfinite(values).all(axis=1)
    values = values[finite_mask]
    if values.size == 0:
        return values, 0
    variances = np.var(values, axis=1)
    effective_mask = variances > 1e-12
    values = values[effective_mask]
    variances = variances[effective_mask]
    if values.shape[0] == 0:
        return values, 0
    if top_n is not None and values.shape[0] > top_n:
        top_indices = np.argsort(variances)[::-1][:top_n]
        values = values[top_indices]
    return values, int(values.shape[0])


def parse_top_variable_gene_limit(value: str) -> int | None:
    """Parse top-variable-gene setting."""
    return None if str(value).upper() == "ALL" else int(value)


def get_sample_group_assignments(grouping_name: str | None, sample_columns: list[str]) -> dict[str, str]:
    """Map samples to saved QC groups, using Unassigned for missing samples."""
    assignments = {sample: "Unassigned" for sample in sample_columns}
    grouping_dict = st.session_state.get("qc_grouping_sets", {}).get(grouping_name or "")
    if not grouping_dict:
        return assignments
    for group_name, samples in grouping_dict.items():
        for sample in samples:
            if sample in assignments:
                assignments[sample] = group_name
    return assignments


def get_group_labels_for_samples(grouping_name: str | None, sample_columns: list[str]) -> list[str]:
    """Return group labels in sample order."""
    assignments = get_sample_group_assignments(grouping_name, sample_columns)
    return [assignments.get(sample, "Unassigned") for sample in sample_columns]


def get_qc_group_color_key(group_name: str) -> str:
    """Return the shared color key for a QC group label."""
    return f"group:{group_name}"


def build_pca_component_pairs(max_pc: int) -> list[str]:
    """Build available sequential PCA component-pair labels."""
    limit = min(max(max_pc - 1, 0), 9)
    return [f"PC{index} vs PC{index + 1}" for index in range(1, limit + 1)]


def parse_pca_component_pair(pair_label: str) -> tuple[int, int]:
    """Parse a component-pair label into zero-based component indices."""
    match = re.match(r"^PC(\d+) vs PC(\d+)$", pair_label)
    if not match:
        return 0, 1
    return int(match.group(1)) - 1, int(match.group(2)) - 1


def compute_cached_pca(matrix_name: str, matrix_info: dict[str, Any], sample_columns: list[str]) -> dict[str, Any]:
    """Compute or return cached PCA scores for one expression matrix."""
    cache_payload = {
        "matrix": matrix_name,
        "signature": matrix_info.get("signature"),
        "top_variable_genes": 500,
    }
    cache_key = hashlib.sha256(json.dumps(cache_payload, sort_keys=True).encode("utf-8")).hexdigest()
    cache = st.session_state.setdefault("qc_pca_cache", {})
    if cache_key in cache:
        return cache[cache_key]
    matrix_df = matrix_info["df"]
    if len(sample_columns) < 3:
        result = {"error": "At least 3 samples are required for PCA component-pair plots."}
        cache[cache_key] = result
        return result
    values, effective_genes = select_effective_expression_values(matrix_df, sample_columns, top_n=500)
    if effective_genes < 2:
        result = {"error": "At least 2 non-constant genes are required for PCA."}
        cache[cache_key] = result
        return result
    sample_by_gene = values.T
    centered = sample_by_gene - np.mean(sample_by_gene, axis=0, keepdims=True)
    try:
        u_matrix, singular_values, _ = np.linalg.svd(centered, full_matrices=False)
    except Exception as exc:
        result = {"error": f"PCA could not be calculated: {exc}"}
        cache[cache_key] = result
        return result
    max_pc = int(min(len(sample_columns) - 1, effective_genes, len(singular_values), 10))
    if max_pc < 2:
        result = {"error": "At least 2 principal components are required for PCA plotting."}
        cache[cache_key] = result
        return result
    scores = u_matrix[:, :max_pc] * singular_values[:max_pc]
    eigenvalues = (singular_values[:max_pc] ** 2) / max(len(sample_columns) - 1, 1)
    total_variance = float(np.sum((singular_values ** 2) / max(len(sample_columns) - 1, 1)))
    explained = (eigenvalues / total_variance * 100) if total_variance > 0 else np.zeros_like(eigenvalues)
    result = {
        "scores": scores,
        "explained_variance": explained,
        "max_pc": max_pc,
        "effective_genes": effective_genes,
        "component_pairs": build_pca_component_pairs(max_pc),
    }
    cache[cache_key] = result
    return result


def compute_cached_sample_correlation(
    matrix_name: str,
    matrix_info: dict[str, Any],
    sample_columns: list[str],
    top_variable_genes: str,
) -> dict[str, Any]:
    """Compute or return cached sample-by-sample Pearson correlation."""
    top_n = parse_top_variable_gene_limit(top_variable_genes)
    cache_payload = {
        "matrix": matrix_name,
        "signature": matrix_info.get("signature"),
        "top_variable_genes": top_variable_genes,
    }
    cache_key = hashlib.sha256(json.dumps(cache_payload, sort_keys=True).encode("utf-8")).hexdigest()
    cache = st.session_state.setdefault("qc_corr_cache", {})
    if cache_key in cache:
        return cache[cache_key]
    values, effective_genes = select_effective_expression_values(matrix_info["df"], sample_columns, top_n=top_n)
    if len(sample_columns) < 2:
        result = {"error": "At least 2 samples are required for sample correlation."}
        cache[cache_key] = result
        return result
    if effective_genes < 2:
        result = {"error": "At least 2 non-constant genes are required for sample correlation."}
        cache[cache_key] = result
        return result
    corr = np.corrcoef(values.T)
    corr = np.nan_to_num(corr, nan=0.0, posinf=0.0, neginf=0.0)
    np.fill_diagonal(corr, 1.0)
    result = {"correlation": corr, "effective_genes": effective_genes}
    cache[cache_key] = result
    return result


def make_qc_pca_plot(
    pca_result: dict[str, Any],
    sample_columns: list[str],
    settings: dict[str, Any],
    grouping_name: str | None,
    width: int,
    height: int,
) -> tuple[go.Figure, list[str], list[str], str]:
    """Create the PCA Plotly scatter figure and return color entities."""
    x_index, y_index = parse_pca_component_pair(settings.get("component_pair", "PC1 vs PC2"))
    scores = pca_result["scores"]
    explained = pca_result["explained_variance"]
    plot_by_group = settings.get("plot_by") == "QC assignment group"
    groups = get_group_labels_for_samples(grouping_name, sample_columns)
    color_labels = groups if plot_by_group else sample_columns
    color_keys = [get_qc_group_color_key(group) if plot_by_group else f"sample:{sample}" for group, sample in zip(groups, sample_columns)]
    unique_keys: list[str] = []
    unique_labels: list[str] = []
    for key, label in zip(color_keys, color_labels):
        if key not in unique_keys:
            unique_keys.append(key)
            unique_labels.append(label)
    color_map = {key: get_qc_color("pca_plot", key, index) for index, key in enumerate(unique_keys)}

    fig = go.Figure()
    for key, label in zip(unique_keys, unique_labels):
        indices = [index for index, color_key in enumerate(color_keys) if color_key == key]
        customdata = [
            [sample_columns[index], groups[index], scores[index, x_index], scores[index, y_index]]
            for index in indices
        ]
        fig.add_trace(
            go.Scatter(
                x=[scores[index, x_index] for index in indices],
                y=[scores[index, y_index] for index in indices],
                mode="markers+text" if settings.get("show_sample_labels") else "markers",
                text=[sample_columns[index] for index in indices] if settings.get("show_sample_labels") else None,
                textposition="top center",
                marker=dict(
                    size=int(settings.get("point_size", 10)),
                    color=color_map[key],
                    line=dict(width=1, color="rgba(17, 24, 39, 0.55)"),
                ),
                customdata=customdata,
                hovertemplate=(
                    "Sample: %{customdata[0]}<br>"
                    + ("Group: %{customdata[1]}<br>" if plot_by_group else "")
                    + f"PC{x_index + 1}: %{{customdata[2]:.4g}}<br>"
                    + f"PC{y_index + 1}: %{{customdata[3]:.4g}}<extra></extra>"
                ),
                name=str(label),
            )
        )
    title = f"PCA - PC{x_index + 1} vs PC{y_index + 1}"
    fig.update_layout(
        title=title,
        width=width,
        height=height,
        plot_bgcolor="white",
        paper_bgcolor="white",
        margin=dict(l=75, r=40, t=70, b=75),
        xaxis=dict(
            title=f"PC{x_index + 1} ({explained[x_index]:.1f}% variance)",
            showgrid=True,
            gridcolor="#E5E7EB",
            zeroline=True,
            title_font=dict(size=int(settings.get("axis_label_font_size", 14))),
            tickfont=dict(size=int(settings.get("axis_label_font_size", 14))),
        ),
        yaxis=dict(
            title=f"PC{y_index + 1} ({explained[y_index]:.1f}% variance)",
            showgrid=True,
            gridcolor="#E5E7EB",
            zeroline=True,
            title_font=dict(size=int(settings.get("axis_label_font_size", 14))),
            tickfont=dict(size=int(settings.get("axis_label_font_size", 14))),
        ),
        legend=dict(title="QC assignment group" if plot_by_group else "Sample"),
        font=dict(size=13),
    )
    return fig, unique_keys, unique_labels, title


def _discrete_colorscale(colors: list[str]) -> list[list[Any]]:
    """Build a Plotly colorscale that maps integer group codes to fixed colors."""
    if not colors:
        return [[0, "#9CA3AF"], [1, "#9CA3AF"]]
    if len(colors) == 1:
        return [[0, colors[0]], [1, colors[0]]]
    scale = []
    denom = len(colors) - 1
    for index, color in enumerate(colors):
        pos = index / denom
        scale.append([pos, color])
    return scale


def make_sample_correlation_plot(
    corr_result: dict[str, Any],
    sample_columns: list[str],
    settings: dict[str, Any],
    grouping_name: str | None,
) -> tuple[go.Figure, list[str], list[str]]:
    """Create the sample-correlation Plotly heatmap."""
    corr = corr_result["correlation"]
    plot_by_group = settings.get("plot_by") == "QC assignment group"
    show_annotation = bool(settings.get("show_group_annotation", True)) and plot_by_group
    groups = get_group_labels_for_samples(grouping_name, sample_columns)
    unique_groups: list[str] = []
    for group in groups:
        if group not in unique_groups:
            unique_groups.append(group)
    group_keys = [get_qc_group_color_key(group) for group in unique_groups]
    group_colors = [get_qc_color("sample_correlation", key, index) for index, key in enumerate(group_keys)]
    group_index = {group: index for index, group in enumerate(unique_groups)}
    group_codes = [group_index[group] for group in groups]

    text = [[f"{value:.2f}" for value in row] for row in corr] if settings.get("show_correlation_values") else None
    fig = go.Figure()
    fig.add_trace(
        go.Heatmap(
            z=corr,
            x=sample_columns,
            y=sample_columns,
            xgap=1,
            ygap=1,
            zmin=-1,
            zmax=1,
            colorscale=[[0, "#1f3b82"], [0.5, "#f8fafc"], [1, "#c62828"]],
            colorbar=dict(
                title="Pearson r",
                thickness=int(settings.get("colorbar_thickness", 16)),
                len=float(settings.get("colorbar_length", 0.72)),
                x=0.88,
            ),
            text=text,
            texttemplate="%{text}" if text else None,
            hovertemplate="Sample 1: %{y}<br>Sample 2: %{x}<br>Pearson r: %{z:.3f}<extra></extra>",
            xaxis="x",
            yaxis="y",
        )
    )
    if show_annotation and unique_groups:
        fig.add_trace(
            go.Heatmap(
                z=[group_codes],
                x=sample_columns,
                y=["Group"],
                xgap=1,
                ygap=1,
                zmin=0,
                zmax=max(len(unique_groups) - 1, 1),
                colorscale=_discrete_colorscale(group_colors),
                showscale=False,
                hovertemplate="Sample: %{x}<br>Group: %{customdata}<extra></extra>",
                customdata=[groups],
                xaxis="x",
                yaxis="y2",
            )
        )
        fig.add_trace(
            go.Heatmap(
                z=[[code] for code in group_codes],
                x=["Group"],
                y=sample_columns,
                xgap=1,
                ygap=1,
                zmin=0,
                zmax=max(len(unique_groups) - 1, 1),
                colorscale=_discrete_colorscale(group_colors),
                showscale=False,
                hovertemplate="Sample: %{y}<br>Group: %{customdata}<extra></extra>",
                customdata=[[group] for group in groups],
                xaxis="x2",
                yaxis="y",
            )
        )
        for group, color in zip(unique_groups, group_colors):
            fig.add_trace(
                go.Scatter(
                    x=[None],
                    y=[None],
                    mode="markers",
                    marker=dict(size=10, color=color),
                    name=group,
                    showlegend=True,
                )
            )
    annotation_fraction = 0.03 if show_annotation else 0.0
    annotation_gap = 0.012 if show_annotation else 0.0
    main_start = annotation_fraction + annotation_gap if show_annotation else 0.0
    main_end = 0.84 if show_annotation else 0.86
    main_y_end = 0.89 if show_annotation else 1.0
    top_y_start = main_y_end + annotation_gap
    top_y_end = min(top_y_start + annotation_fraction, 1.0)
    fig.update_layout(
        title="Sample Correlation",
        width=int(settings.get("plot_size", 860)),
        height=int(settings.get("plot_size", 860)),
        plot_bgcolor="white",
        paper_bgcolor="white",
        margin=dict(l=120, r=145, t=80, b=105),
        xaxis=dict(
            domain=[main_start, main_end],
            tickangle=int(settings.get("x_axis_angle", 45)),
            tickfont=dict(size=int(settings.get("label_size", 12))),
            side="bottom",
        ),
        yaxis=dict(
            domain=[0, main_y_end] if show_annotation else [0, 1],
            autorange="reversed",
            tickfont=dict(size=int(settings.get("label_size", 12))),
        ),
        xaxis2=dict(domain=[0, annotation_fraction], showticklabels=False, showgrid=False, zeroline=False),
        yaxis2=dict(
            domain=[top_y_start, top_y_end] if show_annotation else [1, 1],
            showticklabels=False,
            showgrid=False,
            zeroline=False,
        ),
        legend=dict(title="Group" if show_annotation else None, x=1.04, y=0.98),
        font=dict(size=13),
    )
    return fig, group_keys if show_annotation else [], unique_groups if show_annotation else []


def run_normalization_workflow() -> None:
    """Prepare input, run R normalization, and store output metadata in session state."""
    processed_counts_df = st.session_state.get("processed_counts_df")
    sample_columns = st.session_state.get("sample_columns", [])
    input_signature = get_normalization_input_signature()
    validation = validate_normalization_input(processed_counts_df, sample_columns)
    if not validation["valid"]:
        st.session_state["normalization_run_status"] = {
            "success": False,
            "stderr": "\n".join(validation["errors"]),
            "stdout": "",
            "returncode": None,
            "command": [],
        }
        st.session_state["normalization_input_signature"] = input_signature
        return

    output_dir = create_normalization_output_dir()
    input_path = prepare_normalization_input(processed_counts_df, sample_columns, output_dir)
    st.session_state["normalization_run_status"] = {
        "success": None,
        "stderr": "",
        "stdout": "",
        "returncode": None,
        "command": [],
    }
    run_status = run_r_normalization(input_path, output_dir)
    st.session_state["normalization_run_status"] = run_status
    st.session_state["normalization_output_dir"] = str(output_dir)
    st.session_state["normalization_input_signature"] = input_signature
    st.session_state["normalization_tables"] = {}
    if run_status["success"]:
        clear_qc_expression_analysis_cache()
        report = load_normalization_report(output_dir)
        factor_tables = {}
        for label, filename in NORMALIZATION_FACTOR_FILES.items():
            path = output_dir / filename
            factor_tables[label] = pd.read_csv(path) if path.exists() else pd.DataFrame()
        st.session_state["normalization_report"] = report
        st.session_state["normalization_results"] = {
            "output_dir": str(output_dir),
            "factor_tables": factor_tables,
        }
        st.session_state["normalization_selected_matrix"] = "Raw counts"
        st.session_state["normalization_table_page"] = 1
    else:
        clear_qc_expression_analysis_cache()
        st.session_state["normalization_results"] = None
        st.session_state["normalization_report"] = None


def render_matrix_table_viewer(
    table_df: pd.DataFrame,
    table_name: str,
    download_filename: str,
    key_prefix: str,
) -> None:
    """Render a searchable, paginated matrix preview with CSV download."""
    if table_df.empty:
        st.warning(f"{table_name} is empty.")
        return
    if key_prefix == "normalization":
        search_key = "normalization_table_search"
        rows_key = "normalization_table_rows_per_page"
        page_key = "normalization_table_page"
        selected_gene_key = "normalization_selected_gene"
    else:
        search_key = f"{key_prefix}_search"
        rows_key = f"{key_prefix}_rows_per_page"
        page_key = f"{key_prefix}_page"
        selected_gene_key = f"{key_prefix}_selected_gene"
    previous_search_key = f"{key_prefix}_previous_search"
    previous_rows_key = f"{key_prefix}_previous_rows"
    previous_gene_key = f"{key_prefix}_previous_gene"
    st.session_state.setdefault(search_key, "")
    st.session_state.setdefault(rows_key, 25)
    st.session_state.setdefault(page_key, 1)
    st.session_state.setdefault(selected_gene_key, "")

    top_cols = st.columns([2.2, 0.9, 1.7], gap="small")
    with top_cols[0]:
        search_query = st.text_input(
            "Search genes",
            key=search_key,
            placeholder="Type at least 3 characters...",
        )
        st.caption("Type at least 3 characters to search. Use the matched gene selector to lock one gene.")
    with top_cols[1]:
        rows_per_page = st.selectbox(
            "Show rows",
            [10, 25, 50, 100, 250],
            index=[10, 25, 50, 100, 250].index(int(st.session_state.get(rows_key, 25))),
            key=rows_key,
        )

    gene_series = table_df["Gene"].astype(str) if "Gene" in table_df.columns else pd.Series([], dtype=str)
    search_text = search_query.strip()
    search_active = len(search_text) >= 3
    matched_gene_list: list[str] = []
    if search_active:
        mask = gene_series.str.contains(re.escape(search_text), case=False, na=False)
        matched_df = table_df.loc[mask].copy()
        matched_gene_list = matched_df["Gene"].astype(str).tolist() if "Gene" in matched_df.columns else []
        suggestions = matched_gene_list[:25]
        if suggestions:
            st.caption("Matched genes: " + " / ".join(f"`{gene}`" for gene in suggestions))
            if len(matched_gene_list) > 25:
                st.caption(f"Showing first 25 of {len(matched_gene_list):,} matched genes.")
        else:
            st.caption("No matched genes found.")
    else:
        matched_df = table_df
        st.session_state[selected_gene_key] = ""

    selector_options = [""] + matched_gene_list[:200] if search_active and matched_gene_list else [""]
    if st.session_state.get(selected_gene_key, "") not in selector_options:
        st.session_state[selected_gene_key] = ""
    with top_cols[2]:
        selected_gene = st.selectbox(
            "Jump to matched gene",
            selector_options,
            key=selected_gene_key,
            disabled=not (search_active and matched_gene_list),
        )

    if selected_gene:
        display_df = table_df.loc[gene_series == selected_gene].copy()
    else:
        display_df = matched_df

    if (
        st.session_state.get(previous_search_key) != search_query
        or st.session_state.get(previous_rows_key) != rows_per_page
        or st.session_state.get(previous_gene_key) != selected_gene
    ):
        st.session_state[page_key] = 1
        st.session_state[previous_search_key] = search_query
        st.session_state[previous_rows_key] = rows_per_page
        st.session_state[previous_gene_key] = selected_gene

    total_rows = int(display_df.shape[0])
    rows_per_page = int(rows_per_page)
    total_pages = max((total_rows + rows_per_page - 1) // rows_per_page, 1)
    st.session_state[page_key] = min(max(int(st.session_state.get(page_key, 1)), 1), total_pages)
    page = int(st.session_state[page_key])
    start = (page - 1) * rows_per_page
    end = min(start + rows_per_page, total_rows)
    page_df = display_df.iloc[start:end].copy()

    st.write(f"Matrix shape: `{table_df.shape[0]:,}` genes x `{max(table_df.shape[1] - 1, 0):,}` samples")
    if selected_gene:
        st.write(f"Showing selected gene: `{selected_gene}`")
    elif search_active:
        st.write(f"Search matched `{total_rows:,}` genes")
    st.write(f"Showing rows `{start + 1 if total_rows else 0:,}`-`{end:,}` of `{total_rows:,}`")

    nav_cols = st.columns([0.7, 0.7, 0.6, 1.1, 3.0], gap="small")
    with nav_cols[0]:
        if st.button("Previous", key=f"{key_prefix}_previous", disabled=page <= 1):
            st.session_state[page_key] = max(page - 1, 1)
            st.rerun()
    with nav_cols[1]:
        if st.button("Next", key=f"{key_prefix}_next", disabled=page >= total_pages):
            st.session_state[page_key] = min(page + 1, total_pages)
            st.rerun()
    with nav_cols[2]:
        st.download_button(
            "CSV",
            data=display_df.to_csv(index=False).encode("utf-8"),
            file_name=download_filename,
            mime="text/csv",
            key=f"{key_prefix}_csv",
        )
    with nav_cols[3]:
        st.write(f"Page {page} of {total_pages}")
    st.dataframe(page_df, use_container_width=True, hide_index=True)


def render_normalization_report(report: dict[str, Any] | None) -> None:
    """Render a compact normalization report summary."""
    if not report:
        return
    st.markdown("### Normalization report")
    report_fields = [
        ("Original genes", "original_genes"),
        ("Genes after duplicate merge", "genes_after_duplicate_merge"),
        ("Genes used after zero filtering", "genes_used_after_zero_filtering"),
        ("Samples", "samples"),
        ("R version", "r_version"),
        ("DESeq2 version", "deseq2_version"),
        ("edgeR version", "edger_version"),
    ]
    for label, key in report_fields:
        render_info_line(label, report.get(key, "Not available"))


def render_normalization_tab() -> None:
    """Render the normalization workflow and output matrix viewer."""
    st.subheader("Normalization")
    st.write("Normalization uses the processed raw count matrix after gene-symbol conversion and duplicate-gene merging.")
    if np is None:
        st.error("NumPy is required for normalization input checks. Install with: pip install numpy")
        return

    processed_counts_df = st.session_state.get("processed_counts_df")
    sample_columns = st.session_state.get("sample_columns", [])
    if processed_counts_df is None:
        st.info("Please upload and process a count matrix first.")
        return

    validation = validate_normalization_input(processed_counts_df, sample_columns)
    summary = validation["summary"]
    st.markdown("### Matrix summary")
    render_info_line("Number of genes", f"{summary.get('genes', 0):,}")
    render_info_line("Number of samples", f"{summary.get('samples', 0):,}")
    render_info_line("Duplicated gene IDs count", f"{summary.get('duplicated_gene_ids', 0):,}")
    for warning in validation["warnings"]:
        st.warning(warning)
    for error in validation["errors"]:
        st.error(error)
    if not validation["valid"]:
        return

    current_signature = get_normalization_input_signature()
    stored_signature = st.session_state.get("normalization_input_signature")
    if stored_signature and stored_signature != current_signature:
        clear_normalization_state()
        stored_signature = None

    run_status = st.session_state.get("normalization_run_status")
    results = st.session_state.get("normalization_results")
    has_current_results = bool(results) and st.session_state.get("normalization_input_signature") == current_signature
    failed_current_input = (
        isinstance(run_status, dict)
        and run_status.get("success") is False
        and st.session_state.get("normalization_input_signature") == current_signature
    )

    if not has_current_results and not failed_current_input:
        with st.spinner("Running DESeq2 and edgeR normalization..."):
            run_normalization_workflow()
        run_status = st.session_state.get("normalization_run_status")
        results = st.session_state.get("normalization_results")
        has_current_results = bool(results) and st.session_state.get("normalization_input_signature") == current_signature
        failed_current_input = (
            isinstance(run_status, dict)
            and run_status.get("success") is False
            and st.session_state.get("normalization_input_signature") == current_signature
        )

    if run_status and not run_status.get("success"):
        st.error(run_status.get("stderr", "Normalization failed."))
        with st.expander("Rscript stderr/stdout"):
            st.code(run_status.get("stderr", ""), language="text")
            if run_status.get("stdout"):
                st.code(run_status.get("stdout", ""), language="text")
        if failed_current_input and st.button("Retry normalization", type="primary"):
            clear_normalization_state()
            with st.spinner("Running DESeq2 and edgeR normalization..."):
                run_normalization_workflow()
            st.rerun()
        return

    if not has_current_results:
        return
    st.success("Normalization completed.")

    output_dir = Path(results["output_dir"])
    for label, table in results.get("factor_tables", {}).items():
        st.markdown(f"### {label}")
        st.dataframe(table, use_container_width=True, hide_index=True)
    render_normalization_report(st.session_state.get("normalization_report"))

    st.subheader("Normalized matrix")
    matrix_labels = list(NORMALIZATION_OUTPUTS.keys())
    selected_matrix = st.session_state.get("normalization_selected_matrix", "Raw counts")
    if selected_matrix not in matrix_labels:
        selected_matrix = "Raw counts"
        st.session_state["normalization_selected_matrix"] = selected_matrix
    selected_label = st.selectbox(
        "Matrix",
        matrix_labels,
        index=matrix_labels.index(selected_matrix),
        key="normalization_selected_matrix",
        label_visibility="collapsed",
    )
    try:
        table_df = load_normalization_table(output_dir, selected_label)
    except FileNotFoundError:
        st.error(f"{selected_label} output is missing.")
        return
    render_matrix_table_viewer(
        table_df,
        selected_label,
        NORMALIZATION_OUTPUTS[selected_label],
        "normalization",
    )


def reset_qc_pca_settings() -> None:
    """Reset PCA plot settings and related widget state."""
    st.session_state["qc_pca_settings"] = default_qc_pca_settings()
    st.session_state["qc_pca_reset_nonce"] = int(st.session_state.get("qc_pca_reset_nonce", 0)) + 1
    clear_qc_export_cache("pca_plot")
    st.rerun()


def reset_qc_corr_settings() -> None:
    """Reset sample-correlation plot settings and related widget state."""
    st.session_state["qc_corr_settings"] = default_qc_corr_settings()
    st.session_state["qc_corr_reset_nonce"] = int(st.session_state.get("qc_corr_reset_nonce", 0)) + 1
    clear_qc_export_cache("sample_correlation")
    st.rerun()


def _select_expression_matrix(settings: dict[str, Any], matrices: dict[str, dict[str, Any]], key: str) -> str | None:
    """Render a normalized expression matrix selector with DESeq2 VST preferred."""
    if not matrices:
        return None
    labels = list(matrices.keys())
    selected = settings.get("matrix", "DESeq2 VST")
    if selected not in labels:
        selected = "DESeq2 VST" if "DESeq2 VST" in labels else labels[0]
        settings["matrix"] = selected
    selected = st.selectbox(
        "Matrix",
        labels,
        index=labels.index(selected),
        key=key,
    )
    settings["matrix"] = selected
    return selected


def render_qc_pca_section(sample_columns: list[str]) -> None:
    """Render PCA controls and plot."""
    st.markdown("### PCA")
    st.write("PCA is always computed on individual samples. QC assignment groups only affect coloring.")
    if np is None:
        st.error("NumPy is required for PCA. Install with: pip install numpy")
        return
    matrices = get_available_expression_matrices()
    if not matrices:
        st.warning("Please complete normalization first. DESeq2 VST is required for PCA and sample correlation.")
        return
    settings = st.session_state.setdefault("qc_pca_settings", default_qc_pca_settings())
    nonce = get_qc_reset_nonce("pca_plot")
    grouping_sets = st.session_state.get("qc_grouping_sets", {})
    grouping_names = list(grouping_sets.keys())

    control_cols = st.columns([1.35, 1.25, 1.65, 1.0, 0.7, 0.55, 0.55, 1.6], gap="small")
    with control_cols[0]:
        selected_matrix = _select_expression_matrix(settings, matrices, f"pca_plot*matrix*{nonce}")
    with control_cols[1]:
        settings["plot_by"] = st.selectbox(
            "Plot by",
            ["Sample name", "QC assignment group"],
            index=["Sample name", "QC assignment group"].index(settings.get("plot_by", "Sample name")),
            key=f"pca_plot*plot_by*{nonce}",
        )
    plot_by_group = settings["plot_by"] == "QC assignment group"
    with control_cols[2]:
        if plot_by_group and grouping_names:
            selected_grouping = settings.get("grouping_set") if settings.get("grouping_set") in grouping_names else grouping_names[0]
            settings["grouping_set"] = st.selectbox(
                "QC assignment set",
                grouping_names,
                index=grouping_names.index(selected_grouping),
                key=f"pca_plot*grouping_set*{nonce}",
            )
        else:
            st.selectbox(
                "QC assignment set",
                ["No QC assignment selected"],
                disabled=True,
                key=f"pca_plot*grouping_set_disabled*{nonce}",
            )
            settings["grouping_set"] = None

    if selected_matrix is None:
        return
    pca_result = compute_cached_pca(selected_matrix, matrices[selected_matrix], sample_columns)
    if pca_result.get("error"):
        st.warning(pca_result["error"])
        return
    component_pairs = pca_result.get("component_pairs", [])
    if settings.get("component_pair") not in component_pairs:
        settings["component_pair"] = component_pairs[0] if component_pairs else "PC1 vs PC2"
    with control_cols[3]:
        settings["component_pair"] = st.selectbox(
            "Component pair",
            component_pairs,
            index=component_pairs.index(settings["component_pair"]),
            key=f"pca_plot*component_pair*{nonce}",
            disabled=not component_pairs,
        )
    with control_cols[4]:
        align_button_with_input()
        if st.button("Reset", key=f"pca_plot*reset*{nonce}"):
            reset_qc_pca_settings()

    if plot_by_group and not grouping_names:
        st.warning("Please create and save a QC grouping set first.")

    with st.expander("Advanced settings"):
        adv_cols = st.columns(5)
        with adv_cols[0]:
            settings["width"] = st.slider("Plot width", 700, 1440, int(settings.get("width", 980)), 20, key=f"pca_plot*width*{nonce}")
        with adv_cols[1]:
            settings["height"] = st.slider("Plot height", 420, 900, int(settings.get("height", 560)), 20, key=f"pca_plot*height*{nonce}")
        with adv_cols[2]:
            settings["axis_label_font_size"] = st.slider(
                "Axis label font size",
                10,
                28,
                int(settings.get("axis_label_font_size", 14)),
                1,
                key=f"pca_plot*axis_label_font_size*{nonce}",
            )
        with adv_cols[3]:
            settings["point_size"] = st.slider("Point size", 5, 24, int(settings.get("point_size", 10)), 1, key=f"pca_plot*point_size*{nonce}")
        with adv_cols[4]:
            settings["show_sample_labels"] = st.checkbox(
                "Show sample labels",
                value=bool(settings.get("show_sample_labels", False)),
                key=f"pca_plot*show_sample_labels*{nonce}",
            )

    fig, entity_keys, entity_labels, _ = make_qc_pca_plot(
        pca_result,
        sample_columns,
        settings,
        settings.get("grouping_set"),
        int(settings.get("width", 980)),
        int(settings.get("height", 560)),
    )
    with st.expander("Color settings"):
        render_qc_color_settings("pca_plot", entity_keys, entity_labels)
    export_signature = hashlib.sha256(
        json.dumps(
            {
                "plot": "pca_plot",
                "matrix": selected_matrix,
                "matrix_signature": matrices[selected_matrix].get("signature"),
                "settings": settings,
                "grouping": st.session_state.get("qc_grouping_sets", {}).get(settings.get("grouping_set"), {}),
            },
            sort_keys=True,
            default=str,
        ).encode("utf-8")
    ).hexdigest()[:8]
    render_qc_export_buttons(fig, "pca_plot", f"pca_plot_{APP_VERSION}_{export_signature}", export_signature, control_cols[5], control_cols[6], nonce)
    st.caption(f"Matrix: {selected_matrix}. PCA used `{pca_result.get('effective_genes', 0):,}` top variable non-constant genes.")
    st.plotly_chart(fig, use_container_width=False, key="pca_plotly_chart")


def render_qc_sample_correlation_section(sample_columns: list[str]) -> None:
    """Render sample-correlation controls and heatmap."""
    st.markdown("### Sample Correlation")
    st.write("Sample-by-sample Pearson correlation heatmap using top variable genes.")
    if np is None:
        st.error("NumPy is required for sample correlation. Install with: pip install numpy")
        return
    matrices = get_available_expression_matrices()
    if not matrices:
        st.warning("Please complete normalization first. DESeq2 VST is required for PCA and sample correlation.")
        return
    settings = st.session_state.setdefault("qc_corr_settings", default_qc_corr_settings())
    nonce = get_qc_reset_nonce("sample_correlation")
    grouping_sets = st.session_state.get("qc_grouping_sets", {})
    grouping_names = list(grouping_sets.keys())

    control_cols = st.columns([1.35, 1.25, 1.15, 1.65, 0.7, 0.55, 0.55, 1.6], gap="small")
    with control_cols[0]:
        selected_matrix = _select_expression_matrix(settings, matrices, f"sample_correlation*matrix*{nonce}")
    with control_cols[1]:
        previous_plot_by = settings.get("plot_by", "Sample name")
        settings["plot_by"] = st.selectbox(
            "Plot by",
            ["Sample name", "QC assignment group"],
            index=["Sample name", "QC assignment group"].index(settings.get("plot_by", "Sample name")),
            key=f"sample_correlation*plot_by*{nonce}",
        )
    with control_cols[2]:
        top_options = ["500", "1000", "2000", "5000", "ALL"]
        current_top = settings.get("top_variable_genes", "1000")
        if current_top not in top_options:
            current_top = "1000"
        settings["top_variable_genes"] = st.selectbox(
            "Top variable genes",
            top_options,
            index=top_options.index(current_top),
            key=f"sample_correlation*top_variable_genes*{nonce}",
        )
    plot_by_group = settings["plot_by"] == "QC assignment group"
    if previous_plot_by != settings["plot_by"] and plot_by_group:
        settings["show_group_annotation"] = True
    with control_cols[3]:
        if plot_by_group and grouping_names:
            selected_grouping = settings.get("grouping_set") if settings.get("grouping_set") in grouping_names else grouping_names[0]
            settings["grouping_set"] = st.selectbox(
                "QC assignment set",
                grouping_names,
                index=grouping_names.index(selected_grouping),
                key=f"sample_correlation*grouping_set*{nonce}",
            )
        else:
            st.selectbox(
                "QC assignment set",
                ["No QC assignment selected"],
                disabled=True,
                key=f"sample_correlation*grouping_set_disabled*{nonce}",
            )
            settings["grouping_set"] = None
    with control_cols[4]:
        align_button_with_input()
        if st.button("Reset", key=f"sample_correlation*reset*{nonce}"):
            reset_qc_corr_settings()
    if selected_matrix is None:
        return
    if plot_by_group and not grouping_names:
        st.warning("Please create and save a QC grouping set first.")

    with st.expander("Advanced settings"):
        slider_cols = st.columns(4)
        with slider_cols[0]:
            settings["plot_size"] = st.slider("Plot size", 640, 1200, int(settings.get("plot_size", 860)), 20, key=f"sample_correlation*plot_size*{nonce}")
        with slider_cols[1]:
            settings["label_size"] = st.slider("Label size", 9, 18, int(settings.get("label_size", 12)), 1, key=f"sample_correlation*label_size*{nonce}")
        with slider_cols[2]:
            settings["colorbar_thickness"] = st.slider(
                "Colorbar thickness",
                8,
                34,
                int(settings.get("colorbar_thickness", 16)),
                1,
                key=f"sample_correlation*colorbar_thickness*{nonce}",
            )
        with slider_cols[3]:
            settings["colorbar_length"] = st.slider(
                "Colorbar length",
                0.45,
                1.0,
                float(settings.get("colorbar_length", 0.72)),
                0.05,
                key=f"sample_correlation*colorbar_length*{nonce}",
            )
        option_cols = st.columns(3)
        with option_cols[0]:
            settings["show_correlation_values"] = st.checkbox(
                "Show correlation values",
                value=bool(settings.get("show_correlation_values", False)),
                key=f"sample_correlation*show_values*{nonce}",
            )
        with option_cols[1]:
            settings["show_group_annotation"] = st.checkbox(
                "Show group annotation",
                value=bool(settings.get("show_group_annotation", True)) and plot_by_group,
                disabled=not plot_by_group,
                key=f"sample_correlation*show_annotation*{nonce}",
            )
        with option_cols[2]:
            angle_options = [0, 30, 45, 60, 90]
            settings["x_axis_angle"] = st.selectbox(
                "X-axis angle",
                angle_options,
                index=angle_options.index(int(settings.get("x_axis_angle", 45))),
                key=f"sample_correlation*x_axis_angle*{nonce}",
            )

    corr_result = compute_cached_sample_correlation(
        selected_matrix,
        matrices[selected_matrix],
        sample_columns,
        settings.get("top_variable_genes", "1000"),
    )
    if corr_result.get("error"):
        st.warning(corr_result["error"])
        return
    fig, entity_keys, entity_labels = make_sample_correlation_plot(
        corr_result,
        sample_columns,
        settings,
        settings.get("grouping_set"),
    )
    with st.expander("Color settings"):
        if entity_keys:
            render_qc_color_settings("sample_correlation", entity_keys, entity_labels)
        else:
            st.caption("Group annotation colors are available when plotting by QC assignment group.")
    export_signature = hashlib.sha256(
        json.dumps(
            {
                "plot": "sample_correlation",
                "matrix": selected_matrix,
                "matrix_signature": matrices[selected_matrix].get("signature"),
                "settings": settings,
                "grouping": st.session_state.get("qc_grouping_sets", {}).get(settings.get("grouping_set"), {}),
            },
            sort_keys=True,
            default=str,
        ).encode("utf-8")
    ).hexdigest()[:8]
    render_qc_export_buttons(
        fig,
        "sample_correlation",
        f"sample_correlation_{APP_VERSION}_{export_signature}",
        export_signature,
        control_cols[5],
        control_cols[6],
        nonce,
    )
    st.caption(
        f"Matrix: {selected_matrix}. Correlation used `{corr_result.get('effective_genes', 0):,}` non-constant genes."
    )
    st.plotly_chart(fig, use_container_width=False, key="sample_correlation_plotly_chart")


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

    qc_views = [
        "QC Summary & Grouping",
        "Library Size",
        "Detected Genes",
        "Zero-count Fraction",
        "PCA",
        "Sample Correlation",
    ]
    stored_qc_view = st.session_state.get("qc_active_view", "QC Summary & Grouping")
    if stored_qc_view == "QC Summary":
        stored_qc_view = "QC Summary & Grouping"
        st.session_state["qc_active_view"] = stored_qc_view
    active_view = st.radio(
        "QC view",
        qc_views,
        index=qc_views.index(stored_qc_view) if stored_qc_view in qc_views else 0,
        horizontal=True,
        key="qc_active_view",
        label_visibility="collapsed",
    )

    if active_view == "QC Summary & Grouping":
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
    elif active_view == "Library Size":
        render_qc_barplot_section("library_size", sample_qc_df)
    elif active_view == "Detected Genes":
        render_qc_barplot_section("detected_genes", sample_qc_df)
    elif active_view == "Zero-count Fraction":
        render_qc_barplot_section("zero_fraction", sample_qc_df)
    elif active_view == "PCA":
        render_qc_pca_section(sample_columns)
    elif active_view == "Sample Correlation":
        render_qc_sample_correlation_section(sample_columns)


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
            white-space: nowrap !important;
            text-align: center !important;
            justify-content: center !important;
            display: inline-flex !important;
            align-items: center !important;
            width: auto !important;
            min-width: fit-content !important;
            max-width: none !important;
            overflow: visible !important;
            text-overflow: clip !important;
            font-size: clamp(0.82rem, 0.8vw, 0.95rem) !important;
            padding: 0.45rem 0.9rem !important;
            line-height: 1.2 !important;
        }

        div[data-testid="stButton"] button p,
        div[data-testid="stDownloadButton"] button p,
        div[data-testid="stButton"] button span,
        div[data-testid="stDownloadButton"] button span {
            white-space: nowrap !important;
            overflow: visible !important;
            text-overflow: clip !important;
            margin: 0 !important;
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

    upload_tab, qc_tab, norm_tab = st.tabs(["Upload Count Matrix", "Quality Control", "Normalization"])
    with upload_tab:
        render_upload_count_matrix_tab()
    with qc_tab:
        render_quality_control_tab()
    with norm_tab:
        render_normalization_tab()


if __name__ == "__main__":
    main()
