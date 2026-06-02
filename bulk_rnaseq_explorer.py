"""
Bulk RNA-seq Explorer
Version: bulk_rnaseq_explorer_v2_0

Scope for v2.0:
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
- Add a Normalization workflow driven by DESeq2 and edgeR through Rscript.

To reduce Streamlit toolbar/menu visibility, users may create `.streamlit/config.toml` with:

[client]
toolbarMode = "minimal"

Required for QC plots:
pip install plotly kaleido

Optional for faster gene-map cache:
pip install pyarrow

Explicitly out of scope:
- DEG analysis, PCA, sample correlation, heatmap, volcano plot,
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


APP_VERSION = "bulk_rnaseq_explorer_v2_0"

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
        "axis_label_font_size": 14,
        "axis_title_font_size": 14,
        "colors": {},
    }


def sync_qc_axis_label_font_size(settings: dict[str, Any]) -> int:
    """Keep v1.11 axis-label sizing compatible with v1.10 session state."""
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
    nonce = st.session_state.setdefault("qc_plot_reset_nonce", {}).setdefault(plot_id, 0)

    control_cols = st.columns([1.45, 1.9, 1.15, 0.75, 0.55, 0.55, 2.2], gap="small")
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
    render_qc_export_buttons(fig, plot_id, filename_base, control_cols[4], control_cols[5], nonce)
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
    nonce: int,
) -> None:
    """Render direct PNG/SVG download buttons for one QC plot."""
    png_key = f"{plot_id}:{filename_base}:png"
    svg_key = f"{plot_id}:{filename_base}:svg"
    png_bytes = get_plotly_export_bytes(fig, "png", png_key)
    svg_bytes = get_plotly_export_bytes(fig, "svg", svg_key)
    with png_column:
        align_button_with_input()
        st.download_button(
            "PNG",
            data=png_bytes or b"",
            file_name=f"{filename_base}.png",
            mime="image/png",
            disabled=png_bytes is None,
            key=f"{plot_id}*download*png*{filename_base}*{nonce}",
            help="Download PNG",
        )
    with svg_column:
        align_button_with_input()
        st.download_button(
            "SVG",
            data=svg_bytes or b"",
            file_name=f"{filename_base}.svg",
            mime="image/svg+xml",
            disabled=svg_bytes is None,
            key=f"{plot_id}*download*svg*{filename_base}*{nonce}",
            help="Download SVG",
        )
    if png_bytes is None or svg_bytes is None:
        st.warning("Static image export requires kaleido. Install with: pip install kaleido")


def validate_normalization_input(counts_df: pd.DataFrame | None, sample_columns: list[str]) -> dict[str, Any]:
    """Validate the processed raw count matrix used by normalization."""
    result: dict[str, Any] = {
        "valid": False,
        "errors": [],
        "warnings": [],
        "summary": {},
        "total_counts_df": pd.DataFrame(columns=["Sample", "Total counts"]),
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
    result["total_counts_df"] = pd.DataFrame(
        {"Sample": list(total_counts.index), "Total counts": [float(value) for value in total_counts.values]}
    )
    result["numeric_counts"] = numeric_counts
    return result


def get_normalization_input_signature() -> str:
    """Return a signature for the current processed matrix and selected sample columns."""
    payload = {
        "counts": st.session_state.get("counts_file_signature"),
        "gene_id_column": st.session_state.get("gene_id_column"),
        "samples": st.session_state.get("sample_columns", []),
        "processed_shape": getattr(st.session_state.get("processed_counts_df"), "shape", None),
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


def get_normalization_output_dimensions(output_dir: Path) -> pd.DataFrame:
    """Return generated/missing status and dimensions for normalization matrix outputs."""
    rows: list[dict[str, Any]] = []
    for label, filename in NORMALIZATION_OUTPUTS.items():
        path = output_dir / filename
        if not path.exists():
            rows.append({"Matrix": label, "File": filename, "Genes": None, "Samples": None, "Status": "Missing"})
            continue
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.reader(handle)
            header = next(reader, [])
            genes = sum(1 for _ in reader)
        rows.append(
            {
                "Matrix": label,
                "File": filename,
                "Genes": genes,
                "Samples": max(len(header) - 1, 0),
                "Status": "Generated",
            }
        )
    return pd.DataFrame(rows)


def run_normalization_workflow() -> None:
    """Prepare input, run R normalization, and store output metadata in session state."""
    processed_counts_df = st.session_state.get("processed_counts_df")
    sample_columns = st.session_state.get("sample_columns", [])
    validation = validate_normalization_input(processed_counts_df, sample_columns)
    if not validation["valid"]:
        st.session_state["normalization_run_status"] = {
            "success": False,
            "stderr": "\n".join(validation["errors"]),
            "stdout": "",
            "returncode": None,
            "command": [],
        }
        return

    output_dir = create_normalization_output_dir()
    input_path = prepare_normalization_input(processed_counts_df, sample_columns, output_dir)
    run_status = run_r_normalization(input_path, output_dir)
    st.session_state["normalization_run_status"] = run_status
    st.session_state["normalization_output_dir"] = str(output_dir)
    st.session_state["normalization_input_signature"] = get_normalization_input_signature()
    st.session_state["normalization_tables"] = {}
    if run_status["success"]:
        report = load_normalization_report(output_dir)
        dimensions = get_normalization_output_dimensions(output_dir)
        factor_tables = {}
        for label, filename in NORMALIZATION_FACTOR_FILES.items():
            path = output_dir / filename
            factor_tables[label] = pd.read_csv(path) if path.exists() else pd.DataFrame()
        st.session_state["normalization_report"] = report
        st.session_state["normalization_results"] = {
            "output_dir": str(output_dir),
            "dimensions": dimensions,
            "factor_tables": factor_tables,
        }
        st.session_state["normalization_selected_matrix"] = "Raw counts"
        st.session_state["normalization_table_page"] = 1
    else:
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
    search_key = f"{key_prefix}_search"
    rows_key = f"{key_prefix}_rows_per_page"
    page_key = f"{key_prefix}_page"
    previous_search_key = f"{key_prefix}_previous_search"
    previous_rows_key = f"{key_prefix}_previous_rows"
    st.session_state.setdefault(search_key, "")
    st.session_state.setdefault(rows_key, 25)
    st.session_state.setdefault(page_key, 1)

    control_cols = st.columns([2.4, 1.0, 0.6, 0.6, 0.8, 2.0], gap="small")
    with control_cols[0]:
        search_query = st.text_input(
            "Search genes",
            key=search_key,
            placeholder="Type at least 3 characters...",
        )
    with control_cols[1]:
        rows_per_page = st.selectbox(
            "Show rows",
            [10, 25, 50, 100, 250],
            index=[10, 25, 50, 100, 250].index(int(st.session_state.get(rows_key, 25))),
            key=rows_key,
        )

    if (
        st.session_state.get(previous_search_key) != search_query
        or st.session_state.get(previous_rows_key) != rows_per_page
    ):
        st.session_state[page_key] = 1
        st.session_state[previous_search_key] = search_query
        st.session_state[previous_rows_key] = rows_per_page

    gene_series = table_df["Gene"].astype(str) if "Gene" in table_df.columns else pd.Series([], dtype=str)
    search_active = len(search_query.strip()) >= 3
    if search_active:
        mask = gene_series.str.contains(re.escape(search_query.strip()), case=False, na=False)
        display_df = table_df.loc[mask].copy()
        suggestions = display_df["Gene"].astype(str).head(10).tolist() if "Gene" in display_df.columns else []
        if suggestions:
            st.caption("Suggestions: " + " / ".join(f"`{gene}`" for gene in suggestions))
        else:
            st.caption("No matched genes found.")
    else:
        display_df = table_df

    total_rows = int(display_df.shape[0])
    rows_per_page = int(rows_per_page)
    total_pages = max((total_rows + rows_per_page - 1) // rows_per_page, 1)
    st.session_state[page_key] = min(max(int(st.session_state.get(page_key, 1)), 1), total_pages)
    page = int(st.session_state[page_key])
    start = (page - 1) * rows_per_page
    end = min(start + rows_per_page, total_rows)
    page_df = display_df.iloc[start:end].copy()

    st.write(f"Matrix shape: `{table_df.shape[0]:,}` genes x `{max(table_df.shape[1] - 1, 0):,}` samples")
    if search_active:
        st.write(f"Search matched `{total_rows:,}` genes")
    st.write(f"Showing rows `{start + 1 if total_rows else 0:,}`-`{end:,}` of `{total_rows:,}`")

    with control_cols[2]:
        align_button_with_input()
        if st.button("Previous", key=f"{key_prefix}_previous", disabled=page <= 1):
            st.session_state[page_key] = max(page - 1, 1)
            st.rerun()
    with control_cols[3]:
        align_button_with_input()
        if st.button("Next", key=f"{key_prefix}_next", disabled=page >= total_pages):
            st.session_state[page_key] = min(page + 1, total_pages)
            st.rerun()
    with control_cols[4]:
        align_button_with_input()
        st.download_button(
            "CSV",
            data=display_df.to_csv(index=False).encode("utf-8"),
            file_name=download_filename,
            mime="text/csv",
            key=f"{key_prefix}_csv",
        )
    st.caption(f"Page {page} of {total_pages}")
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
        ("Removed all-zero genes", "all_zero_genes"),
        ("Samples", "samples"),
        ("R version", "r_version"),
        ("DESeq2 version", "deseq2_version"),
        ("edgeR version", "edger_version"),
        ("Timestamp", "timestamp"),
        ("Rounding applied", "rounding_applied"),
        ("Max deviation from integer", "max_deviation_from_integer"),
        ("Non-integer value fraction", "non_integer_value_fraction"),
    ]
    for label, key in report_fields:
        render_info_line(label, report.get(key, "Not available"))
    with st.expander("View normalization report JSON"):
        st.json(report)


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
    render_info_line("All-zero genes count", f"{summary.get('all_zero_genes', 0):,}")
    render_info_line("Non-integer value fraction", f"{summary.get('non_integer_value_fraction', 0):.6g}")
    render_info_line("Max deviation from integer", f"{summary.get('max_deviation_from_integer', 0):.6g}")
    st.markdown("### Total counts per sample")
    st.dataframe(validation["total_counts_df"], use_container_width=True, hide_index=True)
    for warning in validation["warnings"]:
        st.warning(warning)
    for error in validation["errors"]:
        st.error(error)

    if st.button("Run normalization", disabled=not validation["valid"], type="primary"):
        with st.spinner("Running DESeq2 and edgeR normalization..."):
            run_normalization_workflow()
        if st.session_state.get("normalization_run_status", {}).get("success"):
            st.success("Normalization completed.")
        else:
            st.error("Normalization failed.")

    run_status = st.session_state.get("normalization_run_status")
    if run_status and not run_status.get("success"):
        st.error(run_status.get("stderr", "Normalization failed."))
        with st.expander("Rscript stderr/stdout"):
            st.code(run_status.get("stderr", ""), language="text")
            if run_status.get("stdout"):
                st.code(run_status.get("stdout", ""), language="text")

    results = st.session_state.get("normalization_results")
    if not results:
        return

    output_dir = Path(results["output_dir"])
    st.markdown("### Results summary")
    render_info_line("Output directory", output_dir)
    dimensions = results.get("dimensions", pd.DataFrame())
    if isinstance(dimensions, pd.DataFrame) and not dimensions.empty:
        st.dataframe(dimensions, use_container_width=True, hide_index=True)
    for label, table in results.get("factor_tables", {}).items():
        st.markdown(f"### {label}")
        st.dataframe(table, use_container_width=True, hide_index=True)
    render_normalization_report(st.session_state.get("normalization_report"))

    selected_label = st.selectbox(
        "Normalized matrix",
        list(NORMALIZATION_OUTPUTS.keys()),
        index=list(NORMALIZATION_OUTPUTS.keys()).index(st.session_state.get("normalization_selected_matrix", "Raw counts")),
        key="normalization_selected_matrix",
    )
    table_df = load_normalization_table(output_dir, selected_label)
    render_matrix_table_viewer(
        table_df,
        selected_label,
        NORMALIZATION_OUTPUTS[selected_label],
        f"normalization_{re.sub(r'[^a-zA-Z0-9]+', '_', selected_label).strip('_').lower()}",
    )


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
