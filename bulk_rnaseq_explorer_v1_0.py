"""
BL Bulk RNA-seq Explorer
Version: BL_Bulk-seq_v1.0

Scope for this MVP:
- Streamlit application skeleton for a browser-to-Python migration.
- Project setup/resource detection.
- Raw count matrix, sample metadata, and optional gene map upload.
- Basic input validation and project state management.

Explicitly out of scope for this version:
- DESeq2, Rscript calls, DEG analysis, QC plots, PCA, heatmap, volcano plot,
  pathway database parsing, GSEA, ORA, cloud storage, login, Duke DCC, SLURM.
"""

from __future__ import annotations

import csv
from io import StringIO
from pathlib import Path
from typing import Any

import pandas as pd
import streamlit as st


APP_VERSION = "BL_Bulk-seq_v1.0"
PREFERRED_SAMPLE_COLUMNS = [
    "sample",
    "Sample",
    "sample_id",
    "SampleID",
    "sample_name",
    "SampleName",
    "sampleid",
    "name",
]


def init_session_state() -> None:
    """Initialize persistent project state."""
    defaults: dict[str, Any] = {
        "app_version": APP_VERSION,
        "counts_df": None,
        "metadata_df": None,
        "gene_map_df": None,
        "gene_id_column": None,
        "sample_name_column": None,
        "group_column": None,
        "validation_status": "Not run",
        "validation_report": None,
        "downstream_invalidated": False,
        "qc_results": None,
        "deg_results": None,
        "pathway_results": None,
        "plots": None,
        "counts_file_name": None,
        "metadata_file_name": None,
        "gene_map_file_name": None,
    }
    for key, value in defaults.items():
        if key not in st.session_state:
            st.session_state[key] = value


def invalidate_downstream_results() -> None:
    """Clear placeholders for future analyses when inputs or key settings change."""
    st.session_state["downstream_invalidated"] = True
    st.session_state["qc_results"] = None
    st.session_state["deg_results"] = None
    st.session_state["pathway_results"] = None
    st.session_state["plots"] = None
    st.session_state["validation_status"] = "Not run"
    st.session_state["validation_report"] = None


def detect_local_resources() -> pd.DataFrame:
    """Detect local historical/reference resources without parsing them."""
    cwd = Path.cwd()
    resource_specs = [
        ("Historical browser prototype", "index_v5.4.9.html", "File"),
        ("Top-level mouse gene map JS", "mouse_ensembl_to_symbol.js", "File"),
        ("Gene sets JS bundle", "gene_sets_2026_v1.js", "File"),
        ("GMT database directory", "database_raw", "Directory"),
        ("Source mapping directory", "source_mapping", "Directory"),
        ("Source mapping table", "source_mapping/mouse_ensembl_to_symbol", "File"),
        ("Source GTF archive", "source_mapping/Mus_musculus.GRCm39.115.gtf.gz", "File"),
    ]

    rows = []
    for resource, relative_path, resource_type in resource_specs:
        path = cwd / relative_path
        rows.append(
            {
                "Resource": resource,
                "Path": str(path),
                "Detected": path.exists(),
                "Type": resource_type,
            }
        )
    return pd.DataFrame(rows)


def read_table_uploaded_file(uploaded_file) -> pd.DataFrame:
    """Read an uploaded delimited table using suffix-based separator detection."""
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
        try:
            return pd.read_csv(uploaded_file, sep=None, engine="python")
        except Exception:
            uploaded_file.seek(0)
            text = uploaded_file.read().decode("utf-8", errors="replace")
            lines = text.splitlines()
            return pd.DataFrame({"raw_text": lines})


def infer_sample_column(metadata_df: pd.DataFrame | None) -> str | None:
    """Guess the metadata sample-name column from common names, then first column."""
    if metadata_df is None or metadata_df.empty or len(metadata_df.columns) == 0:
        return None
    for candidate in PREFERRED_SAMPLE_COLUMNS:
        if candidate in metadata_df.columns:
            return candidate
    return str(metadata_df.columns[0])


def _make_check(
    section: str,
    check: str,
    level: str,
    message: str,
    value: Any = "",
) -> dict[str, Any]:
    return {
        "Section": section,
        "Check": check,
        "Level": level,
        "Message": message,
        "Value": value,
    }


def _count_severity(checks: list[dict[str, Any]], level: str) -> int:
    return sum(1 for item in checks if item["Level"] == level)


def validate_counts(counts_df: pd.DataFrame | None, gene_id_column: str | None) -> dict[str, Any]:
    """Validate raw count matrix structure and count values."""
    checks: list[dict[str, Any]] = []
    summary: dict[str, Any] = {
        "number_of_genes": 0,
        "number_of_samples": 0,
        "number_of_duplicated_gene_ids": 0,
        "number_of_all_zero_genes": 0,
        "number_of_genes_with_missing_values": 0,
        "number_of_non_integer_count_entries": 0,
        "sample_columns": [],
    }

    if counts_df is None:
        checks.append(_make_check("Counts", "Read count matrix", "Error", "Count matrix is not uploaded."))
        return {"checks": checks, "summary": summary}

    checks.append(_make_check("Counts", "Read count matrix", "Pass", "Count matrix was read successfully."))

    if gene_id_column not in counts_df.columns:
        checks.append(
            _make_check("Counts", "Gene ID column exists", "Error", "Selected gene ID column is not present.")
        )
        return {"checks": checks, "summary": summary}

    checks.append(_make_check("Counts", "Gene ID column exists", "Pass", "Selected gene ID column is present."))

    sample_columns = [column for column in counts_df.columns if column != gene_id_column]
    summary["number_of_genes"] = int(counts_df.shape[0])
    summary["number_of_samples"] = int(len(sample_columns))
    summary["sample_columns"] = [str(column) for column in sample_columns]

    if len(sample_columns) < 2:
        checks.append(
            _make_check(
                "Counts",
                "At least 2 sample columns",
                "Error",
                "Count matrix must contain at least 2 sample columns.",
                len(sample_columns),
            )
        )
    else:
        checks.append(
            _make_check("Counts", "At least 2 sample columns", "Pass", "At least 2 sample columns detected.", len(sample_columns))
        )

    gene_ids = counts_df[gene_id_column]
    missing_gene_ids = int(gene_ids.isna().sum() + (gene_ids.astype("string").str.strip() == "").sum())
    duplicated_gene_ids = int(gene_ids.duplicated(keep=False).sum())
    summary["number_of_duplicated_gene_ids"] = duplicated_gene_ids

    if missing_gene_ids > 0:
        checks.append(
            _make_check("Counts", "Missing gene IDs", "Error", "Gene ID column contains missing values.", missing_gene_ids)
        )
    else:
        checks.append(_make_check("Counts", "Missing gene IDs", "Pass", "No missing gene IDs detected.", 0))

    if duplicated_gene_ids > 0:
        checks.append(
            _make_check("Counts", "Duplicated gene IDs", "Warning", "Duplicated gene IDs detected.", duplicated_gene_ids)
        )
    else:
        checks.append(_make_check("Counts", "Duplicated gene IDs", "Pass", "No duplicated gene IDs detected.", 0))

    duplicated_sample_columns = int(pd.Index(sample_columns).duplicated(keep=False).sum())
    if duplicated_sample_columns > 0:
        checks.append(
            _make_check(
                "Counts",
                "Duplicated sample columns",
                "Error",
                "Duplicated sample column names detected.",
                duplicated_sample_columns,
            )
        )
    else:
        checks.append(_make_check("Counts", "Duplicated sample columns", "Pass", "No duplicated sample columns detected.", 0))

    count_values = counts_df[sample_columns].apply(pd.to_numeric, errors="coerce") if sample_columns else pd.DataFrame()
    non_numeric_columns = []
    for column in sample_columns:
        original_values = counts_df[column]
        original_as_text = original_values.astype("string").str.strip()
        non_missing_original = original_values.notna() & original_as_text.ne("")
        coerced_missing = count_values[column].isna()
        if (non_missing_original & coerced_missing).any():
            non_numeric_columns.append(str(column))
    if non_numeric_columns:
        checks.append(
            _make_check(
                "Counts",
                "Numeric sample columns",
                "Error",
                "One or more sample columns contain non-numeric values.",
                ", ".join(non_numeric_columns[:10]),
            )
        )
    else:
        checks.append(_make_check("Counts", "Numeric sample columns", "Pass", "All sample columns are numeric."))

    missing_entries = int(count_values.isna().sum().sum()) if not count_values.empty else 0
    genes_with_missing = int(count_values.isna().any(axis=1).sum()) if not count_values.empty else 0
    summary["number_of_genes_with_missing_values"] = genes_with_missing
    if missing_entries > 0:
        checks.append(
            _make_check("Counts", "Missing count values", "Warning", "Missing values detected in count matrix.", missing_entries)
        )
    else:
        checks.append(_make_check("Counts", "Missing count values", "Pass", "No missing count values detected.", 0))

    negative_entries = int((count_values < 0).sum().sum()) if not count_values.empty else 0
    if negative_entries > 0:
        checks.append(_make_check("Counts", "Negative counts", "Error", "Negative count values detected.", negative_entries))
    else:
        checks.append(_make_check("Counts", "Negative counts", "Pass", "No negative count values detected.", 0))

    valid_numeric_values = count_values.dropna()
    non_integer_entries = int(((valid_numeric_values % 1) != 0).sum().sum()) if not valid_numeric_values.empty else 0
    summary["number_of_non_integer_count_entries"] = non_integer_entries
    if non_integer_entries > 0:
        checks.append(
            _make_check(
                "Counts",
                "Integer counts",
                "Warning",
                "DESeq2 expects raw integer counts.",
                non_integer_entries,
            )
        )
    else:
        checks.append(_make_check("Counts", "Integer counts", "Pass", "All numeric count entries are integers.", 0))

    all_zero_genes = int((count_values.fillna(0).sum(axis=1) == 0).sum()) if not count_values.empty else 0
    summary["number_of_all_zero_genes"] = all_zero_genes
    if all_zero_genes > 0:
        checks.append(_make_check("Counts", "All-zero genes", "Warning", "All-zero genes detected.", all_zero_genes))
    else:
        checks.append(_make_check("Counts", "All-zero genes", "Pass", "No all-zero genes detected.", 0))

    return {"checks": checks, "summary": summary}


def validate_metadata(
    metadata_df: pd.DataFrame | None,
    sample_name_column: str | None,
    count_sample_names: list[str],
) -> dict[str, Any]:
    """Validate sample metadata and sample-name matching against count columns."""
    checks: list[dict[str, Any]] = []
    summary: dict[str, Any] = {
        "metadata_sample_number": 0,
        "matched_sample_number": 0,
        "samples_in_counts_but_missing_in_metadata": [],
        "samples_in_metadata_but_missing_in_counts": [],
    }

    if metadata_df is None:
        checks.append(_make_check("Metadata", "Read metadata", "Error", "Metadata is not uploaded."))
        return {"checks": checks, "summary": summary}

    checks.append(_make_check("Metadata", "Read metadata", "Pass", "Metadata was read successfully."))

    if metadata_df.shape[1] < 2:
        checks.append(
            _make_check("Metadata", "At least 2 columns", "Error", "Metadata must contain at least 2 columns.", metadata_df.shape[1])
        )
    else:
        checks.append(_make_check("Metadata", "At least 2 columns", "Pass", "At least 2 metadata columns detected.", metadata_df.shape[1]))

    if sample_name_column not in metadata_df.columns:
        checks.append(
            _make_check("Metadata", "Sample name column exists", "Error", "Selected sample name column is not present.")
        )
        return {"checks": checks, "summary": summary}

    checks.append(_make_check("Metadata", "Sample name column exists", "Pass", "Selected sample name column is present."))

    sample_series = metadata_df[sample_name_column].astype("string").str.strip()
    missing_samples = int(sample_series.isna().sum() + (sample_series == "").sum())
    duplicated_samples = int(sample_series.duplicated(keep=False).sum())
    metadata_samples = [str(value) for value in sample_series.dropna().tolist() if str(value) != ""]
    metadata_sample_set = set(metadata_samples)
    count_sample_set = set(map(str, count_sample_names))

    missing_in_metadata = sorted(count_sample_set - metadata_sample_set)
    missing_in_counts = sorted(metadata_sample_set - count_sample_set)
    matched_samples = sorted(count_sample_set & metadata_sample_set)

    summary["metadata_sample_number"] = int(len(metadata_samples))
    summary["matched_sample_number"] = int(len(matched_samples))
    summary["samples_in_counts_but_missing_in_metadata"] = missing_in_metadata
    summary["samples_in_metadata_but_missing_in_counts"] = missing_in_counts

    if missing_samples > 0:
        checks.append(
            _make_check("Metadata", "Missing sample names", "Error", "Sample name column contains missing values.", missing_samples)
        )
    else:
        checks.append(_make_check("Metadata", "Missing sample names", "Pass", "No missing sample names detected.", 0))

    if duplicated_samples > 0:
        checks.append(
            _make_check("Metadata", "Duplicated sample names", "Error", "Duplicated sample names detected.", duplicated_samples)
        )
    else:
        checks.append(_make_check("Metadata", "Duplicated sample names", "Pass", "No duplicated sample names detected.", 0))

    if count_sample_names and len(matched_samples) == 0:
        checks.append(
            _make_check(
                "Metadata",
                "Sample matching",
                "Error",
                "No count sample columns match metadata sample names.",
                "0 matched",
            )
        )
    elif missing_in_metadata or missing_in_counts:
        checks.append(
            _make_check(
                "Metadata",
                "Sample matching",
                "Warning",
                "Count and metadata sample names partially match.",
                f"{len(matched_samples)} matched",
            )
        )
    else:
        checks.append(
            _make_check(
                "Metadata",
                "Sample matching",
                "Pass",
                "Count sample columns and metadata sample names match.",
                f"{len(matched_samples)} matched",
            )
        )

    return {"checks": checks, "summary": summary}


def validate_group_column(metadata_df: pd.DataFrame | None, group_column: str | None) -> dict[str, Any]:
    """Validate selected experimental group column."""
    checks: list[dict[str, Any]] = []
    summary: dict[str, Any] = {"number_of_groups": 0, "group_sample_counts": {}}

    if metadata_df is None or not group_column:
        checks.append(_make_check("Group", "Group column selected", "Warning", "No group column selected."))
        return {"checks": checks, "summary": summary}

    if group_column not in metadata_df.columns:
        checks.append(_make_check("Group", "Group column exists", "Error", "Selected group column is not present."))
        return {"checks": checks, "summary": summary}

    checks.append(_make_check("Group", "Group column exists", "Pass", "Selected group column is present."))

    group_series = metadata_df[group_column].astype("string").str.strip()
    missing_groups = int(group_series.isna().sum() + (group_series == "").sum())
    if missing_groups > 0:
        checks.append(_make_check("Group", "Missing group values", "Warning", "Group column contains missing values.", missing_groups))
    else:
        checks.append(_make_check("Group", "Missing group values", "Pass", "No missing group values detected.", 0))

    valid_groups = group_series.dropna()
    valid_groups = valid_groups[valid_groups != ""]
    group_counts = valid_groups.value_counts().sort_index()
    summary["number_of_groups"] = int(group_counts.shape[0])
    summary["group_sample_counts"] = {str(group): int(count) for group, count in group_counts.items()}

    if group_counts.shape[0] < 2:
        checks.append(_make_check("Group", "At least 2 groups", "Error", "At least 2 groups are required.", group_counts.shape[0]))
    else:
        checks.append(_make_check("Group", "At least 2 groups", "Pass", "At least 2 groups detected.", group_counts.shape[0]))

    single_replicate_groups = [str(group) for group, count in group_counts.items() if count == 1]
    if single_replicate_groups:
        checks.append(
            _make_check(
                "Group",
                "Biological replicates",
                "Warning",
                "DESeq2 requires biological replicates for reliable inference.",
                ", ".join(single_replicate_groups),
            )
        )
    elif not group_counts.empty:
        checks.append(
            _make_check("Group", "Biological replicates", "Pass", "All groups contain at least 2 samples.")
        )

    return {"checks": checks, "summary": summary}


def run_full_validation() -> dict[str, Any]:
    """Run all available validation checks and store report in session state."""
    counts_df = st.session_state["counts_df"]
    metadata_df = st.session_state["metadata_df"]
    gene_id_column = st.session_state["gene_id_column"]
    sample_name_column = st.session_state["sample_name_column"]
    group_column = st.session_state["group_column"]

    if counts_df is None or metadata_df is None:
        st.session_state["validation_status"] = "Not run"
        report = {
            "checks": [],
            "counts_summary": {},
            "metadata_summary": {},
            "group_summary": {},
            "status": "Not run",
        }
        st.session_state["validation_report"] = report
        return report

    counts_report = validate_counts(counts_df, gene_id_column)
    count_sample_names = counts_report["summary"].get("sample_columns", [])
    metadata_report = validate_metadata(metadata_df, sample_name_column, count_sample_names)
    group_report = validate_group_column(metadata_df, group_column)

    all_checks = counts_report["checks"] + metadata_report["checks"] + group_report["checks"]
    if _count_severity(all_checks, "Error") > 0:
        status = "Failed"
    elif _count_severity(all_checks, "Warning") > 0:
        status = "Warning"
    else:
        status = "Passed"

    report = {
        "checks": all_checks,
        "counts_summary": counts_report["summary"],
        "metadata_summary": metadata_report["summary"],
        "group_summary": group_report["summary"],
        "status": status,
    }
    st.session_state["validation_status"] = status
    st.session_state["validation_report"] = report
    return report


def _status_badge(label: str, enabled: bool, locked: bool = False) -> None:
    if enabled:
        st.sidebar.success(label)
    elif locked:
        st.sidebar.caption(f"{label} - Coming soon / Locked")
    else:
        st.sidebar.info(label)


def render_sidebar() -> None:
    """Render workflow sidebar and project status."""
    st.sidebar.title("Workflow")
    _status_badge("1. Project Setup", True)
    _status_badge("2. Upload Data", True)
    _status_badge("3. Input Validation", True)
    _status_badge("4. QC Overview", False, locked=True)
    _status_badge("5. DEG Analysis", False, locked=True)
    _status_badge("6. Visualization", False, locked=True)
    _status_badge("7. Pathway Analysis", False, locked=True)
    _status_badge("8. Export", False, locked=True)

    resources = detect_local_resources()
    gene_map_detected = bool(
        resources.loc[
            resources["Resource"].isin(["Top-level mouse gene map JS", "Source mapping table"]),
            "Detected",
        ].any()
        or st.session_state["gene_map_df"] is not None
    )

    st.sidebar.divider()
    st.sidebar.subheader("Project status")
    st.sidebar.write(f"Counts uploaded: {'Yes' if st.session_state['counts_df'] is not None else 'No'}")
    st.sidebar.write(f"Metadata uploaded: {'Yes' if st.session_state['metadata_df'] is not None else 'No'}")
    st.sidebar.write(f"Gene map detected: {'Yes' if gene_map_detected else 'No'}")
    st.sidebar.write(f"Validation status: {st.session_state['validation_status']}")

    if st.session_state["downstream_invalidated"]:
        st.sidebar.warning("Downstream placeholders invalidated.")


def render_project_setup_tab() -> None:
    """Render project setup and local resource detection."""
    st.subheader("Project Setup")
    st.write(f"Current app version: `{st.session_state['app_version']}`")
    st.write(f"Current working directory: `{Path.cwd()}`")

    resources = detect_local_resources()
    display_resources = resources.copy()
    display_resources["Detected"] = display_resources["Detected"].map(lambda value: "Yes" if value else "No")
    st.dataframe(display_resources, use_container_width=True, hide_index=True)

    st.info(
        "v5.4.9 is used as historical reference only; this Python version will rebuild "
        "the pipeline using modular and reproducible logic."
    )


def _handle_uploaded_table(
    uploaded_file,
    state_key: str,
    file_name_key: str,
    success_label: str,
    invalidate: bool = True,
) -> None:
    if uploaded_file is None:
        return

    if st.session_state.get(file_name_key) == uploaded_file.name and st.session_state.get(state_key) is not None:
        return

    try:
        st.session_state[state_key] = read_table_uploaded_file(uploaded_file)
        st.session_state[file_name_key] = uploaded_file.name
        if invalidate:
            invalidate_downstream_results()
        st.success(success_label)
    except Exception as exc:
        st.session_state[state_key] = None
        st.session_state[file_name_key] = uploaded_file.name
        if invalidate:
            invalidate_downstream_results()
        st.error(f"Failed to read `{uploaded_file.name}`: {exc}")


def render_upload_data_tab() -> None:
    """Render upload controls, key column selectors, and previews."""
    st.subheader("Upload Data")

    counts_file = st.file_uploader("Upload raw count matrix", type=["csv", "tsv", "txt"], key="counts_uploader")
    _handle_uploaded_table(counts_file, "counts_df", "counts_file_name", "Count matrix uploaded.")

    counts_df = st.session_state["counts_df"]
    if counts_df is not None:
        columns = [str(column) for column in counts_df.columns]
        current_gene_col = st.session_state["gene_id_column"] if st.session_state["gene_id_column"] in columns else columns[0]
        selected_gene_col = st.selectbox(
            "Count matrix gene ID column",
            options=columns,
            index=columns.index(current_gene_col),
            key="gene_id_column_selector",
        )
        if selected_gene_col != st.session_state["gene_id_column"]:
            st.session_state["gene_id_column"] = selected_gene_col
            invalidate_downstream_results()

        sample_count = max(len(columns) - 1, 0)
        st.write(f"Detected shape: count matrix = `{counts_df.shape[0]}` genes x `{sample_count}` samples")
        st.dataframe(counts_df.head(5), use_container_width=True)

    st.divider()

    metadata_file = st.file_uploader("Upload sample metadata", type=["csv", "tsv", "txt"], key="metadata_uploader")
    _handle_uploaded_table(metadata_file, "metadata_df", "metadata_file_name", "Metadata uploaded.")

    metadata_df = st.session_state["metadata_df"]
    if metadata_df is not None:
        metadata_columns = [str(column) for column in metadata_df.columns]
        inferred_sample_column = infer_sample_column(metadata_df)
        current_sample_col = (
            st.session_state["sample_name_column"]
            if st.session_state["sample_name_column"] in metadata_columns
            else inferred_sample_column
        )
        selected_sample_col = st.selectbox(
            "Metadata sample name column",
            options=metadata_columns,
            index=metadata_columns.index(current_sample_col),
            key="sample_name_column_selector",
        )
        if selected_sample_col != st.session_state["sample_name_column"]:
            st.session_state["sample_name_column"] = selected_sample_col
            invalidate_downstream_results()

        group_options = ["None"] + metadata_columns
        current_group_col = st.session_state["group_column"] if st.session_state["group_column"] in metadata_columns else "None"
        selected_group_col = st.selectbox(
            "Metadata group column",
            options=group_options,
            index=group_options.index(current_group_col),
            key="group_column_selector",
        )
        normalized_group_col = None if selected_group_col == "None" else selected_group_col
        if normalized_group_col != st.session_state["group_column"]:
            st.session_state["group_column"] = normalized_group_col
            invalidate_downstream_results()

        st.write(f"Detected shape: metadata = `{metadata_df.shape[0]}` samples x `{metadata_df.shape[1]}` columns")
        st.dataframe(metadata_df.head(5), use_container_width=True)

    st.divider()

    gene_map_file = st.file_uploader(
        "Optional: upload gene mapping file",
        type=["csv", "tsv", "txt", "js"],
        key="gene_map_uploader",
    )
    _handle_uploaded_table(gene_map_file, "gene_map_df", "gene_map_file_name", "Gene map uploaded.", invalidate=False)

    resources = detect_local_resources()
    local_gene_map_detected = bool(
        resources.loc[
            resources["Resource"].isin(["Top-level mouse gene map JS", "Source mapping table"]),
            "Detected",
        ].any()
    )
    if st.session_state["gene_map_df"] is not None:
        st.success("Uploaded gene map detected.")
        st.dataframe(st.session_state["gene_map_df"].head(5), use_container_width=True)
    elif local_gene_map_detected:
        st.success("Local default gene map detected.")
    else:
        st.warning("No gene map detected. This is acceptable for Step 1.")


def _render_status_message(status: str) -> None:
    if status == "Passed":
        st.success("Validation passed.")
    elif status == "Warning":
        st.warning("Validation completed with warnings.")
    elif status == "Failed":
        st.error("Validation failed. Resolve errors before downstream analysis.")
    else:
        st.info("Upload both count matrix and metadata to run validation.")


def render_validation_tab() -> None:
    """Render validation controls and reports."""
    st.subheader("Input Validation")

    counts_ready = st.session_state["counts_df"] is not None
    metadata_ready = st.session_state["metadata_df"] is not None
    if not counts_ready or not metadata_ready:
        st.info("Upload both count matrix and metadata before running validation.")
        _render_status_message(st.session_state["validation_status"])
        return

    if st.button("Run validation", type="primary"):
        run_full_validation()

    report = st.session_state["validation_report"]
    _render_status_message(st.session_state["validation_status"])

    if report is None:
        st.info("Validation has not been run for the current inputs.")
        return

    checks_df = pd.DataFrame(report["checks"])
    st.markdown("#### Validation checks")
    st.dataframe(checks_df, use_container_width=True, hide_index=True)

    errors = checks_df[checks_df["Level"] == "Error"] if not checks_df.empty else pd.DataFrame()
    warnings = checks_df[checks_df["Level"] == "Warning"] if not checks_df.empty else pd.DataFrame()
    if not errors.empty:
        st.error("Errors detected:")
        for _, row in errors.iterrows():
            st.write(f"- {row['Section']} / {row['Check']}: {row['Message']} {row['Value']}")
    if not warnings.empty:
        st.warning("Warnings detected:")
        for _, row in warnings.iterrows():
            st.write(f"- {row['Section']} / {row['Check']}: {row['Message']} {row['Value']}")

    st.markdown("#### Count matrix summary")
    count_summary = report["counts_summary"].copy()
    count_summary.pop("sample_columns", None)
    st.json(count_summary)

    st.markdown("#### Metadata summary")
    metadata_summary = report["metadata_summary"]
    st.write(f"Metadata sample number: `{metadata_summary.get('metadata_sample_number', 0)}`")
    st.write(f"Matched sample number: `{metadata_summary.get('matched_sample_number', 0)}`")
    missing_meta = metadata_summary.get("samples_in_counts_but_missing_in_metadata", [])
    missing_counts = metadata_summary.get("samples_in_metadata_but_missing_in_counts", [])
    st.write(f"Samples in counts but missing in metadata: `{len(missing_meta)}`")
    if missing_meta:
        st.code(", ".join(missing_meta), language="text")
    st.write(f"Samples in metadata but missing in counts: `{len(missing_counts)}`")
    if missing_counts:
        st.code(", ".join(missing_counts), language="text")

    st.markdown("#### Group summary")
    group_counts = report["group_summary"].get("group_sample_counts", {})
    if group_counts:
        st.dataframe(
            pd.DataFrame(
                [{"Group": group, "Sample count": count} for group, count in group_counts.items()]
            ),
            use_container_width=True,
            hide_index=True,
        )
    else:
        st.info("No group sample counts available.")


def main() -> None:
    st.set_page_config(
        page_title="BL Bulk RNA-seq Explorer",
        page_icon="",
        layout="wide",
        initial_sidebar_state="expanded",
    )
    init_session_state()

    st.title("BL Bulk RNA-seq Explorer")
    st.caption("Python/Streamlit MVP migrated from validated browser prototype")

    render_sidebar()

    setup_tab, upload_tab, validation_tab = st.tabs(["Project Setup", "Upload Data", "Input Validation"])
    with setup_tab:
        render_project_setup_tab()
    with upload_tab:
        render_upload_data_tab()
    with validation_tab:
        render_validation_tab()


if __name__ == "__main__":
    main()
