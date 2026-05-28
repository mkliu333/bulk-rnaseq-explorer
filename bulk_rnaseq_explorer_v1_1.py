"""
Bulk RNA-seq Explorer
Version: bulk_rnaseq_explorer_v1_1

Scope for v1.1:
- Streamlit MVP workflow migrated toward the validated browser prototype.
- Upload one raw count matrix.
- Detect sample columns from the count matrix.
- Assign samples into two groups in the app.
- Validate raw counts and group assignments before future QC.
- Detect and attempt to parse a local mouse Ensembl-to-symbol gene map.

Explicitly out of scope for this version:
- DESeq2, Rscript, DEG analysis, QC plots, PCA, heatmap, volcano plot,
  GSEA, ORA, pathway database parsing, cloud storage, login, Duke DCC, SLURM.
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


APP_VERSION = "bulk_rnaseq_explorer_v1_1"

TEMPLATE_TEXT = """GeneID\tSample_1\tSample_2\tSample_3\tSample_4
ENSMUSG00000000001\t120\t98\t115\t130
ENSMUSG00000000028\t0\t4\t1\t3
Cxcl1\t50\t80\t320\t400
Actb\t10000\t9800\t10300\t9900
"""


def init_session_state() -> None:
    """Initialize persistent app state."""
    defaults: dict[str, Any] = {
        "app_version": APP_VERSION,
        "counts_df": None,
        "counts_file_name": None,
        "counts_file_signature": None,
        "gene_id_column": None,
        "sample_columns": [],
        "sample_group_assignments": {},
        "group_1_name": "Control",
        "group_2_name": "Treatment",
        "gene_map_df": None,
        "gene_map_status": {
            "source_path": None,
            "detected": False,
            "parsed": False,
            "message": "Not loaded",
            "n_mappings": 0,
        },
        "validation_status": "Not run",
        "validation_report": None,
        "downstream_invalidated": False,
        "qc_results": None,
        "deg_results": None,
        "pathway_results": None,
        "plots": None,
    }
    for key, value in defaults.items():
        if key not in st.session_state:
            st.session_state[key] = value


def invalidate_downstream_results() -> None:
    """Clear future-result placeholders after input or key setting changes."""
    st.session_state["downstream_invalidated"] = True
    st.session_state["validation_status"] = "Not run"
    st.session_state["validation_report"] = None
    st.session_state["qc_results"] = None
    st.session_state["deg_results"] = None
    st.session_state["pathway_results"] = None
    st.session_state["plots"] = None


def detect_local_resources() -> pd.DataFrame:
    """Detect local files and folders used as historical/reference resources."""
    cwd = Path.cwd()
    resource_specs = [
        ("Historical browser prototype", "index_v5.4.9.html", "File"),
        ("Top-level mouse gene map JS", "mouse_ensembl_to_symbol.js", "File"),
        ("Gene sets JS bundle", "gene_sets_2026_v1.js", "File"),
        ("GMT database directory", "database_raw", "Directory"),
        ("Source mapping directory", "source_mapping", "Directory"),
        ("Source mapping directory with space", "source mapping", "Directory"),
        ("Source mapping table", "source_mapping/mouse_ensembl_to_symbol", "File"),
        ("Source mapping table with space", "source mapping/mouse_ensembl_to_symbol", "File"),
        ("Source GTF archive", "source_mapping/Mus_musculus.GRCm39.115.gtf.gz", "File"),
        ("Source GTF archive with space", "source mapping/Mus_musculus.GRCm39.115.gtf.gz", "File"),
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


def _mapping_dataframe_from_object(data: Any) -> pd.DataFrame | None:
    """Convert common gene map JSON shapes into a two-column dataframe."""
    if isinstance(data, dict):
        rows = []
        for key, value in data.items():
            if isinstance(value, str):
                rows.append({"ensembl_id": key, "gene_symbol": value})
            elif isinstance(value, dict):
                symbol = (
                    value.get("gene_symbol")
                    or value.get("symbol")
                    or value.get("gene_name")
                    or value.get("name")
                )
                if symbol is not None:
                    rows.append({"ensembl_id": key, "gene_symbol": symbol})
        return pd.DataFrame(rows) if rows else None

    if isinstance(data, list):
        rows = []
        ensembl_keys = ["ensembl_id", "ensembl", "ensembl_gene_id", "gene_id", "id"]
        symbol_keys = ["gene_symbol", "symbol", "gene_name", "name", "external_gene_name"]
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
    """Extract a JSON-like object from a JavaScript assignment or declaration."""
    assignment_match = re.search(r"=\s*(\{.*\})\s*;?\s*$", text, flags=re.DOTALL)
    object_text = assignment_match.group(1) if assignment_match else text.strip()
    object_text = re.sub(r"^\s*(?:const|let|var)\s+\w+\s*=\s*", "", object_text, flags=re.DOTALL)
    object_text = object_text.rstrip(";").strip()
    return json.loads(object_text)


def parse_gene_map_file(path: Path) -> tuple[pd.DataFrame | None, dict[str, Any]]:
    """Try to parse a local gene map file without making it required for analysis."""
    status = {
        "source_path": str(path),
        "detected": path.exists(),
        "parsed": False,
        "message": "",
        "n_mappings": 0,
    }

    if not path.exists():
        status["message"] = "File not detected"
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

        gene_map_df = _mapping_dataframe_from_object(data)
        if gene_map_df is None or gene_map_df.empty:
            status["message"] = "Detected but not parsed"
            return None, status

        gene_map_df = gene_map_df.dropna().drop_duplicates()
        status["parsed"] = True
        status["message"] = "Parsed successfully"
        status["n_mappings"] = int(gene_map_df.shape[0])
        return gene_map_df, status
    except Exception as exc:
        status["message"] = f"Detected but not parsed: {exc}"
        return None, status


def load_local_gene_map() -> tuple[pd.DataFrame | None, dict[str, Any]]:
    """Load the first available local gene map using the requested priority."""
    cwd = Path.cwd()
    candidates = [
        cwd / "source_mapping" / "mouse_ensembl_to_symbol",
        cwd / "source mapping" / "mouse_ensembl_to_symbol",
        cwd / "source_mapping" / "mouse_ensembl_to_symbol.json",
        cwd / "source mapping" / "mouse_ensembl_to_symbol.json",
        cwd / "mouse_ensembl_to_symbol.js",
    ]

    first_detected_status: dict[str, Any] | None = None
    for path in candidates:
        if not path.exists():
            continue
        gene_map_df, status = parse_gene_map_file(path)
        if first_detected_status is None:
            first_detected_status = status
        if gene_map_df is not None:
            return gene_map_df, status

    if first_detected_status is not None:
        return None, first_detected_status

    return None, {
        "source_path": None,
        "detected": False,
        "parsed": False,
        "message": "No local gene map detected",
        "n_mappings": 0,
    }


def read_count_matrix_file(uploaded_file) -> pd.DataFrame:
    """Read an uploaded count matrix, preferring tab-delimited files."""
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
    if counts_df is None or gene_id_column not in counts_df.columns:
        return []
    return [str(column) for column in counts_df.columns if str(column) != str(gene_id_column)]


def _make_check(section: str, check: str, level: str, message: str, value: Any = "") -> dict[str, Any]:
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
    """Validate raw count matrix structure and values."""
    checks: list[dict[str, Any]] = []
    summary: dict[str, Any] = {
        "number_of_genes": 0,
        "number_of_samples": 0,
        "number_of_duplicated_gene_ids": 0,
        "number_of_duplicated_sample_names": 0,
        "number_of_all_zero_genes": 0,
        "number_of_low_count_genes": 0,
        "number_of_constant_genes": 0,
        "number_of_genes_with_missing_values": 0,
        "number_of_non_integer_count_entries": 0,
        "sample_columns": [],
    }

    if counts_df is None:
        checks.append(_make_check("Counts", "Count matrix uploaded", "Error", "Count matrix is not uploaded."))
        return {"checks": checks, "summary": summary}

    checks.append(_make_check("Counts", "Count matrix uploaded", "Pass", "Count matrix was uploaded."))

    column_names = [str(column) for column in counts_df.columns]
    if gene_id_column not in column_names:
        checks.append(
            _make_check("Counts", "Gene ID column exists", "Error", "Selected gene ID column is not present.")
        )
        return {"checks": checks, "summary": summary}

    checks.append(_make_check("Counts", "Gene ID column exists", "Pass", "Selected gene ID column is present."))

    sample_columns = [column for column in column_names if column != str(gene_id_column)]
    summary["number_of_genes"] = int(counts_df.shape[0])
    summary["number_of_samples"] = int(len(sample_columns))
    summary["sample_columns"] = sample_columns

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

    gene_ids = counts_df[str(gene_id_column)]
    gene_id_text = gene_ids.astype("string").str.strip()
    missing_gene_ids = int(gene_ids.isna().sum() + gene_id_text.eq("").sum())
    duplicated_gene_ids = int(gene_id_text.duplicated(keep=False).sum())
    summary["number_of_duplicated_gene_ids"] = duplicated_gene_ids

    if missing_gene_ids > 0:
        checks.append(_make_check("Counts", "Missing gene IDs", "Error", "Gene ID column contains missing values.", missing_gene_ids))
    else:
        checks.append(_make_check("Counts", "Missing gene IDs", "Pass", "No missing gene IDs detected.", 0))

    if duplicated_gene_ids > 0:
        checks.append(_make_check("Counts", "Duplicated gene IDs", "Warning", "Duplicated gene IDs detected.", duplicated_gene_ids))
    else:
        checks.append(_make_check("Counts", "Duplicated gene IDs", "Pass", "No duplicated gene IDs detected.", 0))

    duplicated_sample_names = int(pd.Index(sample_columns).duplicated(keep=False).sum())
    summary["number_of_duplicated_sample_names"] = duplicated_sample_names
    if duplicated_sample_names > 0:
        checks.append(
            _make_check("Counts", "Duplicated sample names", "Error", "Duplicated sample names detected.", duplicated_sample_names)
        )
    else:
        checks.append(_make_check("Counts", "Duplicated sample names", "Pass", "No duplicated sample names detected.", 0))

    raw_counts = counts_df[sample_columns] if sample_columns else pd.DataFrame()
    count_values = raw_counts.apply(pd.to_numeric, errors="coerce") if sample_columns else pd.DataFrame()
    non_numeric_columns = []
    for column in sample_columns:
        original_text = raw_counts[column].astype("string").str.strip()
        non_missing_original = raw_counts[column].notna() & original_text.ne("")
        coerced_missing = count_values[column].isna()
        if (non_missing_original & coerced_missing).any():
            non_numeric_columns.append(column)

    if non_numeric_columns:
        checks.append(
            _make_check(
                "Counts",
                "Sample columns numeric",
                "Error",
                "One or more sample columns contain non-numeric values.",
                ", ".join(non_numeric_columns[:10]),
            )
        )
    else:
        checks.append(_make_check("Counts", "Sample columns numeric", "Pass", "All sample columns are numeric."))

    missing_entries = int(count_values.isna().sum().sum()) if not count_values.empty else 0
    genes_with_missing = int(count_values.isna().any(axis=1).sum()) if not count_values.empty else 0
    summary["number_of_genes_with_missing_values"] = genes_with_missing
    if missing_entries > 0:
        checks.append(_make_check("Counts", "Missing count values", "Warning", "Missing values detected in count matrix.", missing_entries))
    else:
        checks.append(_make_check("Counts", "Missing count values", "Pass", "No missing count values detected.", 0))

    negative_entries = int((count_values < 0).sum().sum()) if not count_values.empty else 0
    if negative_entries > 0:
        checks.append(_make_check("Counts", "Negative counts", "Error", "Negative count values detected.", negative_entries))
    else:
        checks.append(_make_check("Counts", "Negative counts", "Pass", "No negative count values detected.", 0))

    if not count_values.empty:
        valid_mask = count_values.notna()
        non_integer_mask = valid_mask & ((count_values % 1) != 0)
        non_integer_entries = int(non_integer_mask.sum().sum())
    else:
        non_integer_entries = 0
    summary["number_of_non_integer_count_entries"] = non_integer_entries
    if non_integer_entries > 0:
        checks.append(_make_check("Counts", "Non-integer counts", "Warning", "DESeq2 expects raw integer counts.", non_integer_entries))
    else:
        checks.append(_make_check("Counts", "Non-integer counts", "Pass", "All numeric count entries are integers.", 0))

    if not count_values.empty:
        row_sums = count_values.fillna(0).sum(axis=1)
        all_zero_genes = int(row_sums.eq(0).sum())
        low_count_genes = int(row_sums.lt(10).sum())
        constant_genes = int(count_values.nunique(axis=1, dropna=False).eq(1).sum())
    else:
        all_zero_genes = 0
        low_count_genes = 0
        constant_genes = 0

    summary["number_of_all_zero_genes"] = all_zero_genes
    summary["number_of_low_count_genes"] = low_count_genes
    summary["number_of_constant_genes"] = constant_genes

    if all_zero_genes > 0:
        checks.append(_make_check("Counts", "All-zero genes", "Warning", "All-zero genes detected.", all_zero_genes))
    else:
        checks.append(_make_check("Counts", "All-zero genes", "Pass", "No all-zero genes detected.", 0))

    if low_count_genes > 0:
        checks.append(_make_check("Counts", "Low-count genes", "Warning", "Genes with total count < 10 detected.", low_count_genes))
    else:
        checks.append(_make_check("Counts", "Low-count genes", "Pass", "No low-count genes detected.", 0))

    if constant_genes > 0:
        checks.append(_make_check("Counts", "Constant genes", "Warning", "Genes with identical counts across all samples detected.", constant_genes))
    else:
        checks.append(_make_check("Counts", "Constant genes", "Pass", "No constant genes detected.", 0))

    return {"checks": checks, "summary": summary}


def validate_group_assignments(
    sample_columns: list[str],
    sample_group_assignments: dict[str, str],
    group_1_name: str,
    group_2_name: str,
) -> dict[str, Any]:
    """Validate two-group sample assignment."""
    checks: list[dict[str, Any]] = []

    if group_1_name.strip() == "" or group_2_name.strip() == "":
        checks.append(_make_check("Groups", "Group names provided", "Error", "Both group names must be non-empty."))
    else:
        checks.append(_make_check("Groups", "Group names provided", "Pass", "Both group names are non-empty."))

    if group_1_name == group_2_name:
        checks.append(_make_check("Groups", "Group names distinct", "Error", "Group 1 and Group 2 names must be different."))
    else:
        checks.append(_make_check("Groups", "Group names distinct", "Pass", "Group names are distinct."))

    valid_groups = [group_1_name, group_2_name]
    assigned_groups = {
        sample: group
        for sample, group in sample_group_assignments.items()
        if group in valid_groups and sample in sample_columns
    }
    group_counts = {group_1_name: 0, group_2_name: 0}
    for group in assigned_groups.values():
        group_counts[group] += 1

    unassigned_samples = [
        sample
        for sample in sample_columns
        if sample_group_assignments.get(sample, "Unassigned") not in valid_groups
    ]
    assignment_missing_from_counts = sorted(set(sample_group_assignments) - set(sample_columns))
    counts_missing_from_assignment = sorted(set(sample_columns) - set(sample_group_assignments))
    selected_groups = sorted(set(assigned_groups.values()))

    summary = {
        "group_counts": group_counts,
        "unassigned_samples": unassigned_samples,
        "samples_in_assignment_but_missing_from_counts": assignment_missing_from_counts,
        "samples_in_counts_but_missing_from_assignment": counts_missing_from_assignment,
        "selected_groups": selected_groups,
    }

    if unassigned_samples:
        checks.append(
            _make_check(
                "Groups",
                "All samples assigned",
                "Warning",
                "Some count matrix samples are unassigned.",
                len(unassigned_samples),
            )
        )
    else:
        checks.append(_make_check("Groups", "All samples assigned", "Pass", "All samples are assigned.", 0))

    if len(selected_groups) != 2:
        checks.append(
            _make_check(
                "Groups",
                "Exactly two groups selected",
                "Error",
                "Both configured groups must contain at least one assigned sample.",
                len(selected_groups),
            )
        )
    else:
        checks.append(_make_check("Groups", "Exactly two groups selected", "Pass", "Exactly two groups selected.", 2))

    empty_groups = [group for group, count in group_counts.items() if count < 1]
    if empty_groups:
        checks.append(_make_check("Groups", "Each group has at least 1 sample", "Error", "At least one group has no samples.", ", ".join(empty_groups)))
    else:
        checks.append(_make_check("Groups", "Each group has at least 1 sample", "Pass", "Both groups have at least 1 sample."))

    single_replicate_groups = [group for group, count in group_counts.items() if count == 1]
    if single_replicate_groups:
        checks.append(
            _make_check(
                "Groups",
                "Each group has at least 2 samples",
                "Warning",
                "DESeq2 requires biological replicates for reliable inference.",
                ", ".join(single_replicate_groups),
            )
        )
    elif not empty_groups:
        checks.append(_make_check("Groups", "Each group has at least 2 samples", "Pass", "Both groups have at least 2 samples."))

    if assignment_missing_from_counts:
        checks.append(
            _make_check(
                "Groups",
                "No assignment sample missing from count matrix",
                "Error",
                "Group assignment contains samples not present in the count matrix.",
                ", ".join(assignment_missing_from_counts[:10]),
            )
        )
    else:
        checks.append(_make_check("Groups", "No assignment sample missing from count matrix", "Pass", "No extra assignment samples detected.", 0))

    if counts_missing_from_assignment:
        checks.append(
            _make_check(
                "Groups",
                "No count matrix sample missing from assignment",
                "Warning",
                "Some count matrix samples do not have assignment records.",
                len(counts_missing_from_assignment),
            )
        )
    else:
        checks.append(_make_check("Groups", "No count matrix sample missing from assignment", "Pass", "All count matrix samples have assignment records.", 0))

    return {"checks": checks, "summary": summary}


def run_full_validation() -> dict[str, Any]:
    """Run count and group validation, then store report in session state."""
    counts_df = st.session_state["counts_df"]
    if counts_df is None:
        report = {
            "checks": [],
            "counts_summary": {},
            "group_summary": {},
            "status": "Not run",
        }
        st.session_state["validation_status"] = "Not run"
        st.session_state["validation_report"] = report
        return report

    counts_report = validate_counts(counts_df, st.session_state["gene_id_column"])
    sample_columns = counts_report["summary"].get("sample_columns", [])
    group_report = validate_group_assignments(
        sample_columns,
        st.session_state["sample_group_assignments"],
        st.session_state["group_1_name"],
        st.session_state["group_2_name"],
    )

    all_checks = counts_report["checks"] + group_report["checks"]
    if _count_severity(all_checks, "Error") > 0:
        status = "Failed"
    elif _count_severity(all_checks, "Warning") > 0:
        status = "Warning"
    else:
        status = "Passed"

    report = {
        "checks": all_checks,
        "counts_summary": counts_report["summary"],
        "group_summary": group_report["summary"],
        "status": status,
    }
    st.session_state["validation_status"] = status
    st.session_state["validation_report"] = report
    return report


def _display_gene_map_status() -> None:
    status = st.session_state["gene_map_status"]
    st.write(f"Gene map source path: `{status.get('source_path') or 'None'}`")
    st.write(f"Gene map detected: {'Yes' if status.get('detected') else 'No'}")
    st.write(f"Gene map parsed: {'Yes' if status.get('parsed') else 'No'}")
    st.write(f"Number of gene mappings: `{status.get('n_mappings', 0)}`")
    message = status.get("message", "")
    if status.get("parsed"):
        st.success(message)
    elif status.get("detected"):
        st.warning(message or "Detected but not parsed")
    else:
        st.info(message or "No local gene map detected")


def render_sidebar() -> None:
    """Render workflow and project status sidebar."""
    st.sidebar.title("Workflow")
    st.sidebar.success("1. Upload Count Matrix")
    st.sidebar.success("2. Assign Sample Groups")
    st.sidebar.success("3. Input Validation")
    st.sidebar.caption("4. QC Overview - Coming soon / Locked")
    st.sidebar.caption("5. DEG Analysis - Coming soon / Locked")
    st.sidebar.caption("6. Visualization - Coming soon / Locked")
    st.sidebar.caption("7. Pathway Analysis - Coming soon / Locked")
    st.sidebar.caption("8. Export - Coming soon / Locked")

    sample_columns = st.session_state.get("sample_columns", [])
    assignments = st.session_state.get("sample_group_assignments", {})
    group_1 = st.session_state.get("group_1_name", "Control")
    group_2 = st.session_state.get("group_2_name", "Treatment")
    groups_assigned = (
        sum(1 for sample in sample_columns if assignments.get(sample) == group_1) >= 1
        and sum(1 for sample in sample_columns if assignments.get(sample) == group_2) >= 1
    )

    st.sidebar.divider()
    st.sidebar.subheader("Project status")
    st.sidebar.write(f"Count matrix uploaded: {'Yes' if st.session_state['counts_df'] is not None else 'No'}")
    st.sidebar.write(f"Number of samples: {len(sample_columns)}")
    st.sidebar.write(f"Groups assigned: {'Yes' if groups_assigned else 'No'}")
    st.sidebar.write(f"Gene map detected: {'Yes' if st.session_state['gene_map_status'].get('detected') else 'No'}")
    st.sidebar.write(f"Validation status: {st.session_state['validation_status']}")

    if st.session_state["downstream_invalidated"]:
        st.sidebar.warning("Downstream placeholders invalidated.")

    with st.sidebar.expander("Local resources", expanded=False):
        resources = detect_local_resources().copy()
        resources["Detected"] = resources["Detected"].map(lambda value: "Yes" if value else "No")
        st.dataframe(resources, use_container_width=True, hide_index=True)
        _display_gene_map_status()


def _reset_counts_dependent_state(counts_df: pd.DataFrame, file_name: str) -> None:
    st.session_state["counts_df"] = counts_df
    st.session_state["counts_file_name"] = file_name
    st.session_state["gene_id_column"] = str(counts_df.columns[0]) if len(counts_df.columns) else None
    st.session_state["sample_columns"] = get_sample_columns(counts_df, st.session_state["gene_id_column"])
    st.session_state["sample_group_assignments"] = {}
    invalidate_downstream_results()


def _uploaded_file_signature(uploaded_file) -> str:
    """Create a stable signature so same-name replacement uploads are detected."""
    uploaded_file.seek(0)
    content = uploaded_file.getvalue()
    uploaded_file.seek(0)
    return hashlib.sha256(content).hexdigest()


def render_upload_count_matrix_tab() -> None:
    """Render count matrix upload, template, examples, and preview."""
    st.subheader("Upload Count Matrix")
    st.write(
        "Please upload a raw count matrix as a tab-delimited file (.tsv or .txt). "
        "The first column should contain Ensembl IDs or gene symbols, and the first row should contain sample names."
    )

    st.download_button(
        "Download sample template",
        data=TEMPLATE_TEXT,
        file_name="bulk_rnaseq_count_matrix_template.tsv",
        mime="text/tab-separated-values",
    )

    st.markdown("#### Sample template")
    st.dataframe(pd.read_csv(StringIO(TEMPLATE_TEXT), sep="\t"), use_container_width=True, hide_index=True)

    st.markdown("#### Example input layout")
    example_layout = pd.DataFrame(
        {
            "EnsemblID / Gene symbols": ["ENSMUSG00000000001", "ENSMUSG00000000028", "Cxcl1"],
            "Sample 1": [120, 0, 50],
            "Sample 2": [98, 4, 80],
            "Sample 3": [115, 1, 320],
        }
    )
    st.dataframe(example_layout, use_container_width=True, hide_index=True)

    uploaded_file = st.file_uploader(
        "Upload raw count matrix",
        type=["tsv", "txt", "csv"],
        key="count_matrix_uploader",
        help="Tab-delimited .tsv or .txt is recommended. CSV is also accepted.",
    )

    if uploaded_file is not None:
        uploaded_signature = _uploaded_file_signature(uploaded_file)
    else:
        uploaded_signature = None

    if uploaded_file is not None and st.session_state.get("counts_file_signature") != uploaded_signature:
        try:
            counts_df = read_count_matrix_file(uploaded_file)
            _reset_counts_dependent_state(counts_df, uploaded_file.name)
            st.session_state["counts_file_signature"] = uploaded_signature
            st.success("Count matrix uploaded.")
        except Exception as exc:
            st.error(f"Failed to read `{uploaded_file.name}`: {exc}")

    counts_df = st.session_state["counts_df"]
    if counts_df is not None:
        columns = [str(column) for column in counts_df.columns]
        st.success(f"Uploaded file: `{st.session_state['counts_file_name']}`")

        current_gene_id = st.session_state["gene_id_column"] if st.session_state["gene_id_column"] in columns else columns[0]
        selected_gene_id = st.selectbox(
            "Detected gene ID column",
            options=columns,
            index=columns.index(current_gene_id),
            key="gene_id_column_selector_v1_1",
        )
        if selected_gene_id != st.session_state["gene_id_column"]:
            st.session_state["gene_id_column"] = selected_gene_id
            st.session_state["sample_columns"] = get_sample_columns(counts_df, selected_gene_id)
            st.session_state["sample_group_assignments"] = {}
            invalidate_downstream_results()

        sample_columns = st.session_state["sample_columns"]
        st.write(f"Detected shape: `{counts_df.shape[0]}` genes x `{len(sample_columns)}` samples")
        st.write(f"Detected gene ID column: `{st.session_state['gene_id_column']}`")
        st.markdown("#### Preview")
        st.dataframe(counts_df.head(5), use_container_width=True)

        st.markdown("#### Sample columns")
        if sample_columns:
            st.dataframe(pd.DataFrame({"Sample": sample_columns}), use_container_width=True, hide_index=True)
        else:
            st.warning("No sample columns detected after selecting the gene ID column.")

    with st.expander("Local resource and gene map detection", expanded=False):
        resources = detect_local_resources().copy()
        resources["Detected"] = resources["Detected"].map(lambda value: "Yes" if value else "No")
        st.dataframe(resources, use_container_width=True, hide_index=True)
        _display_gene_map_status()


def _rename_group_assignments(old_name: str, new_name: str) -> None:
    if old_name == new_name:
        return
    assignments = st.session_state["sample_group_assignments"]
    st.session_state["sample_group_assignments"] = {
        sample: (new_name if group == old_name else group)
        for sample, group in assignments.items()
    }


def render_assign_sample_groups_tab() -> None:
    """Render in-app two-group sample assignment controls."""
    st.subheader("Assign Sample Groups")

    counts_df = st.session_state["counts_df"]
    if counts_df is None:
        st.info("Please upload a count matrix first.")
        return

    sample_columns = st.session_state["sample_columns"]
    if not sample_columns:
        st.warning("No sample columns are available. Check the selected gene ID column.")
        return

    col1, col2 = st.columns(2)
    with col1:
        group_1_input = st.text_input("Group 1 name", value=st.session_state["group_1_name"])
    with col2:
        group_2_input = st.text_input("Group 2 name", value=st.session_state["group_2_name"])

    if group_1_input != st.session_state["group_1_name"]:
        _rename_group_assignments(st.session_state["group_1_name"], group_1_input)
        st.session_state["group_1_name"] = group_1_input
        invalidate_downstream_results()
    if group_2_input != st.session_state["group_2_name"]:
        _rename_group_assignments(st.session_state["group_2_name"], group_2_input)
        st.session_state["group_2_name"] = group_2_input
        invalidate_downstream_results()

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
            invalidate_downstream_results()
            st.rerun()
    with action_col2:
        if st.button("Clear all assignments"):
            st.session_state["sample_group_assignments"] = {sample: "Unassigned" for sample in sample_columns}
            for sample in sample_columns:
                st.session_state[f"group_assignment_{sample}"] = "Unassigned"
            invalidate_downstream_results()
            st.rerun()

    st.markdown("#### Sample assignments")
    options = ["Unassigned", group_1, group_2]
    for sample in sample_columns:
        current_value = st.session_state["sample_group_assignments"].get(sample, "Unassigned")
        if current_value not in options:
            current_value = "Unassigned"
        widget_key = f"group_assignment_{sample}"
        if st.session_state.get(widget_key) not in options:
            st.session_state[widget_key] = current_value
        selected_group = st.selectbox(
            sample,
            options=options,
            index=options.index(current_value),
            key=widget_key,
        )
        if selected_group != st.session_state["sample_group_assignments"].get(sample, "Unassigned"):
            st.session_state["sample_group_assignments"][sample] = selected_group
            invalidate_downstream_results()

    assignments = st.session_state["sample_group_assignments"]
    summary_rows = [
        {"Sample": sample, "Assigned group": assignments.get(sample, "Unassigned")}
        for sample in sample_columns
    ]
    summary_df = pd.DataFrame(summary_rows)

    group_1_count = int((summary_df["Assigned group"] == group_1).sum())
    group_2_count = int((summary_df["Assigned group"] == group_2).sum())
    unassigned_count = int((summary_df["Assigned group"] == "Unassigned").sum())

    st.markdown("#### Group assignment summary")
    st.dataframe(summary_df, use_container_width=True, hide_index=True)

    metric_col1, metric_col2, metric_col3 = st.columns(3)
    metric_col1.metric(f"{group_1}", group_1_count)
    metric_col2.metric(f"{group_2}", group_2_count)
    metric_col3.metric("Unassigned", unassigned_count)

    if group_1_count >= 1 and group_2_count >= 1:
        st.success("Groups assigned: Yes")
    else:
        st.warning("Groups assigned: No. Assign at least one sample to each group.")


def _render_validation_status(status: str) -> None:
    if status == "Passed":
        st.success("Validation passed.")
    elif status == "Warning":
        st.warning("Validation completed with warnings.")
    elif status == "Failed":
        st.error("Validation failed. Resolve errors before continuing to QC.")
    else:
        st.info("Upload a count matrix and assign sample groups to run validation.")


def render_validation_tab() -> None:
    """Render validation controls and reports."""
    st.subheader("Input Validation")

    if st.session_state["counts_df"] is None:
        st.info("Please upload a count matrix first.")
        _render_validation_status(st.session_state["validation_status"])
        return

    if st.button("Run validation", type="primary"):
        run_full_validation()

    _render_validation_status(st.session_state["validation_status"])
    report = st.session_state["validation_report"]
    if report is None:
        st.info("Validation has not been run for the current inputs.")
        return

    checks_df = pd.DataFrame(report["checks"])
    st.markdown("#### Validation checks")
    st.dataframe(checks_df, use_container_width=True, hide_index=True)

    if not checks_df.empty:
        errors = checks_df[checks_df["Level"] == "Error"]
        warnings = checks_df[checks_df["Level"] == "Warning"]
        if not errors.empty:
            st.error("Errors detected:")
            for _, row in errors.iterrows():
                st.write(f"- {row['Section']} / {row['Check']}: {row['Message']} {row['Value']}")
        if not warnings.empty:
            st.warning("Warnings detected:")
            for _, row in warnings.iterrows():
                st.write(f"- {row['Section']} / {row['Check']}: {row['Message']} {row['Value']}")

    st.markdown("#### Count matrix summary")
    count_summary = dict(report.get("counts_summary", {}))
    count_summary.pop("sample_columns", None)
    st.json(count_summary)

    st.markdown("#### Group assignment summary")
    group_summary = report.get("group_summary", {})
    group_counts = group_summary.get("group_counts", {})
    if group_counts:
        st.dataframe(
            pd.DataFrame(
                [{"Group": group, "Sample count": count} for group, count in group_counts.items()]
            ),
            use_container_width=True,
            hide_index=True,
        )

    unassigned = group_summary.get("unassigned_samples", [])
    if unassigned:
        st.write(f"Unassigned samples: `{len(unassigned)}`")
        st.code(", ".join(unassigned), language="text")

    extra_assignments = group_summary.get("samples_in_assignment_but_missing_from_counts", [])
    if extra_assignments:
        st.write(f"Samples in assignment but missing from count matrix: `{len(extra_assignments)}`")
        st.code(", ".join(extra_assignments), language="text")

    missing_assignments = group_summary.get("samples_in_counts_but_missing_from_assignment", [])
    if missing_assignments:
        st.write(f"Count matrix samples missing from assignment records: `{len(missing_assignments)}`")
        st.code(", ".join(missing_assignments), language="text")


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

    upload_tab, groups_tab, validation_tab = st.tabs(
        ["Upload Count Matrix", "Assign Sample Groups", "Input Validation"]
    )
    with upload_tab:
        render_upload_count_matrix_tab()
    with groups_tab:
        render_assign_sample_groups_tab()
    with validation_tab:
        render_validation_tab()


if __name__ == "__main__":
    main()
