# Bulk RNA-seq Explorer

A Python/Streamlit-based MVP for interactive bulk RNA-seq count matrix exploration, gene identifier handling, quality control, sample grouping, and normalization.

This project is currently in active prototype development. It was migrated from an earlier browser-only HTML prototype into a more reproducible Python/Streamlit architecture.

## Current Status

**Current app version:** `bulk_rnaseq_explorer_v2_0`

**Current development stage:** Streamlit MVP / active internal prototype

The app currently supports:

* Uploading a raw bulk RNA-seq count matrix.
* Detecting sample columns from the uploaded matrix.
* Detecting whether gene identifiers are mouse Ensembl IDs, gene symbols, mixed, or unclear.
* Converting mouse Ensembl IDs to gene symbols using local mapping resources.
* Merging duplicated processed gene symbols by summing raw counts.
* Creating a processed raw count matrix for downstream analysis.
* Computing sample-level quality control metrics.
* Creating reusable QC grouping sets.
* Generating configurable QC bar plots.
* Exporting QC plots as PNG and true SVG through Plotly/Kaleido.
* Running a first-pass normalization workflow through R/Bioconductor.
* Previewing and exporting normalized matrices with search and pagination.

Planned future features include:

* PCA and sample correlation.
* DESeq2-based differential expression analysis.
* Volcano plots.
* Heatmaps.
* ORA pathway analysis.
* GSEA / fgsea-style pathway analysis.
* Exportable result reports.
* Modular project structure for future integration with single-cell analysis tools.

## Background

The original browser-only HTML prototype was developed as a fast, interactive proof-of-concept for bulk RNA-seq analysis. It included count matrix upload/parsing, sample grouping, gene mapping, QC, normalization, exploratory DEG-style analysis, volcano plots, heatmaps, ORA, and GSEA-style pathway analysis.

The current Python/Streamlit version treats the HTML app as a validated product reference, but rebuilds the workflow with more reproducible Python/R-compatible logic.

The long-term goal is to build a broader scientific analysis platform that can eventually integrate:

* Bulk RNA-seq analysis.
* Single-cell RNA-seq analysis.
* Spatial, flow cytometry, and image analysis modules.
* AI-assisted scientific workflow organization.
* Local and potentially cloud/server-based analysis execution.

## Repository Structure

Current structure:

```text
bulk-rnaseq-explorer/
├── archive/
│   ├── bulk_rnaseq_explorer_v1_0.py
│   ├── bulk_rnaseq_explorer_v1_1.py
│   ├── bulk_rnaseq_explorer_v1_2.py
│   ├── bulk_rnaseq_explorer_v1_3.py
│   ├── bulk_rnaseq_explorer_v1_4.py
│   └── bulk_rnaseq_explorer_v1_5.py
├── assets/
│   └── generated cache files, ignored by Git when appropriate
├── database_raw/
│   └── local GMT pathway databases
├── source_mapping/
│   └── mouse Ensembl-to-gene-symbol mapping resources
├── r_scripts/
│   └── normalize_counts.R
├── outputs/
│   └── normalization/
├── bulk_rnaseq_explorer.py
├── index_v5.4.9.html
├── .gitignore
└── README.md
```

The active app entry point is:

```text
bulk_rnaseq_explorer.py
```

The archived `bulk_rnaseq_explorer_v1_0.py` to `bulk_rnaseq_explorer_v1_5.py` files preserve early migration history.

The current project intentionally still uses one main Python file while the core workflow is being stabilized. A more modular structure can be introduced later after the upload, QC, normalization, DEG, and pathway workflows are stable.

## Installation

Recommended environment: Python 3.10+ in a dedicated Conda environment.

Example:

```bash
conda create -n rnaseq python=3.10
conda activate rnaseq
pip install streamlit pandas numpy plotly kaleido
```

Optional for faster local gene-map cache:

```bash
pip install pyarrow
```

Normalization requires R and Bioconductor packages:

```r
install.packages("BiocManager")
BiocManager::install(c("DESeq2", "edgeR"))
install.packages("jsonlite")
```

`Rscript` must be available from the terminal PATH.

## Running the App

From the project folder:

```bash
conda activate rnaseq
streamlit run bulk_rnaseq_explorer.py
```

The Streamlit app will open in a local browser window.

## Input Format

The app expects a raw bulk RNA-seq count matrix with genes as rows and samples as columns.

Recommended format: tab-delimited `.txt` or `.tsv`.

CSV is also accepted.

Example:

```text
EnsemblID / Gene_symbol    Sample_1    Sample_2    Sample_3    Sample_4
ENSMUSG00000000001         120         98          115         130
ENSMUSG00000000028         0           4           1           3
Cxcl1                      50          80          320         400
Actb                       10000       9800        10300       9900
```

The first column should contain gene identifiers. The remaining columns should contain raw count values for samples.

## Main Workflow

### 1. Upload Count Matrix

The upload workflow supports:

* Raw count matrix upload.
* Automatic sample detection.
* Gene ID column selection.
* Gene ID mode detection:

  * Mouse Ensembl ID.
  * Gene symbol.
  * Mixed.
  * Unknown.
* Mouse Ensembl-to-gene-symbol conversion using local mapping resources.
* Duplicate gene-symbol detection.
* Duplicate gene merging by summing raw counts.
* Processed count matrix generation for downstream analysis.

The processed matrix is the main input for QC and normalization.

### 2. Quality Control

The QC workflow provides:

* Sample-level summary table.
* Reusable QC grouping sets.
* Library size plot.
* Detected genes plot.
* Zero-count fraction plot.
* Plotting by individual samples or saved QC groups.
* Group-level aggregation:

  * Mean.
  * Median.
  * Sum, when appropriate.
* Overlay of individual sample dots on group-level bar plots.
* Plot customization:

  * Width.
  * Height.
  * X-axis label angle.
  * Axis label font size.
  * Bar colors.
* Plot export:

  * PNG.
  * True SVG through Plotly/Kaleido.

### 3. Normalization

The v2.0 normalization workflow adds the first R/Bioconductor-backed normalization module.

Normalization uses the processed raw count matrix after gene-symbol conversion and duplicate-gene merging.

The module generates:

* Raw counts.
* CPM.
* log2(CPM + 1).
* DESeq2 normalized counts.
* DESeq2 VST matrix.
* edgeR TMM-normalized CPM.
* edgeR TMM-normalized logCPM.

DESeq2 and edgeR normalization are run through `Rscript` using:

```text
r_scripts/normalize_counts.R
```

The output files are written to session-specific folders under:

```text
outputs/normalization/
```

Expected output files include:

```text
raw_counts.csv
cpm.csv
log2_cpm_plus1.csv
deseq2_size_factors.csv
deseq2_normalized_counts.csv
deseq2_vst.csv
edger_tmm_norm_factors.csv
edger_tmm_cpm.csv
edger_tmm_logcpm.csv
normalization_report.json
```

The Normalization page includes:

* Matrix summary.
* Total counts per sample.
* Non-integer count warnings.
* Raw-count sanity warnings when applicable.
* Rscript execution status.
* DESeq2 size factors table.
* edgeR TMM normalization factors table.
* Output matrix dimension summary.
* Searchable normalized matrix preview.
* Previous / Next pagination.
* Show rows selector.
* CSV export.

Search behavior:

* Gene search starts only after at least 3 typed characters.
* Search is case-insensitive.
* The table viewer shows gene suggestions when matches are found.
* CSV export downloads the currently search-filtered matrix; without search, it downloads the full selected matrix.

## Output Philosophy

The app separates exploratory UI convenience from reproducible backend outputs.

Current principles:

* Raw count parsing and UI interaction are handled in Python/Streamlit.
* DESeq2 and edgeR normalization are delegated to R/Bioconductor.
* Duplicate genes are merged by summing raw counts.
* All normalization matrices preserve the `Gene` column.
* Sample order is preserved from the uploaded matrix.
* R package errors are surfaced clearly in the UI instead of crashing the app.

## Version History

### v1.0-v1.2: Build initial Streamlit migration skeleton

Main changes:

* Created the first Streamlit MVP from the validated HTML prototype.
* Added count matrix upload and basic input parsing.
* Added early workflow navigation and local resource detection.

### v1.3: Add gene symbol conversion and duplicate merging

Main changes:

* Added mouse Ensembl-to-gene-symbol conversion.
* Added duplicate gene-symbol detection.
* Merged duplicated processed genes by summing raw counts.

### v1.4-v1.5: Add Quality Control foundation

Main changes:

* Added sample-level QC summary.
* Added Library Size, Detected Genes, and Zero-count Fraction plots.
* Added early QC grouping and plot export support.

### v1.6-v1.9: Improve QC grouping and plot controls

Main changes:

* Improved reusable QC grouping sets.
* Added group-level plot aggregation with overlaid sample dots.
* Added advanced plot settings, color settings, reset behavior, and direct PNG/SVG export.

### v1.10-v1.15: Fix QC action layout and stabilize UI state

Main changes:

* Fixed QC grouping save/reset behavior.
* Improved Streamlit widget state handling.
* Cleaned up QC button layout, export buttons, and related UI logic.
* Preserved true SVG export through Plotly/Kaleido.

### v2.0: Add R-backed normalization workflow

Main changes:

* Added a standalone Normalization section.
* Added Rscript-based DESeq2 and edgeR normalization.
* Generated raw counts, CPM, logCPM, DESeq2 normalized counts, DESeq2 VST, edgeR TMM CPM, and edgeR TMM logCPM.
* Added searchable, paginated normalized matrix preview with CSV export.
* Added normalization report output with package/version and filtering metadata.

## Roadmap

Near-term priorities:

* Add PCA and sample correlation after normalization.
* Add DESeq2-based DEG analysis through R.
* Add DEG result tables and volcano plots.
* Add heatmap visualization.
* Add ORA and GSEA / fgsea-style pathway analysis.
* Improve project/session output organization.
* Gradually modularize the codebase after the core workflow stabilizes.

Long-term priorities:

* Integrate with single-cell RNA-seq analysis tools.
* Support additional omics and image-based workflows.
* Add project-level data management.
* Support local, server, and cloud execution modes.
* Build toward a broader scientific analysis platform.
