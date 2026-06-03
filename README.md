# Bulk RNA-seq Explorer

A Python/Streamlit-based software for interactive bulk RNA-seq count matrix exploration, quality control, sample grouping, normalization and DEG & Enriched pathway analysis.

This project was migrated from an earlier browser-only HTML prototype into a more reproducible Python/Streamlit + R/Bioconductor architecture.

## Current Status

Current app version: `bulk_rnaseq_explorer_v2_5`

The app currently supports:

* Automatic mouse Ensembl ID to gene symbol conversion using local mapping resources.
* Duplicated gene-symbol detection and merging by summing raw counts.
* Sample-level quality control with reusable QC grouping sets.
* Configurable QC visualizations, including library size, detected genes, zero-count fraction, PCA, and sample correlation.
* PNG and true SVG plot export through Plotly/Kaleido.
* R/Bioconductor-backed normalization using CPM, log2(CPM + 1), DESeq2 normalized counts, DESeq2 VST, edgeR TMM CPM, and edgeR TMM logCPM.
* Searchable and paginated normalized matrix preview with manual CSV export.

Planned future features include:

* DESeq2-based differential expression analysis.
* DEG result tables and volcano plots.
* DEG heatmaps and normalized-expression visualization.
* ORA pathway enrichment analysis.
* GSEA / fgsea-style pathway analysis.
* Exportable analysis reports.
* Future modular integration with single-cell analysis tools.

## Background

The current Python/Streamlit version treats the previous HTML-based UI version as a validated product reference, but rebuilds the workflow with more reproducible Python/R-compatible logic.

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
│   └── local GMT pathway database resources
├── source_mapping/
│   └── mouse Ensembl-to-gene-symbol mapping resources
├── r_scripts/
│   └── normalize_counts.R
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

The current project intentionally still uses one main Python file while the core workflow is being stabilized. A more modular structure can be introduced later after upload, QC, normalization, DEG, and pathway workflows are stable.

## Installation

Recommended environment: Python 3.10+ in a dedicated Conda environment.

Create and activate the environment:

```bash
conda create -n rnaseq python=3.10
conda activate rnaseq
```

Install Python dependencies:

```bash
pip install streamlit pandas numpy plotly kaleido
```

Optional, for faster local gene-map cache:

```bash
pip install pyarrow
```

Normalization requires local R installation and Rscript availability from terminal PATH.

Install R packages from R or RStudio:

```r
install.packages("BiocManager")
install.packages("jsonlite")
BiocManager::install(c("DESeq2", "edgeR"))
```

Check that Rscript is available:

```bash
Rscript --version
```

Check that required R packages can load:

```bash
Rscript -e "library(DESeq2); library(edgeR); library(jsonlite); sessionInfo()"
```

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

The upload workflow parses the raw count matrix, detects sample columns, identifies the gene ID mode, converts mouse Ensembl IDs to gene symbols when possible, and merges duplicated processed gene symbols by summing raw counts.

The resulting processed raw count matrix is used for QC and normalization.

### 2. Quality Control

The QC workflow provides sample-level QC summaries, reusable QC grouping sets, and configurable QC plots.

Current QC visualizations include:

* Library size.
* Detected genes.
* Zero-count fraction.
* PCA using normalized expression matrices.
* Sample correlation heatmap using normalized expression matrices.

QC plots support sample-level or saved-group coloring where appropriate, advanced visual settings, color customization, PNG export, and true SVG export.

### 3. Normalization

The normalization workflow uses the processed raw count matrix after gene-symbol conversion and duplicate-gene merging.

Normalization is performed through R/Bioconductor using:

```text
r_scripts/normalize_counts.R
```

The module generates:

* Raw counts.
* CPM.
* log2(CPM + 1).
* DESeq2 normalized counts.
* DESeq2 VST.
* edgeR TMM CPM.
* edgeR TMM logCPM.

As of v2.5, normalization results are loaded back into Streamlit memory instead of being automatically exported into a permanent project-level output folder. Users can inspect matrices in the UI and export selected tables manually as CSV.

### 4. Normalized Matrix Viewer

The normalized matrix viewer supports:

* Matrix selection.
* Gene search after at least 3 typed characters.
* Matched-gene selection.
* Previous / Next pagination.
* Show rows control.
* CSV export.

## Version History

### v1.0-v1.2: Initial Streamlit migration

Main changes:

* Created the first Streamlit MVP from the validated HTML prototype.
* Added count matrix upload, basic parsing, workflow navigation, and local resource detection.

### v1.3: Gene symbol conversion and duplicate merging

Main changes:

* Added mouse Ensembl-to-gene-symbol conversion.
* Added duplicate gene-symbol detection and merging by summing raw counts.
* Generated a processed raw count matrix for downstream workflows.

### v1.4-v1.5: Initial Quality Control workflow

Main changes:

* Added sample-level QC summary.
* Added Library Size, Detected Genes, and Zero-count Fraction plots.
* Added early QC grouping and plot export support.

### v1.6-v1.9: QC plot interaction and export refinement

Main changes:

* Migrated more QC barplot behavior from the previous HTML version.
* Added sample/group plotting modes, group aggregation, sample-dot overlays, advanced plot settings, and color settings.
* Improved reset behavior, export caching, and PNG/SVG export.

### v1.10-v1.15: QC layout, grouping stability, and UI polish

Main changes:

* Improved QC grouping save/reset behavior.
* Stabilized Streamlit widget state handling using reset nonces and cleanup logic.
* Reworked responsive button layouts and reduced stale state issues.
* Improved barplot responsiveness and visual consistency.

### v2.0: R/Bioconductor-backed normalization

Main changes:

* Added normalization as an independent workflow section.
* Added Rscript-based DESeq2 and edgeR normalization.
* Added normalized matrix viewer with search, pagination, and CSV export.

### v2.1: Normalization UI simplification and QC performance cleanup

Main changes:

* Removed redundant normalization summary/output displays.
* Changed normalization to auto-run when valid processed counts are available.
* Improved matrix viewer controls and reduced repeated QC loading overhead.

### v2.2: PCA and sample correlation

Main changes:

* Added PCA QC plot using normalized expression matrices.
* Added sample correlation heatmap.
* Added QC grouping-based coloring and annotation support for PCA and correlation plots.

### v2.3-v2.4: QC view layout and grouping-order fixes

Main changes:

* Reorganized Quality Control into compact selectable views.
* Removed approximate VST behavior inherited from the HTML prototype.
* Improved saved QC grouping order preservation in PCA and sample correlation.
* Improved sample correlation annotation layout and value-label coloring.

### v2.5: Normalization memory workflow and heatmap layout refinement

Main changes:

* Removed permanent automatic normalization output folder generation.
* Loaded normalization outputs into Streamlit session memory for UI preview and manual CSV export.
* Further hardened QC grouping editor state to reduce browser autofill/widget-state contamination.
* Refined sample correlation heatmap layout, annotation bars, and correlation-value overlay behavior.

## License

No license has been selected yet.

This repository is currently private and intended for internal development.
::: 
