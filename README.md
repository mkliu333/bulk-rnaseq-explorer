# Bulk RNA-seq Explorer

A Python/Streamlit-based MVP for interactive bulk RNA-seq count matrix exploration, quality control, sample grouping, gene identifier handling, and future downstream differential expression/pathway analysis.

This project is currently in active prototype development. It was migrated from an earlier browser-only HTML prototype into a more reproducible Python/Streamlit architecture.

## Project Status

Current stage: **early Streamlit MVP / active development**

The current prototype supports:

* Uploading a raw bulk RNA-seq count matrix
* Detecting sample columns from the count matrix
* Detecting whether gene identifiers are Ensembl IDs or gene symbols
* Converting mouse Ensembl IDs to gene symbols using local mapping resources
* Merging duplicated gene symbols by summing raw counts
* Creating a processed count matrix for downstream analysis
* Computing basic sample-level QC metrics
* Building early Quality Control visualizations
* Managing versioned development through Git

Planned features include:

* Normalization module
* PCA and sample correlation
* DESeq2-based differential expression analysis through R backend
* Volcano plot and heatmap visualization
* ORA and GSEA pathway analysis
* Exportable reports and result tables
* Modular project structure for future integration with single-cell analysis tools

## Background

The original browser-only prototype was developed as a fast, interactive proof-of-concept for bulk RNA-seq analysis. It included upload/parsing, sample grouping, gene mapping, QC, normalization, DEG-style exploratory analysis, volcano/heatmap visualization, ORA, and GSEA-style pathway analysis.

The current Python/Streamlit version treats the HTML prototype as a product and UI reference, but rebuilds the backend using more reproducible Python/R-compatible logic.

The long-term goal is to develop a reliable scientific analysis platform that can eventually integrate:

* Bulk RNA-seq analysis
* Single-cell RNA-seq analysis
* Spatial/flow/image analysis modules
* AI-assisted scientific workflow organization
* Local and potentially cloud/server-based analysis execution

## Current Repository Structure

```text
bulk-rnaseq-explorer/
├── assets/
│   └── generated cache files, ignored by Git when appropriate
├── database_raw/
│   └── local GMT pathway databases
├── source_mapping/
│   └── mouse Ensembl-to-gene-symbol mapping resources
├── bulk_rnaseq_explorer_v1_0.py
├── bulk_rnaseq_explorer_v1_1.py
├── bulk_rnaseq_explorer_v1_2.py
├── bulk_rnaseq_explorer_v1_3.py
├── bulk_rnaseq_explorer_v1_4.py
├── bulk_rnaseq_explorer_v1_5.py
├── index_v5.4.9.html
├── .gitignore
└── README.md
```

The multiple versioned `.py` files are intentionally preserved during the prototype stage to make early migration history easy to inspect. Later, the project should be refactored into a cleaner modular structure with a single app entry point.

Planned future structure:

```text
bulk-rnaseq-explorer/
├── app.py
├── core/
│   ├── io.py
│   ├── state.py
│   ├── gene_mapping.py
│   └── validation.py
├── analysis/
│   ├── qc.py
│   ├── normalization.py
│   └── deg.py
├── plotting/
│   ├── qc_plots.py
│   ├── pca.py
│   └── heatmap.py
├── r_scripts/
│   ├── run_deseq2.R
│   ├── run_fgsea.R
│   └── run_ora.R
├── assets/
├── database_raw/
├── source_mapping/
├── outputs/
├── requirements.txt
└── README.md
```

## Version History

### v1.0

Initial Streamlit skeleton.

Main features:

* Basic Streamlit app layout
* Count matrix upload
* Metadata upload
* Optional gene map upload
* Basic input validation
* Sidebar workflow structure

### v1.1

Refactored the workflow toward the HTML prototype.

Main changes:

* Removed mandatory metadata workflow
* Switched to single raw count matrix upload
* Detected sample columns directly from the count matrix
* Added in-app sample grouping
* Improved count matrix validation

### v1.2

Added migration planning and validation layer.

Main changes:

* Added legacy HTML migration map
* Treated the HTML prototype as a product specification
* Improved group assignment validation
* Added QC readiness summary
* Improved local gene map detection

### v1.3

Cleaned product UI and added gene identifier processing.

Main changes:

* Removed developer-facing validation and migration tabs from the main UI
* Added gene ID mode detection
* Added Ensembl ID to gene symbol conversion
* Added duplicate gene symbol detection
* Added duplicate merging by summing raw counts
* Generated processed count matrix for downstream QC

### v1.4

Added QC data summary layer.

Main changes:

* Simplified sidebar
* Removed stale project status display
* Added sample-level QC table
* Added basic dataset summary from processed counts
* Fixed stale state when uploaded file is removed
* Prepared app for QC visualization

### v1.5

Added early Quality Control plotting workflow.

Main changes:

* Renamed QC Overview to Quality Control
* Moved QC grouping into the Quality Control workflow
* Added QC grouping set concept
* Added interactive QC bar plots:

  * Library size
  * Detected genes
  * Zero-count fraction
* Began using Plotly-based visualization

## Installation

Recommended environment: Python with Streamlit and common data analysis packages.

Install required Python packages:

```bash
pip install streamlit pandas plotly kaleido
```

Optional package for faster cached gene mapping:

```bash
pip install pyarrow
```

## Running the App

Run a specific prototype version:

```bash
streamlit run bulk_rnaseq_explorer_v1_5.py
```

Future versions may consolidate the app into:

```bash
streamlit run app.py
```

## Input Format

The app expects a raw bulk RNA-seq count matrix.

Recommended format: tab-delimited `.tsv` or `.txt`

Example:

```text
EnsemblID / Gene_symbol    Sample_1    Sample_2    Sample_3    Sample_4
ENSMUSG00000000001         120         98          115         130
ENSMUSG00000000028         0           4           1           3
Cxcl1                      50          80          320         400
Actb                       10000       9800        10300       9900
```

Rules:

* First row should contain sample names
* First column should contain Ensembl IDs or gene symbols
* Count values should be raw counts
* DESeq2-compatible downstream analysis will require integer raw counts

## Gene Identifier Handling

The app currently supports:

* Mouse Ensembl ID detection
* Mouse gene symbol detection
* Local Ensembl-to-symbol conversion
* Duplicate gene symbol detection
* Duplicate gene merging by summing raw counts

Gene mapping resources are loaded from local mapping files when available.

A generated `.parquet` cache may be created for faster future loading. This cache is treated as a generated file and may be ignored by Git.

## Development Notes

This repository is currently private and used for active development. The early versioned files are kept intentionally to preserve the migration record from HTML prototype to Streamlit MVP.

As the project matures, development should move toward:

* One main app entry point
* Modular Python files
* Reusable backend analysis functions
* Clear input/output contracts
* R backend integration for publication-grade DESeq2, fgsea, and ORA analyses
* Cleaner deployment-ready repo structure

## Git Workflow

Typical workflow:

```bash
git status
git add <changed_files>
git commit -m "Describe the update"
git push
```

For work across multiple computers:

```bash
git pull
```

should be run before starting new work.

## License

No license has been selected yet. This repository is currently private and intended for internal development.
