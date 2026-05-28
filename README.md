# Bulk RNA-seq Explorer

A Python/Streamlit-based MVP for interactive bulk RNA-seq count matrix exploration, gene identifier handling, quality control, sample grouping, and future downstream RNA-seq analysis.

This project is currently in active prototype development. It was migrated from an earlier browser-only HTML prototype into a more reproducible Python/Streamlit architecture.

## Current Status

Current app version: **bulk_rnaseq_explorer_v1_10**

Current development stage: **early Streamlit MVP / active internal prototype**

The app currently supports:

* Uploading a raw bulk RNA-seq count matrix
* Detecting sample columns from the uploaded matrix
* Detecting whether gene identifiers are mouse Ensembl IDs, gene symbols, mixed, or unclear
* Converting mouse Ensembl IDs to gene symbols using local mapping resources
* Merging duplicated processed gene symbols by summing raw counts
* Creating a processed count matrix for downstream analysis
* Computing sample-level QC metrics
* Displaying a clean Quality Control workflow
* Creating reusable QC grouping sets
* Generating configurable QC bar plots:

  * Library Size
  * Detected Genes
  * Zero-count Fraction
* Plotting by individual samples or saved QC groups
* Supporting group-level aggregation with overlaid sample dots
* Customizing plot width, height, x-axis angle, axis-title font size, and colors
* Exporting QC plots as PNG and true SVG through Plotly/Kaleido
* Tracking development through Git/GitHub

Planned future features include:

* Normalization module
* PCA and sample correlation
* DESeq2-based differential expression analysis through an R backend
* Volcano plots
* Heatmaps
* ORA pathway analysis
* GSEA / fgsea-style pathway analysis
* Exportable result tables and reports
* Modular project structure for future integration with single-cell analysis tools

## Background

The original browser-only HTML prototype was developed as a fast, interactive proof-of-concept for bulk RNA-seq analysis. It included count matrix upload/parsing, sample grouping, gene mapping, QC, normalization, exploratory DEG-style analysis, volcano plots, heatmaps, ORA, and GSEA-style pathway analysis.

The current Python/Streamlit version treats the HTML app as a validated product reference, but rebuilds the backend using more reproducible Python/R-compatible logic.

The long-term goal is to build a scientific analysis platform that can eventually integrate:

* Bulk RNA-seq analysis
* Single-cell RNA-seq analysis
* Spatial/flow/image analysis modules
* AI-assisted scientific workflow organization
* Local and potentially cloud/server-based analysis execution

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
├── bulk_rnaseq_explorer.py
├── index_v5.4.9.html
├── .gitignore
└── README.md
```

The old `bulk_rnaseq_explorer_v1_0.py` to `bulk_rnaseq_explorer_v1_5.py` files are preserved in `archive/` as early migration history.

The active app entry point is now:

```text
bulk_rnaseq_explorer.py
```

Future intended structure:

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

This modular structure has not yet been implemented. The project is intentionally still kept as a single main Streamlit file while the core workflow is being stabilized.

## Installation

Recommended environment: Python with Streamlit and common data analysis packages.

Minimum required packages:

```bash
pip install streamlit pandas plotly kaleido
```

Optional package for faster cached gene mapping:

```bash
pip install pyarrow
```

`kaleido` is required for static PNG/SVG export from Plotly.

## Running the App

Run the current active app:

```bash
streamlit run bulk_rnaseq_explorer.py
```

To reduce Streamlit toolbar/menu visibility, optionally create:

```text
.streamlit/config.toml
```

with:

```toml
[client]
toolbarMode = "minimal"
```

## Input Format

The app expects a raw bulk RNA-seq count matrix.

Recommended format: tab-delimited `.tsv` or `.txt`

CSV is also accepted.

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
* Publication-grade downstream differential expression analysis should use integer raw counts

## Gene Identifier Handling

The app currently supports:

* Mouse Ensembl ID detection
* Mouse gene symbol detection
* Mixed/unclear gene ID warning
* Local mouse Ensembl-to-symbol conversion
* Duplicate gene symbol detection
* Duplicate gene merging by summing raw counts

Gene mapping resources are loaded from local mapping files when available.

A generated `.parquet` cache may be created for faster future loading. This cache is treated as a generated file and should usually be ignored by Git.

## Quality Control Workflow

The current Quality Control workflow includes:

### Dataset Summary

The Dataset Summary table currently reports:

* Sample
* Library size
* Detected genes
* Zero-count genes
* Zero fraction
* Mean count

### QC Grouping

Users can create named QC grouping sets.

Current behavior:

* Users can assign samples into groups
* A sample should only belong to one group within a grouping set
* Saved grouping sets can be selected for group-level QC plotting
* Saving a grouping set resets the editor back to the default empty grouping state
* Grouping placeholder names are shown as placeholders rather than pre-filled values

### QC Bar Plots

Current QC bar plots:

* Library Size
* Detected Genes
* Zero-count Fraction

Each plot supports:

* Plot by sample name
* Plot by QC assignment group
* Mean / median / sum aggregation where appropriate
* Sample-dot overlay on grouped bar plots
* Adjustable plot width
* Adjustable plot height
* Adjustable x-axis angle
* Adjustable axis-title font size
* Custom sample/group colors
* Reset to default plot settings
* PNG export
* True SVG export through Plotly/Kaleido

SVG export should use:

```python
fig.to_image(format="svg")
```

and should not use screenshot/canvas-based pseudo-SVG export.

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

Refactored toward the original HTML prototype workflow.

Main changes:

* Removed mandatory metadata upload
* Switched to single raw count matrix upload
* Detected sample columns directly from the count matrix
* Added early in-app sample grouping
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
* Fixed stale state when uploaded file was removed
* Prepared app for QC visualization

### v1.5

Added early Quality Control plotting workflow.

Main changes:

* Renamed QC Overview to Quality Control
* Moved QC grouping into the Quality Control workflow
* Added QC grouping set concept
* Added interactive QC bar plots:

  * Library Size
  * Detected Genes
  * Zero-count Fraction
* Began using Plotly-based visualization

### v1.6

Migrated more of the original HTML QC barplot behavior.

Main changes:

* Improved QC barplot controls
* Added Plot by Sample name / QC assignment group
* Disabled QC assignment set and aggregation controls when plotting by sample
* Added group-level aggregation
* Added sample dots over group-level bars
* Added advanced plot settings
* Added color settings
* Added PNG/SVG export logic
* Added per-plot reset behavior
* Improved QC grouping draft handling

### v1.7

Improved QC grouping and plot performance behavior.

Main changes:

* Removed unnecessary Apply grouping draft button
* Reorganized QC grouping editor buttons
* Reduced unnecessary QC recomputation
* Added cached QC summary behavior
* Added cached QC plot data behavior
* Improved plot setting responsiveness
* Removed Plotly legend from QC barplots
* Improved reset/export separation

### v1.8

Stabilized QC grouping editor behavior.

Main changes:

* Removed `Median count` from the displayed Dataset Summary table
* Moved Add group into the QC grouping editor card
* Added Clear grouping info
* Changed grouping editor from form-based delayed updates to immediate session-state updates
* Added stable group IDs
* Prevented already-assigned samples from appearing in other group options
* Preserved saved QC grouping sets for plotting
* Improved single-row button layout

### v1.9

Improved QC barplot reset/export/color behavior.

Main changes:

* Added axis-title font-size setting
* Made color settings more compact
* Improved plot reset cleanup
* Added direct PNG/SVG download buttons
* Kept SVG export based on Plotly/Kaleido
* Improved download button layout
* Added deterministic visual export cache keys

### v1.10

Focused on QC grouping reset behavior and robust barplot state handling.

Main changes:

* Save QC grouping now saves the grouping set and then resets the grouping editor
* Grouping set names and group names now use placeholders instead of pre-filled example text
* Axis title font size is intended to control both x-axis and y-axis title font size
* Added stronger reset troubleshooting strategy using plot-specific reset nonce / widget-key versioning
* Improved cleanup of old Streamlit widget keys after reset
* Improved cleanup of plot-specific export and plot-data caches
* Centered Download PNG/SVG button text
* Preserved true SVG export through Plotly/Kaleido

## Development Notes

This repository is currently private and used for active development.

The early versioned files are preserved as a migration archive. Current development should happen in:

```text
bulk_rnaseq_explorer.py
```

Going forward, the preferred workflow is:

* Keep one active app entry file
* Let Git track version history
* Avoid creating new `v1_11.py`, `v1_12.py`, etc.
* Use clear commit messages for each functional update
* Refactor into modules only after the core workflow stabilizes

## Git Workflow

Typical workflow:

```bash
git status
git add bulk_rnaseq_explorer.py
git commit -m "Describe the update"
git push
```

For work across multiple computers:

```bash
git pull
```

should be run before starting new work.

Avoid committing generated cache files such as:

```text
assets/*.parquet
```

unless there is a specific reason to version them.

## Current Near-Term Roadmap

Next likely development steps:

1. Finish QC barplot interaction polish
2. Add normalization setup
3. Add normalized matrix generation
4. Add PCA
5. Add sample correlation
6. Add DEG setup
7. Add R/DESeq2 backend runner
8. Add volcano plot
9. Add heatmap
10. Add ORA/GSEA pathway analysis

## License

No license has been selected yet.

This repository is currently private and intended for internal development.
