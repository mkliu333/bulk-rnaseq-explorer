#!/usr/bin/env Rscript

parse_args <- function(args) {
  parsed <- list()
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) {
      stop(paste("Unexpected argument:", key), call. = FALSE)
    }
    if (i == length(args)) {
      stop(paste("Missing value for argument:", key), call. = FALSE)
    }
    parsed[[substring(key, 3)]] <- args[[i + 1]]
    i <- i + 2
  }
  parsed
}

require_package <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(paste("Missing R package:", pkg), call. = FALSE)
  }
}

write_matrix_csv <- function(matrix_obj, path, gene_col = "Gene") {
  output_df <- data.frame(Gene = rownames(matrix_obj), matrix_obj, check.names = FALSE)
  names(output_df)[1] <- gene_col
  write.csv(output_df, path, row.names = FALSE, quote = FALSE)
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
counts_path <- args[["counts"]]
outdir <- args[["outdir"]]
gene_col <- ifelse(is.null(args[["gene_col"]]), "Gene", args[["gene_col"]])
prior_count <- ifelse(is.null(args[["prior_count"]]), 1, as.numeric(args[["prior_count"]]))

if (is.null(counts_path) || is.null(outdir)) {
  stop("Required arguments: --counts and --outdir", call. = FALSE)
}

require_package("DESeq2")
require_package("edgeR")
require_package("jsonlite")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

counts_df <- read.csv(counts_path, check.names = FALSE, stringsAsFactors = FALSE)
if (!(gene_col %in% names(counts_df))) {
  stop(paste("Gene column not found:", gene_col), call. = FALSE)
}
if (ncol(counts_df) < 3) {
  stop("At least two sample columns are required.", call. = FALSE)
}

gene_ids <- as.character(counts_df[[gene_col]])
sample_names <- names(counts_df)[names(counts_df) != gene_col]
count_values <- counts_df[, sample_names, drop = FALSE]
count_matrix <- as.matrix(sapply(count_values, as.numeric))
if (ncol(count_matrix) == 1) {
  count_matrix <- matrix(count_matrix, ncol = 1)
  colnames(count_matrix) <- sample_names
}
rownames(count_matrix) <- gene_ids
colnames(count_matrix) <- sample_names

if (any(is.na(count_matrix))) {
  stop("Non-numeric count values detected.", call. = FALSE)
}
if (any(count_matrix < 0)) {
  stop("Negative count values detected.", call. = FALSE)
}
if (any(is.na(gene_ids) | trimws(gene_ids) == "")) {
  stop("Empty gene IDs detected.", call. = FALSE)
}

original_genes <- nrow(count_matrix)
count_matrix <- rowsum(count_matrix, group = gene_ids, reorder = FALSE)
genes_after_duplicate_merge <- nrow(count_matrix)
write_matrix_csv(count_matrix, file.path(outdir, "raw_counts.csv"), gene_col)

deviation <- abs(count_matrix - round(count_matrix))
max_deviation <- ifelse(length(deviation) > 0, max(deviation), 0)
non_integer_fraction <- ifelse(length(deviation) > 0, sum(deviation > 1e-6) / length(deviation), 0)
rounding_applied <- max_deviation > 1e-6 || non_integer_fraction > 0

all_zero_mask <- rowSums(count_matrix) == 0
all_zero_genes <- sum(all_zero_mask)
filtered_counts <- count_matrix[!all_zero_mask, , drop = FALSE]
if (nrow(filtered_counts) == 0) {
  stop("No non-zero genes remain after all-zero gene filtering.", call. = FALSE)
}
rounded_counts <- round(filtered_counts)
rounded_counts[rounded_counts < 0] <- 0

y_raw <- edgeR::DGEList(counts = rounded_counts)
cpm_raw <- edgeR::cpm(y_raw, normalized.lib.sizes = FALSE, log = FALSE)
log2_cpm_plus1 <- log2(cpm_raw + 1)
write_matrix_csv(cpm_raw, file.path(outdir, "cpm.csv"), gene_col)
write_matrix_csv(log2_cpm_plus1, file.path(outdir, "log2_cpm_plus1.csv"), gene_col)

col_data <- data.frame(row.names = colnames(rounded_counts), condition = rep("all", ncol(rounded_counts)))
dds <- DESeq2::DESeqDataSetFromMatrix(countData = rounded_counts, colData = col_data, design = ~1)
dds <- DESeq2::estimateSizeFactors(dds)
deseq2_norm_counts <- DESeq2::counts(dds, normalized = TRUE)
vsd <- tryCatch(
  DESeq2::vst(dds, blind = TRUE),
  error = function(vst_error) {
    DESeq2::varianceStabilizingTransformation(dds, blind = TRUE)
  }
)
deseq2_vst <- SummarizedExperiment::assay(vsd)
size_factors <- data.frame(
  Sample = names(DESeq2::sizeFactors(dds)),
  SizeFactor = as.numeric(DESeq2::sizeFactors(dds)),
  check.names = FALSE
)
write.csv(size_factors, file.path(outdir, "deseq2_size_factors.csv"), row.names = FALSE, quote = FALSE)
write_matrix_csv(deseq2_norm_counts, file.path(outdir, "deseq2_normalized_counts.csv"), gene_col)
write_matrix_csv(deseq2_vst, file.path(outdir, "deseq2_vst.csv"), gene_col)

y_tmm <- edgeR::DGEList(counts = rounded_counts)
y_tmm <- edgeR::calcNormFactors(y_tmm, method = "TMM")
tmm_cpm <- edgeR::cpm(y_tmm, normalized.lib.sizes = TRUE, log = FALSE)
tmm_logcpm <- edgeR::cpm(y_tmm, normalized.lib.sizes = TRUE, log = TRUE, prior.count = prior_count)
tmm_factors <- data.frame(
  Sample = colnames(rounded_counts),
  LibSize = as.numeric(y_tmm$samples$lib.size),
  NormFactor = as.numeric(y_tmm$samples$norm.factors),
  EffectiveLibSize = as.numeric(y_tmm$samples$lib.size * y_tmm$samples$norm.factors),
  check.names = FALSE
)
write.csv(tmm_factors, file.path(outdir, "edger_tmm_norm_factors.csv"), row.names = FALSE, quote = FALSE)
write_matrix_csv(tmm_cpm, file.path(outdir, "edger_tmm_cpm.csv"), gene_col)
write_matrix_csv(tmm_logcpm, file.path(outdir, "edger_tmm_logcpm.csv"), gene_col)

report <- list(
  original_genes = as.integer(original_genes),
  genes_after_duplicate_merge = as.integer(genes_after_duplicate_merge),
  all_zero_genes = as.integer(all_zero_genes),
  genes_used_after_zero_filtering = as.integer(nrow(filtered_counts)),
  samples = as.integer(ncol(filtered_counts)),
  sample_names = sample_names,
  prior_count = prior_count,
  timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"),
  r_version = R.version.string,
  deseq2_version = as.character(utils::packageVersion("DESeq2")),
  edger_version = as.character(utils::packageVersion("edgeR")),
  rounding_applied = rounding_applied,
  max_deviation_from_integer = max_deviation,
  non_integer_value_fraction = non_integer_fraction
)
jsonlite::write_json(report, file.path(outdir, "normalization_report.json"), pretty = TRUE, auto_unbox = TRUE)

cat("Normalization completed successfully.\n")
