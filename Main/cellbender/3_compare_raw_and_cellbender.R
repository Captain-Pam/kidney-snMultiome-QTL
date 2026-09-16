#!/usr/bin/env Rscript
#
# Purpose:
#   Compare extracted Gene Expression counts before and after CellBender,
#   summarize per-sample changes, and identify features most affected by
#   ambient-RNA removal.
#
# Input:
#   1. A root directory containing one CellBender *_filtered.h5 file per
#      sample, stored within a sample-named subdirectory.
#   2. A root directory containing one extracted 10x Gene Expression matrix
#      directory per sample.
#
# Output:
#   A dual-assay Seurat object, per-sample median QC statistics, a table of
#   feature-level count changes, and a feature-comparison plot.

suppressPackageStartupMessages({
  library(Seurat)
  library(scCustomize)
})

# Set the CellBender, extracted Gene Expression, and QC output directories.
cellbender_dir <- "/path/to/cellbender_outputs"
raw_gex_dir <- "/path/to/extracted_gex_matrices"
output_dir <- "/path/to/cellbender_qc"

# Retain the filtering and plotting thresholds used in the analysis.
min_cells <- 3L
min_features <- 200L
pct_diff_thresholds <- c(25, 50)

if (!dir.exists(cellbender_dir)) {
  stop("CellBender input directory does not exist: ", cellbender_dir)
}
if (!dir.exists(raw_gex_dir)) {
  stop("Raw Gene Expression input directory does not exist: ", raw_gex_dir)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Identify CellBender files and use their parent directory names as sample IDs.
cellbender_files <- list.files(
  cellbender_dir,
  pattern = "_filtered\\.h5$",
  full.names = TRUE,
  recursive = TRUE
)
if (length(cellbender_files) == 0L) {
  stop("No CellBender *_filtered.h5 files were found in: ", cellbender_dir)
}

cellbender_samples <- basename(dirname(cellbender_files))
if (anyDuplicated(cellbender_samples)) {
  duplicated_samples <- unique(cellbender_samples[duplicated(cellbender_samples)])
  stop(
    "Multiple CellBender files were found for sample(s): ",
    paste(duplicated_samples, collapse = ", ")
  )
}
names(cellbender_files) <- cellbender_samples

# Identify sample directories containing complete extracted 10x matrices.
raw_sample_dirs <- list.dirs(raw_gex_dir, full.names = TRUE, recursive = FALSE)
required_10x_files <- c("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz")
has_complete_matrix <- vapply(
  raw_sample_dirs,
  function(sample_dir) {
    all(file.exists(file.path(sample_dir, required_10x_files)))
  },
  logical(1)
)
raw_sample_dirs <- raw_sample_dirs[has_complete_matrix]
if (length(raw_sample_dirs) == 0L) {
  stop("No complete extracted 10x matrices were found in: ", raw_gex_dir)
}

raw_samples <- basename(raw_sample_dirs)
if (anyDuplicated(raw_samples)) {
  stop("Extracted Gene Expression sample directory names must be unique.")
}
names(raw_sample_dirs) <- raw_samples

# Match samples available before and after CellBender processing.
sample_names <- sort(intersect(cellbender_samples, raw_samples))
if (length(sample_names) == 0L) {
  stop("No matching sample names were found between the two input directories.")
}

unmatched_cellbender <- setdiff(cellbender_samples, sample_names)
unmatched_raw <- setdiff(raw_samples, sample_names)
if (length(unmatched_cellbender) > 0L) {
  warning(
    "Skipping CellBender sample(s) without extracted matrices: ",
    paste(sort(unmatched_cellbender), collapse = ", ")
  )
}
if (length(unmatched_raw) > 0L) {
  warning(
    "Skipping extracted matrix sample(s) without CellBender output: ",
    paste(sort(unmatched_raw), collapse = ", ")
  )
}

cellbender_files <- cellbender_files[sample_names]
raw_sample_dirs <- raw_sample_dirs[sample_names]
message("Comparing ", length(sample_names), " matched sample(s).")

# Read and merge the CellBender-filtered matrices with sample-prefixed barcodes.
message("Reading CellBender-filtered matrices.")
cellbender_matrices <- lapply(
  cellbender_files,
  scCustomize::Read_CellBender_h5_Mat
)
cellbender_merged <- scCustomize::Merge_Sparse_Data_All(
  matrix_list = cellbender_matrices,
  add_cell_ids = sample_names
)

# Read the extracted pre-CellBender Gene Expression matrices.
message("Reading extracted Gene Expression matrices.")
raw_matrices <- lapply(raw_sample_dirs, function(sample_dir) {
  counts <- Seurat::Read10X(data.dir = sample_dir)
  if (is.list(counts)) {
    if (!"Gene Expression" %in% names(counts)) {
      stop("Gene Expression data were not found in: ", sample_dir)
    }
    counts <- counts[["Gene Expression"]]
  }
  counts
})
raw_merged <- scCustomize::Merge_Sparse_Data_All(
  matrix_list = raw_matrices,
  add_cell_ids = sample_names
)

# Preserve an explicit sample assignment for every prefixed raw-data barcode.
sample_map <- Map(
  function(counts, sample_name) {
    merged_barcodes <- paste(sample_name, colnames(counts), sep = "_")
    stats::setNames(rep(sample_name, length(merged_barcodes)), merged_barcodes)
  },
  raw_matrices,
  sample_names
)
cell_to_sample <- do.call(c, unname(sample_map))

# Create a dual-assay object after applying the original feature and cell filters.
message("Creating the dual-assay Seurat object.")
dual_seurat <- scCustomize::Create_CellBender_Merged_Seurat(
  raw_cell_bender_matrix = cellbender_merged,
  raw_counts_matrix = raw_merged,
  raw_assay_name = "RAW",
  min_cells = min_cells,
  min_features = min_features
)

sample_assignment <- unname(cell_to_sample[colnames(dual_seurat)])
if (anyNA(sample_assignment)) {
  stop("Sample assignments could not be recovered for all retained barcodes.")
}
dual_seurat$orig.ident <- sample_assignment

# Quantify CellBender-associated changes per cell and summarize them per sample.
message("Calculating per-cell and per-sample count differences.")
dual_seurat <- scCustomize::Add_CellBender_Diff(
  seurat_object = dual_seurat,
  raw_assay_name = "RAW",
  cell_bender_assay_name = "RNA"
)
median_stats <- scCustomize::Median_Stats(
  seurat_object = dual_seurat,
  group_by_var = "orig.ident",
  default_var = TRUE,
  median_var = c("nCount_Diff", "nFeature_Diff")
)

qs::qsave(
  dual_seurat,
  file = file.path(output_dir, "dual_seurat.qs"),
  nthreads = 2
)
utils::write.csv(
  median_stats,
  file = file.path(output_dir, "Raw_counts_vs_cellbender_median_stats.csv"),
  quote = FALSE,
  row.names = FALSE
)

# Summarize feature-level count changes across all matched samples.
message("Calculating feature-level count differences.")
feature_diff <- scCustomize::CellBender_Feature_Diff(
  seurat_object = dual_seurat,
  raw_assay = "RAW",
  cell_bender_assay = "RNA"
)
feature_diff_output <- data.frame(
  Feature = rownames(feature_diff),
  feature_diff,
  row.names = NULL,
  check.names = FALSE
)
utils::write.csv(
  feature_diff_output,
  file = file.path(
    output_dir,
    "Raw_counts_vs_cellbender_most_changed_features.csv"
  ),
  quote = FALSE,
  row.names = FALSE
)

# Plot feature-level comparisons at the two thresholds used in the analysis.
comparison_plots <- lapply(pct_diff_thresholds, function(threshold) {
  scCustomize::CellBender_Diff_Plot(
    feature_diff_df = feature_diff,
    pct_diff_threshold = threshold
  )
})
combined_plot <- patchwork::wrap_plots(comparison_plots, ncol = 2)
ggplot2::ggsave(
  filename = file.path(output_dir, "Raw_counts_vs_cellbender.jpg"),
  plot = combined_plot,
  width = 15,
  height = 10,
  dpi = 300
)

message("CellBender comparison completed.")
