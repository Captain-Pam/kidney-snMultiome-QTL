#!/usr/bin/env Rscript
#
# Purpose:
#   Extract the Gene Expression feature type from a Cell Ranger ARC multiome
#   matrix to create a gene-expression-only input for CellBender.
#
# Input:
#   A Cell Ranger ARC raw_feature_bc_matrix directory containing
#   barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz.
#
# Output:
#   A gene-expression-only 10x matrix directory containing
#   barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz.

# Set the input and output directories for one sample.
input_dir <- "/path/to/cellranger_arc/outs/raw_feature_bc_matrix"
output_dir <- "/path/to/output/gex_matrix"

# Check the required Cell Ranger matrix files.
required_files <- c("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz")
input_files <- file.path(input_dir, required_files)
if (!all(file.exists(input_files))) {
  stop("Missing required input files in: ", input_dir)
}

message("Reading the Cell Ranger ARC matrix.")

barcodes <- data.table::fread(
  file.path(input_dir, "barcodes.tsv.gz"),
  header = FALSE
)
features <- data.table::fread(
  file.path(input_dir, "features.tsv.gz"),
  header = FALSE,
  sep = "\t"
)

if (ncol(features) < 3L) {
  stop("features.tsv.gz must contain a feature-type column.")
}

matrix_file <- file.path(input_dir, "matrix.mtx.gz")
matrix_connection <- gzfile(matrix_file, open = "rt")
matrix_header <- readLines(matrix_connection, n = 3L)
close(matrix_connection)

if (length(matrix_header) != 3L) {
  stop("matrix.mtx.gz does not contain a valid three-line header.")
}

matrix_entries <- data.table::fread(
  matrix_file,
  skip = 3L,
  header = FALSE,
  col.names = c("feature_index", "barcode_index", "count")
)

# Retain Gene Expression features and remap their matrix row indices.
gex_feature_indices <- which(features[[3L]] == "Gene Expression")
if (length(gex_feature_indices) == 0L) {
  stop("No Gene Expression features were found.")
}

gex_features <- features[gex_feature_indices, 1:3, with = FALSE]
gex_row_map <- integer(nrow(features))
gex_row_map[gex_feature_indices] <- seq_along(gex_feature_indices)

gex_matrix <- matrix_entries[feature_index %in% gex_feature_indices]
gex_matrix[, feature_index := gex_row_map[feature_index]]

# Update the Matrix Market dimensions for the filtered matrix.
matrix_header[[3L]] <- sprintf(
  "%d %d %d",
  nrow(gex_features),
  nrow(barcodes),
  nrow(gex_matrix)
)

# Write the gene-expression-only matrix in 10x format.
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

output_matrix <- file.path(output_dir, "matrix.mtx")
output_matrix_gz <- paste0(output_matrix, ".gz")
if (file.exists(output_matrix)) {
  file.remove(output_matrix)
}
if (file.exists(output_matrix_gz)) {
  file.remove(output_matrix_gz)
}

writeLines(matrix_header, output_matrix)
data.table::fwrite(
  gex_matrix,
  output_matrix,
  append = TRUE,
  sep = " ",
  col.names = FALSE
)

gzip_bin <- Sys.which("gzip")
if (gzip_bin == "") {
  stop("gzip was not found on PATH.")
}
gzip_status <- system2(gzip_bin, c("-f", shQuote(output_matrix)))
if (gzip_status != 0L) {
  stop("gzip failed while compressing matrix.mtx.")
}

data.table::fwrite(
  barcodes,
  file.path(output_dir, "barcodes.tsv.gz"),
  sep = "\t",
  col.names = FALSE,
  compress = "gzip"
)
data.table::fwrite(
  gex_features,
  file.path(output_dir, "features.tsv.gz"),
  sep = "\t",
  col.names = FALSE,
  compress = "gzip"
)

message("Gene Expression matrix extraction completed.")
