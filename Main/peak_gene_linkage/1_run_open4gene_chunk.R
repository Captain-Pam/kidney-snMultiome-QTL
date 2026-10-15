#!/usr/bin/env Rscript
# Purpose: Run original Open4Gene fasthurdle tests for one peak-gene chunk.
# Inputs: Chunk RNA/ATAC H5AD files, cell metadata, and precomputed TSS +/-250 kb pairs.
# Outputs: res_<chunk>.tsv and status_<chunk>.tsv in the output directory.
# Example: Rscript Main/peak_gene_linkage/1_run_open4gene_chunk.R \
#   --chunk-name chunk_001 --input-dir data --output-dir results --batches 20 \
#   --cores 20 --python /path/to/python

suppressPackageStartupMessages({
  library(Open4Gene)
  library(fasthurdle)
  library(argparser)
  library(data.table)
  library(Matrix)
  library(reticulate)
  library(parallel)
})

# Parse the chunk inputs and resource settings.
parser <- arg_parser("Run Open4Gene for one chunk")
parser <- add_argument(parser, "--chunk-name", help = "Chunk name", default = NA_character_)
parser <- add_argument(parser, "--input-dir", help = "Input directory", default = NA_character_)
parser <- add_argument(parser, "--output-dir", help = "Output directory", default = NA_character_)
parser <- add_argument(parser, "--python", help = "Python executable", default = "")
parser <- add_argument(parser, "--batches", help = "Batches per cell type", default = 20L, type = "integer")
parser <- add_argument(parser, "--cores", help = "Worker cores", default = 20L, type = "integer")
args <- parse_args(parser)
if (anyNA(c(args$chunk_name, args$input_dir, args$output_dir)) ||
    any(!nzchar(c(args$chunk_name, args$input_dir, args$output_dir)))) {
  stop("--chunk-name, --input-dir, and --output-dir are required")
}
if (is.na(args$batches) || args$batches < 1L || is.na(args$cores) || args$cores < 1L) {
  stop("--batches and --cores must be integers >= 1")
}

covariates <- c("log_nCount_RNA", "percent_mito", "log_nFrags")
min_gene_expr <- 0.01
min_peak_access <- 0.005
min_gene_cells <- 50
min_peak_cells <- 50

# Convert AnnData to feature-by-cell sparse matrices, preserving feature/cell names.
anndata_matrix <- function(adata) {
  x <- as(t(adata$X), "CsparseMatrix")
  rownames(x) <- as.character(adata$var_names$to_list())
  colnames(x) <- as.character(adata$obs_names$to_list())
  if (any(is.infinite(x@x))) {
    finite_values <- x@x[is.finite(x@x)]
    if (length(finite_values) == 0L) stop("Infinite matrix entries have no finite replacement")
    x@x[is.infinite(x@x)] <- max(finite_values)
  }
  x
}

# Apply detection fraction and minimum detected-cell filters within each cell type.
celltype_qc <- function(rna, atac) {
  gene_cells <- Matrix::rowSums(rna > 0)
  peak_cells <- Matrix::rowSums(atac > 0)
  list(
    rna = rna[gene_cells / ncol(rna) > min_gene_expr & gene_cells >= min_gene_cells, , drop = FALSE],
    atac = atac[peak_cells / ncol(atac) > min_peak_access & peak_cells >= min_peak_cells, , drop = FALSE]
  )
}

# Use the original package implementation without patches.
run_open4gene <- function(rna, atac, metadata, pairs, celltype) {
  obj <- CreateOpen4GeneObj(
    RNA = rna, ATAC = atac, Meta.data = metadata,
    Peak2Gene.Pairs = pairs, Covariates = covariates, Celltypes = "cell_state"
  )
  obj <- Open4Gene(object = obj, Celltype = celltype, Binary = FALSE,
                   Method = "fasthurdle", MinNum.Cells = 5)
  obj@Res
}

# Prepare outputs and remove a stale successful result before starting.
dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)
result_path <- file.path(args$output_dir, paste0("res_", args$chunk_name, ".tsv"))
status_path <- file.path(args$output_dir, paste0("status_", args$chunk_name, ".tsv"))
if (file.exists(result_path)) unlink(result_path)
empty_result <- data.table(Peak = character(), Gene = character(), Celltype = character(),
                           hurdle.Res.zero.p = numeric(), hurdle.Res.count.p = numeric())
empty_status <- data.table(chunk_name = args$chunk_name, Celltype = "", batch_id = 0L,
                           n_pairs = 0L, status = "success", error = "")

# The pair table already contains the TSS +/-250 kb candidate window.
pairs <- fread(file.path(args$input_dir, "peak_gene_pairs_tss_250kb_with_chunk.tsv"))
if (!all(c("Peak", "Gene", "chunk_name") %in% names(pairs))) {
  stop("Pair table requires Peak, Gene, and chunk_name")
}
pairs <- pairs[chunk_name == args$chunk_name, .(peak = Peak, gene = Gene)]
if (nrow(pairs) == 0L) {
  fwrite(empty_result, result_path, sep = "\t")
  fwrite(empty_status, status_path, sep = "\t")
  message("No candidate pairs for this chunk; wrote empty successful output")
  quit(status = 0L)
}

# Read the paired matrices and align metadata to their shared cells.
if (nzchar(args$python)) use_python(args$python, required = TRUE)
scanpy <- import("scanpy")
chunk_dir <- file.path(args$input_dir, "split_chunks_add_injPT", args$chunk_name)
message("Reading RNA, ATAC, and metadata for ", args$chunk_name)
rna <- anndata_matrix(scanpy$read_h5ad(file.path(chunk_dir, "snRNA_chunk.h5ad")))
atac <- anndata_matrix(scanpy$read_h5ad(file.path(chunk_dir, "snATAC_chunk.h5ad")))
metadata <- as.data.frame(fread(file.path(args$input_dir, "cell_rna_atac_metadata.tsv.gz")))
rownames(metadata) <- as.character(metadata[[1L]])
metadata <- metadata[, -1L, drop = FALSE]
if (!all(c("cell_state", "nCount_RNA", "percent_mito", "nFrags") %in% names(metadata))) {
  stop("Metadata requires cell_state, nCount_RNA, percent_mito, and nFrags")
}
metadata$log_nCount_RNA <- log1p(metadata$nCount_RNA)
metadata$log_nFrags <- log1p(metadata$nFrags)
shared_cells <- Reduce(intersect, list(colnames(rna), colnames(atac), rownames(metadata)))
rna <- rna[, shared_cells, drop = FALSE]
atac <- atac[, shared_cells, drop = FALSE]
metadata <- metadata[shared_cells, , drop = FALSE]
pairs <- pairs[peak %in% rownames(atac) & gene %in% rownames(rna)]

# Build flat tasks after independent cell-type QC.
celltypes <- unique(as.character(metadata$cell_state))
celltypes <- celltypes[!is.na(celltypes) & nzchar(trimws(celltypes))]
tasks <- list()
for (celltype in celltypes) {
  cells <- rownames(metadata)[!is.na(metadata$cell_state) & metadata$cell_state == celltype]
  qc <- celltype_qc(rna[, cells, drop = FALSE], atac[, cells, drop = FALSE])
  ct_pairs <- pairs[peak %in% rownames(qc$atac) & gene %in% rownames(qc$rna)]
  if (nrow(ct_pairs) == 0L) next
  n_batches <- min(args$batches, nrow(ct_pairs))
  batch_ids <- if (n_batches == 1L) rep(1L, nrow(ct_pairs)) else
    cut(seq_len(nrow(ct_pairs)), breaks = n_batches, labels = FALSE)
  for (batch_id in seq_len(n_batches)) {
    tasks[[length(tasks) + 1L]] <- list(
      celltype = celltype, batch_id = batch_id,
      pairs = as.data.frame(ct_pairs[batch_ids == batch_id, .(peak, gene)]),
      rna = qc$rna, atac = qc$atac, metadata = metadata[cells, , drop = FALSE]
    )
  }
}
if (length(tasks) == 0L) {
  fwrite(empty_result, result_path, sep = "\t")
  fwrite(empty_status, status_path, sep = "\t")
  message("No candidate pairs remain after feature matching and cell-type QC; wrote empty successful output")
  quit(status = 0L)
}

# Record every batch outcome and publish results only when all batches succeed.
message("Running ", length(tasks), " cell-type batches")
outcomes <- mclapply(tasks, function(task) {
  status <- data.table(chunk_name = args$chunk_name, Celltype = task$celltype,
                       batch_id = task$batch_id, n_pairs = nrow(task$pairs),
                       status = "success", error = "")
  tryCatch({
    result <- as.data.table(run_open4gene(task$rna, task$atac, task$metadata,
                                        task$pairs, task$celltype))
    if (!"Celltype" %in% names(result)) result[, Celltype := rep(task$celltype, .N)]
    list(result = result, status = status)
  }, error = function(e) {
    status[, `:=`(status = "failed", error = conditionMessage(e))]
    list(result = NULL, status = status)
  })
}, mc.cores = min(args$cores, length(tasks)))
statuses <- rbindlist(lapply(outcomes, `[[`, "status"))
fwrite(statuses, status_path, sep = "\t")
if (any(statuses$status == "failed")) {
  if (file.exists(result_path)) unlink(result_path)
  message("One or more batches failed; see ", status_path)
  quit(status = 1L)
}
results <- rbindlist(lapply(outcomes, `[[`, "result"), use.names = TRUE, fill = TRUE)
if (ncol(results) == 0L) results <- empty_result
fwrite(results, result_path, sep = "\t")
message("Completed chunk ", args$chunk_name, "; wrote ", result_path)
