#!/usr/bin/env Rscript
# Purpose: Build donor-level SCT objects for all selected cell types.
# Inputs:
#   - 1_stage1_weights/2_data/celltypes.tsv with counts file names.
#   - Per-cell-type raw-count pseudobulk CSV files, genes x donors.
#   - QCed WGS PLINK .fam file for genotyped donor IDs.
# Outputs:
#   - 1_stage1_weights/2_data/sct_objects/<celltype>.sct.rds
#   - 1_stage1_weights/3_outputs/summaries/sct_objects/<celltype>.summary.tsv

suppressMessages({
  library(argparser)
  library(Seurat)
  library(Matrix)
})

p <- arg_parser("Build SCT objects for scTWAS Stage 1.")
p <- add_argument(p, "--celltype", default = "all",
                  help = "Cell type name from celltypes.tsv, or 'all'.")
p <- add_argument(p, "--workflow-dir", default = "",
                  help = "Path to the public Main/sctwas directory.")
p <- add_argument(p, "--counts-dir", default = "",
                  help = "Directory containing per-cell-type pseudobulk count CSV files.")
p <- add_argument(p, "--fam", default = "",
                  help = "QCed genotype PLINK .fam file.")
argv <- parse_args(p)

if (argv$workflow_dir == "") stop("--workflow-dir is required.")
if (argv$counts_dir == "") stop("--counts-dir is required.")
if (argv$fam == "") stop("--fam is required.")

stage_dir <- file.path(argv$workflow_dir, "1_stage1_weights")
manifest_path <- file.path(stage_dir, "2_data/celltypes.tsv")
out_dir <- file.path(stage_dir, "2_data/sct_objects")
summary_dir <- file.path(stage_dir, "3_outputs/summaries/sct_objects")

if (!file.exists(manifest_path)) stop("Missing cell type manifest: ", manifest_path)
if (!file.exists(argv$fam)) stop("Missing genotype .fam file: ", argv$fam)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

cell_manifest <- read.table(manifest_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
if (argv$celltype != "all") {
  cell_manifest <- cell_manifest[cell_manifest$celltype == argv$celltype, , drop = FALSE]
  if (nrow(cell_manifest) == 0) stop("Requested cell type is not in manifest: ", argv$celltype)
}

fam_iid <- read.table(argv$fam, stringsAsFactors = FALSE)[[2]]

for (i in seq_len(nrow(cell_manifest))) {
  ct_name <- cell_manifest$celltype[i]
  counts_csv <- file.path(argv$counts_dir, cell_manifest$counts_file[i])
  out_rds <- file.path(out_dir, paste0(ct_name, ".sct.rds"))
  summary_tsv <- file.path(summary_dir, paste0(ct_name, ".summary.tsv"))

  if (!file.exists(counts_csv)) stop("Missing counts file for ", ct_name, ": ", counts_csv)

  cat(sprintf("[%s] reading counts: %s\n", ct_name, counts_csv))
  counts <- read.csv(counts_csv, row.names = 1, check.names = FALSE)
  keep <- intersect(colnames(counts), fam_iid)
  if (length(keep) == 0) stop("No genotyped donor overlap for ", ct_name)

  mat <- as.matrix(counts[, keep, drop = FALSE])
  mat <- round(mat)
  storage.mode(mat) <- "integer"

  cat(sprintf("[%s] donors: counts=%d genotype=%d overlap=%d genes=%d\n",
              ct_name, ncol(counts), length(fam_iid), length(keep), nrow(mat)))
  obj <- CreateSeuratObject(counts = as(mat, "dgCMatrix"))
  obj <- SCTransform(obj, return.only.var.genes = FALSE, verbose = FALSE)
  saveRDS(obj, out_rds)

  feature_count <- nrow(Seurat::SCTResults(obj, slot = "feature.attributes"))
  summary <- data.frame(
    celltype = ct_name,
    counts_genes = nrow(counts),
    counts_donors = ncol(counts),
    genotype_donors = length(fam_iid),
    overlap_donors = length(keep),
    sct_features = feature_count,
    output = out_rds,
    stringsAsFactors = FALSE
  )
  write.table(summary, summary_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("[%s] saved: %s\n", ct_name, out_rds))
}
