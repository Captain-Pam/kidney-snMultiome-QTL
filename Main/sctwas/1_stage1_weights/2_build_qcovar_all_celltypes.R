#!/usr/bin/env Rscript
# Purpose: Build scTWAS quantitative covariate files for selected cell types.
# Inputs:
#   - 1_stage1_weights/2_data/celltypes.tsv with covariate and PEER5 file names.
#   - eQTL covariates with Gender and Age rows.
#   - PEER5 covariates from the matched eQTL workflow.
# Outputs:
#   - 1_stage1_weights/2_data/qcovar/<celltype>.qcovar.txt with FID, IID, Gender, Age, PEER1-PEER5.
#   - 1_stage1_weights/3_outputs/summaries/qcovar/<celltype>.summary.tsv.

suppressMessages(library(argparser))

p <- arg_parser("Build qcovar files for scTWAS Stage 1.")
p <- add_argument(p, "--celltype", default = "all",
                  help = "Cell type name from celltypes.tsv, or 'all'.")
p <- add_argument(p, "--workflow-dir", default = "",
                  help = "Path to the public Main/sctwas directory.")
p <- add_argument(p, "--eqtl-dir", default = "",
                  help = "Root directory of the eQTL workflow with covariate and PEER files.")
argv <- parse_args(p)

if (argv$workflow_dir == "") stop("--workflow-dir is required.")
if (argv$eqtl_dir == "") stop("--eqtl-dir is required.")

stage_dir <- file.path(argv$workflow_dir, "1_stage1_weights")
manifest_path <- file.path(stage_dir, "2_data/celltypes.tsv")
out_dir <- file.path(stage_dir, "2_data/qcovar")
summary_dir <- file.path(stage_dir, "3_outputs/summaries/qcovar")

if (!file.exists(manifest_path)) stop("Missing cell type manifest: ", manifest_path)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

cell_manifest <- read.table(manifest_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
if (argv$celltype != "all") {
  cell_manifest <- cell_manifest[cell_manifest$celltype == argv$celltype, , drop = FALSE]
  if (nrow(cell_manifest) == 0) stop("Requested cell type is not in manifest: ", argv$celltype)
}

for (i in seq_len(nrow(cell_manifest))) {
  ct_name <- cell_manifest$celltype[i]
  cov_path <- file.path(argv$eqtl_dir, "3_covariates", cell_manifest$covariate_file[i])
  peer_path <- file.path(argv$eqtl_dir, "5_peer", cell_manifest$peer5_file[i])
  out_path <- file.path(out_dir, paste0(ct_name, ".qcovar.txt"))
  summary_path <- file.path(summary_dir, paste0(ct_name, ".summary.tsv"))

  if (!file.exists(cov_path)) stop("Missing covariate file for ", ct_name, ": ", cov_path)
  if (!file.exists(peer_path)) stop("Missing PEER5 file for ", ct_name, ": ", peer_path)

  covs <- read.table(cov_path, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
  if (!all(c("Gender", "Age") %in% rownames(covs))) {
    stop("Covariate file lacks Gender/Age rows for ", ct_name)
  }
  ga <- as.data.frame(t(covs[c("Gender", "Age"), , drop = FALSE]))
  ga$IID <- rownames(ga)

  peer <- read.table(peer_path, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
  pe <- as.data.frame(t(peer))
  if (ncol(pe) < 5) stop("PEER file has fewer than 5 factors for ", ct_name)
  pe <- pe[, seq_len(5), drop = FALSE]
  colnames(pe) <- paste0("PEER", seq_len(5))
  pe$IID <- rownames(pe)

  merged <- merge(ga, pe, by = "IID")
  merged$FID <- merged$IID
  out <- merged[, c("FID", "IID", "Gender", "Age", paste0("PEER", seq_len(5)))]
  write.table(out, out_path, quote = FALSE, row.names = FALSE)

  summary <- data.frame(
    celltype = ct_name,
    donors = nrow(out),
    covariates = ncol(out) - 2,
    any_na = anyNA(out),
    output = out_path,
    stringsAsFactors = FALSE
  )
  write.table(summary, summary_path, sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("[%s] qcovar donors=%d covariates=%d any_NA=%s -> %s\n",
              ct_name, nrow(out), ncol(out) - 2, anyNA(out), out_path))
}
