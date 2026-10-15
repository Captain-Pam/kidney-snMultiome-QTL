#!/usr/bin/env Rscript
# Purpose: Combine per-chromosome Stage 1 summaries and build weights_index.txt files.
# Inputs:
#   - 1_stage1_weights/3_outputs/summaries/stage1/<celltype>.chr<chr>.stage1.tsv.
#   - 1_stage1_weights/3_outputs/weights/<celltype>/<gene>.wgt.RDat.
# Outputs:
#   - 1_stage1_weights/3_outputs/summaries/stage1/<celltype>.stage1_summary.tsv.
#   - 1_stage1_weights/3_outputs/weights/<celltype>/weights_index.txt for Stage 2.
#   - 1_stage1_weights/3_outputs/summaries/stage1/all_celltypes_stage1_summary.tsv.

suppressMessages(library(argparser))

p <- arg_parser("Combine scTWAS Stage 1 summaries and build weights indices.")
p <- add_argument(p, "--celltype", default = "all",
                  help = "Cell type name from celltypes.tsv, or 'all'.")
p <- add_argument(p, "--workflow-dir", default = "",
                  help = "Path to the public Main/sctwas directory.")
argv <- parse_args(p)

if (argv$workflow_dir == "") stop("--workflow-dir is required.")

stage_dir <- file.path(argv$workflow_dir, "1_stage1_weights")
manifest_path <- file.path(stage_dir, "2_data/celltypes.tsv")
summary_dir <- file.path(stage_dir, "3_outputs/summaries/stage1")
weight_root <- file.path(stage_dir, "3_outputs/weights")

if (!file.exists(manifest_path)) stop("Missing cell type manifest: ", manifest_path)
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

cell_manifest <- read.table(manifest_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
celltypes <- cell_manifest$celltype
if (argv$celltype != "all") celltypes <- argv$celltype

overview <- list()
for (ct_name in celltypes) {
  files <- file.path(summary_dir, paste0(ct_name, ".chr", 1:22, ".stage1.tsv"))
  files <- files[file.exists(files)]
  if (length(files) == 0) {
    warning("No Stage 1 summaries found for ", ct_name)
    next
  }

  res <- do.call(rbind, lapply(files, function(f) {
    read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  }))
  res <- res[order(res$chr, res$gene_start, res$gene), ]
  summary_path <- file.path(summary_dir, paste0(ct_name, ".stage1_summary.tsv"))
  write.table(res, summary_path, sep = "\t", quote = FALSE, row.names = FALSE)

  keep <- res[res$status == "trained" & !is.na(res$cv_R2_max) & res$cv_R2_max > 0, ]
  weight_dir <- file.path(weight_root, ct_name)
  idx <- data.frame(
    WGT = paste0(keep$gene, ".wgt.RDat"),
    ID = keep$gene,
    CHR = keep$chr,
    P0 = keep$gene_start,
    P1 = keep$gene_end,
    PANEL = ct_name,
    stringsAsFactors = FALSE
  )
  idx <- idx[file.exists(file.path(weight_dir, idx$WGT)), ]
  write.table(idx, file.path(weight_dir, "weights_index.txt"),
              quote = FALSE, row.names = FALSE)

  overview[[ct_name]] <- data.frame(
    celltype = ct_name,
    summary_files = length(files),
    genes_attempted = nrow(res),
    genes_trained = sum(res$status == "trained"),
    genes_cv_R2_gt0 = nrow(idx),
    genes_cv_R2_gt0_01 = sum(res$status == "trained" & res$cv_R2_max > 0.01, na.rm = TRUE),
    genes_cv_R2_gt0_05 = sum(res$status == "trained" & res$cv_R2_max > 0.05, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  cat(sprintf("[%s] trained=%d kept_cv_R2_gt0=%d -> %s\n",
              ct_name, sum(res$status == "trained"), nrow(idx), file.path(weight_dir, "weights_index.txt")))
}

if (length(overview) > 0) {
  overview_df <- do.call(rbind, overview)
  write.table(overview_df, file.path(summary_dir, "all_celltypes_stage1_summary.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
}
