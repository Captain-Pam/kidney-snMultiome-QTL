#!/usr/bin/env Rscript
# Purpose: Summarize scTWAS Stage 2 and FOCUS outputs for trait x cell type combinations.
# Inputs:
#   - 2_stage2_focus/3_outputs/<trait>/stage2/<celltype>/<celltype>.twas.tsv.
#   - 2_stage2_focus/3_outputs/<trait>/focus/<celltype>/<celltype>.focus_all.tsv.
# Outputs:
#   - 2_stage2_focus/3_outputs/<trait>/summary/<celltype>.twas_focus_summary.tsv.
#   - 2_stage2_focus/3_outputs/summaries/all_trait_celltype_twas_focus_summary.tsv.

suppressMessages({
  library(argparser)
  library(data.table)
})

p <- arg_parser("Summarize scTWAS Stage 2 and FOCUS outputs.")
p <- add_argument(p, "--trait", default = "all",
                  help = "Trait name, comma-separated traits, or 'all'.")
p <- add_argument(p, "--celltype", default = "all",
                  help = "Cell type name, comma-separated cell types, or 'all'.")
p <- add_argument(p, "--workflow-dir", default = "",
                  help = "Path to the public Main/sctwas directory.")
argv <- parse_args(p)

if (argv$workflow_dir == "") stop("--workflow-dir is required.")

stage1_dir <- file.path(argv$workflow_dir, "1_stage1_weights")
stage2_dir <- file.path(argv$workflow_dir, "2_stage2_focus")
gwas_manifest <- fread(file.path(stage2_dir, "2_data/gwas_manifest.tsv"), data.table = FALSE)
celltypes_manifest <- fread(file.path(stage1_dir, "2_data/celltypes.tsv"), data.table = FALSE)

traits <- if (argv$trait == "all") gwas_manifest$trait else unlist(strsplit(argv$trait, ",", fixed = TRUE))
celltypes <- if (argv$celltype == "all") celltypes_manifest$celltype else unlist(strsplit(argv$celltype, ",", fixed = TRUE))

pick_col <- function(x, candidates) {
  hit <- candidates[candidates %in% names(x)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

summaries <- list()
for (trait_name in traits) {
  for (ct_name in celltypes) {
    twas_path <- file.path(stage2_dir, "3_outputs", trait_name, "stage2", ct_name, paste0(ct_name, ".twas.tsv"))
    focus_path <- file.path(stage2_dir, "3_outputs", trait_name, "focus", ct_name, paste0(ct_name, ".focus_all.tsv"))
    summary_dir <- file.path(stage2_dir, "3_outputs", trait_name, "summary")
    dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

    twas_rows <- NA_integer_
    twas_non_na <- NA_integer_
    twas_R2_gt0_01 <- NA_integer_
    twas_fdr_005 <- NA_integer_
    twas_p_col <- NA_character_
    twas_r2_col <- NA_character_
    twas_annotated_path <- NA_character_
    if (file.exists(twas_path)) {
      twas <- fread(twas_path)
      twas_rows <- nrow(twas)
      twas_p_col <- pick_col(twas, c("TWAS.P_max", "TWAS.P"))
      twas_z_col <- pick_col(twas, c("TWAS.Z_max", "TWAS.Z"))
      twas_r2_col <- pick_col(twas, c("MODELCV.R2_max", "MODELCV.R2"))
      if (!is.na(twas_p_col) && !is.na(twas_r2_col)) {
        p_used <- suppressWarnings(as.numeric(twas[[twas_p_col]]))
        r2_used <- suppressWarnings(as.numeric(twas[[twas_r2_col]]))
        z_used <- if (!is.na(twas_z_col)) suppressWarnings(as.numeric(twas[[twas_z_col]])) else rep(NA_real_, nrow(twas))
        model_quality_pass <- !is.na(r2_used) & r2_used > 0.01

        twas$TWAS.P_used <- p_used
        twas$TWAS.Z_used <- z_used
        twas$MODELCV.R2_used <- r2_used
        twas$model_R2_gt0_01 <- model_quality_pass
        twas$q.BH <- NA_real_
        ok <- !is.na(p_used) & model_quality_pass
        twas_non_na <- sum(!is.na(p_used))
        twas_R2_gt0_01 <- sum(ok)
        twas$q.BH[ok] <- p.adjust(p_used[ok], method = "BH")
        twas_fdr_005 <- sum(twas$q.BH < 0.05, na.rm = TRUE)
        twas_annotated_path <- file.path(stage2_dir, "3_outputs", trait_name, "stage2", ct_name,
                                         paste0(ct_name, ".twas.qBH.tsv"))
        fwrite(twas, twas_annotated_path, sep = "\t")
      }
    }

    focus_rows <- NA_integer_
    focus_pip_05 <- NA_integer_
    focus_pip_08 <- NA_integer_
    focus_pip_09 <- NA_integer_
    if (file.exists(focus_path)) {
      focus <- fread(focus_path)
      focus <- focus[!is.na(mol_name) & mol_name != "NULL"]
      pip_col <- if ("pips_pop1" %in% names(focus)) "pips_pop1" else if ("pip" %in% names(focus)) "pip" else NA_character_
      focus_rows <- nrow(focus)
      if (!is.na(pip_col)) {
        pip <- suppressWarnings(as.numeric(focus[[pip_col]]))
        focus_pip_05 <- sum(pip > 0.5, na.rm = TRUE)
        focus_pip_08 <- sum(pip > 0.8, na.rm = TRUE)
        focus_pip_09 <- sum(pip > 0.9, na.rm = TRUE)
      }
    }

    summary <- data.frame(
      trait = trait_name,
      celltype = ct_name,
      twas_file = twas_path,
      twas_qBH_file = twas_annotated_path,
      twas_rows = twas_rows,
      twas_non_na = twas_non_na,
      twas_R2_gt0_01 = twas_R2_gt0_01,
      twas_p_col_for_FDR = twas_p_col,
      twas_R2_col_for_filter = twas_r2_col,
      twas_qBH_lt_0_05 = twas_fdr_005,
      focus_file = focus_path,
      focus_rows = focus_rows,
      focus_pip_gt_0_5 = focus_pip_05,
      focus_pip_gt_0_8 = focus_pip_08,
      focus_pip_gt_0_9 = focus_pip_09,
      stringsAsFactors = FALSE
    )
    out_path <- file.path(summary_dir, paste0(ct_name, ".twas_focus_summary.tsv"))
    fwrite(summary, out_path, sep = "\t")
    summaries[[paste(trait_name, ct_name, sep = "__")]] <- summary
    cat(sprintf("[%s %s] TWAS rows=%s R2>0.01=%s FDR<0.05=%s p_col=%s FOCUS rows=%s\n",
                trait_name, ct_name, twas_rows, twas_R2_gt0_01, twas_fdr_005, twas_p_col, focus_rows))
  }
}

if (length(summaries) > 0) {
  all_summary <- rbindlist(summaries, fill = TRUE)
  out_all <- file.path(stage2_dir, "3_outputs/summaries/all_trait_celltype_twas_focus_summary.tsv")
  dir.create(dirname(out_all), recursive = TRUE, showWarnings = FALSE)
  fwrite(all_summary, out_all, sep = "\t")
}
