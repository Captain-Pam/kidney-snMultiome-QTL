#!/usr/bin/env Rscript
# Purpose: Run scTWAS Stage 2 for one GWAS trait and one cell type across autosomes.
# Inputs:
#   - 2_stage2_focus/2_data/sumstats_by_trait/<trait>/<trait>.chr<chr>.sumstats.tsv.
#   - 1_stage1_weights/3_outputs/weights/<celltype>/weights_index.txt and .wgt.RDat files.
#   - 2_stage2_focus/2_data/ldref/chr<chr>.bed/.bim/.fam.
# Outputs:
#   - 2_stage2_focus/3_outputs/<trait>/stage2/<celltype>/<celltype>.chr<chr>.twas.tsv.
#   - 2_stage2_focus/3_outputs/<trait>/stage2/<celltype>/<celltype>.twas.tsv.

suppressMessages({
  library(argparser)
  library(data.table)
  library(parallel)
  library(scTWAS)
})

p <- arg_parser("Run scTWAS Stage 2 for one trait and one cell type.")
p <- add_argument(p, "--trait", help = "GWAS trait name from gwas_manifest.tsv.")
p <- add_argument(p, "--celltype", help = "Cell type name from celltypes.tsv.")
p <- add_argument(p, "--cores", default = 12, help = "Number of chromosomes to process in parallel.")
p <- add_argument(p, "--workflow-dir", default = "",
                  help = "Path to the public Main/sctwas directory.")
argv <- parse_args(p)
argv$cores <- as.integer(argv$cores)

if (argv$workflow_dir == "") stop("--workflow-dir is required.")

stage1_dir <- file.path(argv$workflow_dir, "1_stage1_weights")
stage2_dir <- file.path(argv$workflow_dir, "2_stage2_focus")
sumstat_dir <- file.path(stage2_dir, "2_data/sumstats_by_trait", argv$trait)
weights_dir <- file.path(stage1_dir, "3_outputs/weights", argv$celltype)
weights_index <- file.path(weights_dir, "weights_index.txt")
ld_dir <- file.path(stage2_dir, "2_data/ldref")
out_dir <- file.path(stage2_dir, "3_outputs", argv$trait, "stage2", argv$celltype)

if (!file.exists(weights_index)) stop("Missing weights index: ", weights_index)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

run_chr <- function(cc) {
  sumstat_path <- file.path(sumstat_dir, sprintf("%s.chr%d.sumstats.tsv", argv$trait, cc))
  ld_prefix <- file.path(ld_dir, paste0("chr", cc))
  out_path <- file.path(out_dir, sprintf("%s.chr%d.twas.tsv", argv$celltype, cc))

  if (!file.exists(sumstat_path)) return(NULL)
  if (!file.exists(paste0(ld_prefix, ".bed"))) return(NULL)
  preview <- fread(sumstat_path, nrows = 1, showProgress = FALSE)
  if (nrow(preview) == 0) return(NULL)

  res <- tryCatch(
    suppressMessages(suppressWarnings(
      twas_stage2_sumstats(
        sumstats = sumstat_path,
        weights = weights_index,
        weights_dir = weights_dir,
        ref_ld_chr = ld_prefix,
        chr = cc,
        out = NULL,
        verbose = 0
      )
    )),
    error = function(e) {
      data.frame(
        trait = argv$trait,
        celltype = argv$celltype,
        chr = cc,
        error = conditionMessage(e),
        stringsAsFactors = FALSE
      )
    }
  )
  if (is.null(res) || nrow(res) == 0) return(NULL)
  res$trait <- argv$trait
  res$celltype <- argv$celltype
  res$chr <- cc
  fwrite(as.data.table(res), out_path, sep = "\t")
  res
}

cat(sprintf("[%s %s] Stage 2 start, cores=%d\n", argv$trait, argv$celltype, argv$cores))
res_list <- mclapply(1:22, run_chr, mc.cores = min(argv$cores, 22), mc.preschedule = FALSE)
res_list <- Filter(Negate(is.null), res_list)
if (length(res_list) == 0) stop("No Stage 2 results were produced for ", argv$trait, " / ", argv$celltype)

combined <- rbindlist(lapply(res_list, as.data.table), fill = TRUE)
if ("error" %in% names(combined)) {
  error_rows <- combined[!is.na(error) & error != ""]
  if (nrow(error_rows) > 0) {
    error_path <- file.path(out_dir, paste0(argv$celltype, ".stage2_errors.tsv"))
    fwrite(error_rows, error_path, sep = "\t")
    stop("Stage 2 failed for ", nrow(error_rows), " chromosome(s). See: ", error_path)
  }
  combined[, error := NULL]
}
combined_path <- file.path(out_dir, paste0(argv$celltype, ".twas.tsv"))
fwrite(combined, combined_path, sep = "\t")

n_non_na <- if ("TWAS.P" %in% names(combined)) sum(!is.na(combined$TWAS.P)) else NA_integer_
cat(sprintf("[%s %s] rows=%d nonNA_TWAS.P=%s -> %s\n",
            argv$trait, argv$celltype, nrow(combined), n_non_na, combined_path))
