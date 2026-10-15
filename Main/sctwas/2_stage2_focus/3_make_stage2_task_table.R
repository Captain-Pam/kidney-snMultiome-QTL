#!/usr/bin/env Rscript
# Purpose: Build a trait x cell type task table for scTWAS Stage 2 and FOCUS jobs.
# Inputs:
#   - 2_stage2_focus/2_data/gwas_manifest.tsv.
#   - 1_stage1_weights/2_data/celltypes.tsv.
# Outputs:
#   - 2_stage2_focus/2_data/stage2_tasks.tsv.

suppressMessages({
  library(argparser)
  library(data.table)
})

p <- arg_parser("Build scTWAS Stage 2 and FOCUS task table.")
p <- add_argument(p, "--trait", default = "all",
                  help = "Trait name, comma-separated traits, or 'all'.")
p <- add_argument(p, "--workflow-dir", default = "",
                  help = "Path to the public Main/sctwas directory.")
argv <- parse_args(p)

if (argv$workflow_dir == "") stop("--workflow-dir is required.")

stage1_dir <- file.path(argv$workflow_dir, "1_stage1_weights")
stage2_dir <- file.path(argv$workflow_dir, "2_stage2_focus")
gwas_manifest_path <- file.path(stage2_dir, "2_data/gwas_manifest.tsv")
celltype_path <- file.path(stage1_dir, "2_data/celltypes.tsv")
task_path <- file.path(stage2_dir, "2_data/stage2_tasks.tsv")

if (!file.exists(gwas_manifest_path)) stop("Missing GWAS manifest: ", gwas_manifest_path)
if (!file.exists(celltype_path)) stop("Missing cell type manifest: ", celltype_path)

gwas_manifest <- fread(gwas_manifest_path, sep = "\t", header = TRUE, data.table = FALSE)
celltypes <- fread(celltype_path, sep = "\t", header = TRUE, data.table = FALSE)

traits <- if (argv$trait == "all") {
  gwas_manifest$trait
} else {
  unlist(strsplit(argv$trait, ",", fixed = TRUE))
}
missing_traits <- setdiff(traits, gwas_manifest$trait)
if (length(missing_traits) > 0) stop("Traits not found in manifest: ", paste(missing_traits, collapse = ", "))

tasks <- expand.grid(
  trait = traits,
  celltype = celltypes$celltype,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
tasks$task_id <- seq_len(nrow(tasks))
tasks$stage2_sumstats_ready <- file.exists(file.path(
  stage2_dir, "2_data/sumstats_by_trait", tasks$trait,
  sprintf("%s.chr1.sumstats.tsv", tasks$trait)
))
tasks$weights_index_ready <- file.exists(file.path(
  stage1_dir, "3_outputs/weights", tasks$celltype, "weights_index.txt"
))
tasks <- tasks[, c("task_id", "trait", "celltype", "stage2_sumstats_ready", "weights_index_ready")]

write.table(tasks, task_path, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("Task table rows=%d -> %s\n", nrow(tasks), task_path))
cat(sprintf("Submit Stage 2 with: sbatch --array=1-%d 5_submit_stage2_array.sh\n", nrow(tasks)))
