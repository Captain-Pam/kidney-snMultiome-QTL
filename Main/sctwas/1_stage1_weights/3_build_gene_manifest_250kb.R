#!/usr/bin/env Rscript
# Purpose: Build autosomal gene manifests and cis-range files for 250 kb scTWAS Stage 1.
# Inputs:
#   - 1_stage1_weights/2_data/sct_objects/<celltype>.sct.rds.
#   - GENCODE gene annotation with TSS and gene symbols.
# Outputs:
#   - 1_stage1_weights/2_data/gene_manifests/<celltype>.genes_autosome.tsv.
#   - 1_stage1_weights/2_data/gene_manifests/<celltype>.cis_ranges.txt for PLINK --extract range.

suppressMessages({
  library(argparser)
  library(Seurat)
})

p <- arg_parser("Build per-cell-type 250 kb gene manifests.")
p <- add_argument(p, "--celltype", default = "all",
                  help = "Cell type name from celltypes.tsv, or 'all'.")
p <- add_argument(p, "--window-bp", default = 250000,
                  help = "Cis half-window around TSS in base pairs.")
p <- add_argument(p, "--workflow-dir", default = "",
                  help = "Path to the public Main/sctwas directory.")
p <- add_argument(p, "--gencode", default = "",
                  help = "GENCODE gene info file with TSS.")
argv <- parse_args(p)
argv$window_bp <- as.integer(argv$window_bp)

if (argv$workflow_dir == "") stop("--workflow-dir is required.")
if (argv$gencode == "") stop("--gencode is required.")

stage_dir <- file.path(argv$workflow_dir, "1_stage1_weights")
manifest_path <- file.path(stage_dir, "2_data/celltypes.tsv")
sct_dir <- file.path(stage_dir, "2_data/sct_objects")
out_dir <- file.path(stage_dir, "2_data/gene_manifests")
summary_dir <- file.path(stage_dir, "3_outputs/summaries/manifests")

if (!file.exists(manifest_path)) stop("Missing cell type manifest: ", manifest_path)
if (!file.exists(argv$gencode)) stop("Missing GENCODE file: ", argv$gencode)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

cell_manifest <- read.table(manifest_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
if (argv$celltype != "all") {
  cell_manifest <- cell_manifest[cell_manifest$celltype == argv$celltype, , drop = FALSE]
  if (nrow(cell_manifest) == 0) stop("Requested cell type is not in manifest: ", argv$celltype)
}

gencode <- read.table(argv$gencode, header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)
gencode <- gencode[gencode$seqnames %in% paste0("chr", 1:22) &
                     gencode$gene_type %in% c("protein_coding", "lncRNA"), ]
gencode <- gencode[!duplicated(gencode$gene_symbol), ]

for (i in seq_len(nrow(cell_manifest))) {
  ct_name <- cell_manifest$celltype[i]
  sct_path <- file.path(sct_dir, paste0(ct_name, ".sct.rds"))
  out_path <- file.path(out_dir, paste0(ct_name, ".genes_autosome.tsv"))
  range_path <- file.path(out_dir, paste0(ct_name, ".cis_ranges.txt"))
  summary_path <- file.path(summary_dir, paste0(ct_name, ".summary.tsv"))

  if (!file.exists(sct_path)) stop("Missing SCT object for ", ct_name, ": ", sct_path)
  sct_genes <- rownames(readRDS(sct_path))
  gc <- gencode[gencode$gene_symbol %in% sct_genes, ]

  man <- data.frame(
    gene = gc$gene_symbol,
    chr = as.integer(sub("chr", "", gc$seqnames)),
    tss = gc$TSS,
    gene_start = gc$start,
    gene_end = gc$end,
    cis_start = pmax(1, gc$TSS - argv$window_bp),
    cis_end = gc$TSS + argv$window_bp,
    stringsAsFactors = FALSE
  )
  man <- man[order(man$chr, man$cis_start), ]
  write.table(man, out_path, sep = "\t", quote = FALSE, row.names = FALSE)

  range_df <- data.frame(man$chr, man$cis_start, man$cis_end, man$gene)
  write.table(range_df, range_path, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)

  summary <- data.frame(
    celltype = ct_name,
    sct_genes = length(sct_genes),
    manifest_genes = nrow(man),
    window_bp = argv$window_bp,
    output = out_path,
    range_output = range_path,
    stringsAsFactors = FALSE
  )
  write.table(summary, summary_path, sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("[%s] manifest genes=%d window=%d -> %s\n",
              ct_name, nrow(man), argv$window_bp, out_path))
}
