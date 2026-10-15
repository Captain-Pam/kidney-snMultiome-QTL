#!/usr/bin/env Rscript
# Purpose: Train scTWAS Stage 1 weights for one cell type and one chromosome.
# Inputs:
#   - 1_stage1_weights/2_data/sct_objects/<celltype>.sct.rds.
#   - 1_stage1_weights/2_data/qcovar/<celltype>.qcovar.txt.
#   - 1_stage1_weights/2_data/gene_manifests/<celltype>.genes_autosome.tsv.
#   - 1_stage1_weights/2_data/genotype/chr<chr>.bed/.bim/.fam.
# Outputs:
#   - 1_stage1_weights/3_outputs/weights/<celltype>/<gene>.wgt.RDat.
#   - 1_stage1_weights/3_outputs/summaries/stage1/<celltype>.chr<chr>.stage1.tsv.

suppressMessages({
  library(argparser)
  library(scTWAS)
  library(parallel)
})

p <- arg_parser("Run scTWAS Stage 1 for one cell type and chromosome.")
p <- add_argument(p, "--celltype", help = "Cell type name from celltypes.tsv.")
p <- add_argument(p, "--chr", help = "Autosome chromosome number, 1-22.")
p <- add_argument(p, "--cores", default = 24, help = "Number of parallel cores.")
p <- add_argument(p, "--crossval", default = 5, help = "scTWAS cross-validation folds.")
p <- add_argument(p, "--niter", default = 5, help = "scTWAS SCT iteration count.")
p <- add_argument(p, "--alpha", default = 0.5, help = "Elastic-net alpha.")
p <- add_argument(p, "--workflow-dir", default = "",
                  help = "Path to the public Main/sctwas directory.")
p <- add_argument(p, "--plink",
                  default = Sys.getenv("PLINK_BIN", "plink"),
                  help = "PLINK binary.")
argv <- parse_args(p)

argv$chr <- as.integer(argv$chr)
argv$cores <- as.integer(argv$cores)
argv$crossval <- as.integer(argv$crossval)
argv$niter <- as.integer(argv$niter)
argv$alpha <- as.numeric(argv$alpha)

if (is.na(argv$chr) || !(argv$chr %in% 1:22)) stop("--chr must be 1-22")
if (argv$workflow_dir == "") stop("--workflow-dir is required.")

stage_dir <- file.path(argv$workflow_dir, "1_stage1_weights")
sct_path <- file.path(stage_dir, "2_data/sct_objects", paste0(argv$celltype, ".sct.rds"))
qcov_path <- file.path(stage_dir, "2_data/qcovar", paste0(argv$celltype, ".qcovar.txt"))
manifest_path <- file.path(stage_dir, "2_data/gene_manifests", paste0(argv$celltype, ".genes_autosome.tsv"))
geno_prefix <- file.path(stage_dir, "2_data/genotype", paste0("chr", argv$chr))
weight_dir <- file.path(stage_dir, "3_outputs/weights", argv$celltype)
tmp_dir <- file.path(stage_dir, "3_outputs/tmp", argv$celltype, paste0("chr", argv$chr))
summary_path <- file.path(stage_dir, "3_outputs/summaries/stage1",
                          paste0(argv$celltype, ".chr", argv$chr, ".stage1.tsv"))

for (path in c(sct_path, qcov_path, manifest_path, paste0(geno_prefix, ".bed"), paste0(geno_prefix, ".bim"), paste0(geno_prefix, ".fam"))) {
  if (!file.exists(path)) stop("Missing required input: ", path)
}
dir.create(weight_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

sct_obj <- readRDS(sct_path)
manifest <- read.table(manifest_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
sub_manifest <- manifest[manifest$chr == argv$chr, , drop = FALSE]
if (nrow(sub_manifest) == 0) stop("No genes in manifest for ", argv$celltype, " chr", argv$chr)

train_one <- function(i) {
  gene <- sub_manifest$gene[i]
  safe_gene <- gsub("[^A-Za-z0-9_.-]", "_", gene)
  cis_prefix <- file.path(tmp_dir, paste0("cis_", safe_gene, "_", Sys.getpid()))
  train_tmp <- file.path(tmp_dir, paste0(".train_", safe_gene, "_", Sys.getpid(), "_"))

  plink_status <- suppressWarnings(system2(
    argv$plink,
    c("--bfile", geno_prefix,
      "--chr", argv$chr,
      "--from-bp", sub_manifest$cis_start[i],
      "--to-bp", sub_manifest$cis_end[i],
      "--maf", "0.05",
      "--make-bed",
      "--out", cis_prefix),
    stdout = FALSE,
    stderr = FALSE
  ))
  on.exit(unlink(Sys.glob(paste0(cis_prefix, ".*"))), add = TRUE)

  if (plink_status != 0 || !file.exists(paste0(cis_prefix, ".bed"))) {
    return(data.frame(
      celltype = argv$celltype, gene = gene, chr = argv$chr,
      status = "plink_failed", cv_R2_max = NA_real_, N = NA_integer_,
      nwgt_max = NA_integer_, gene_start = sub_manifest$gene_start[i],
      gene_end = sub_manifest$gene_end[i], error = "PLINK cis subset failed",
      stringsAsFactors = FALSE
    ))
  }
  if (file.exists(paste0(cis_prefix, ".bim")) && length(readLines(paste0(cis_prefix, ".bim"), n = 1)) == 0) {
    return(data.frame(
      celltype = argv$celltype, gene = gene, chr = argv$chr,
      status = "no_cis_snps", cv_R2_max = NA_real_, N = NA_integer_,
      nwgt_max = NA_integer_, gene_start = sub_manifest$gene_start[i],
      gene_end = sub_manifest$gene_end[i], error = "No cis SNPs after filtering",
      stringsAsFactors = FALSE
    ))
  }

  fit <- tryCatch(
    train_sctwas(
      bfile = cis_prefix,
      sctobj = sct_obj,
      gene = gene,
      plink_bin = argv$plink,
      qcovar = qcov_path,
      crossval = argv$crossval,
      niter = argv$niter,
      alpha = argv$alpha,
      out_prefix = file.path(weight_dir, gene),
      tmp_prefix = train_tmp,
      verbose = 0
    ),
    error = function(e) e
  )

  if (inherits(fit, "error")) {
    return(data.frame(
      celltype = argv$celltype, gene = gene, chr = argv$chr,
      status = "train_failed", cv_R2_max = NA_real_, N = NA_integer_,
      nwgt_max = NA_integer_, gene_start = sub_manifest$gene_start[i],
      gene_end = sub_manifest$gene_end[i], error = conditionMessage(fit),
      stringsAsFactors = FALSE
    ))
  }

  cv_values <- suppressWarnings(as.numeric(fit$cv_performance[1, ]))
  best_idx <- if (all(is.na(cv_values))) NA_integer_ else which.max(cv_values)
  r2 <- if (is.na(best_idx)) NA_real_ else cv_values[best_idx]
  if (!is.finite(r2)) r2 <- NA_real_
  nwgt <- if (length(best_idx) == 1 && is.finite(best_idx)) {
    sum(fit$wgt_matrix[, best_idx] != 0, na.rm = TRUE)
  } else {
    NA_integer_
  }

  data.frame(
    celltype = argv$celltype, gene = gene, chr = argv$chr,
    status = "trained", cv_R2_max = r2, N = fit$N, nwgt_max = nwgt,
    gene_start = sub_manifest$gene_start[i], gene_end = sub_manifest$gene_end[i],
    error = "",
    stringsAsFactors = FALSE
  )
}

cat(sprintf("[%s chr%d] genes=%d cores=%d crossval=%d niter=%d alpha=%.2f\n",
            argv$celltype, argv$chr, nrow(sub_manifest), argv$cores,
            argv$crossval, argv$niter, argv$alpha))
res_list <- mclapply(seq_len(nrow(sub_manifest)), train_one,
                     mc.cores = argv$cores, mc.preschedule = TRUE)
res <- do.call(rbind, res_list)
write.table(res, summary_path, sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf("[%s chr%d] trained=%d failed=%d CV_R2>0=%d summary=%s\n",
            argv$celltype, argv$chr,
            sum(res$status == "trained"),
            sum(res$status != "trained"),
            sum(res$status == "trained" & res$cv_R2_max > 0, na.rm = TRUE),
            summary_path))
