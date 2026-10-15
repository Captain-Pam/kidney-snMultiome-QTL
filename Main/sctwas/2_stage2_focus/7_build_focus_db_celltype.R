#!/usr/bin/env Rscript
# Purpose: Build a FOCUS SQLite weight database from scTWAS Stage 1 .wgt.RDat files for one cell type.
# Inputs:
#   - 1_stage1_weights/3_outputs/weights/<celltype>/weights_index.txt.
#   - 1_stage1_weights/2_data/gene_manifests/<celltype>.genes_autosome.tsv.
# Outputs:
#   - 2_stage2_focus/3_outputs/focus_db/<celltype>.db.

suppressMessages({
  library(argparser)
  library(data.table)
  library(DBI)
  library(RSQLite)
})

p <- arg_parser("Build FOCUS DB from scTWAS weights.")
p <- add_argument(p, "--celltype", help = "Cell type name from celltypes.tsv.")
p <- add_argument(p, "--workflow-dir", default = "",
                  help = "Path to the public Main/sctwas directory.")
argv <- parse_args(p)

if (argv$workflow_dir == "") stop("--workflow-dir is required.")

stage1_dir <- file.path(argv$workflow_dir, "1_stage1_weights")
stage2_dir <- file.path(argv$workflow_dir, "2_stage2_focus")
weights_dir <- file.path(stage1_dir, "3_outputs/weights", argv$celltype)
weights_index <- file.path(weights_dir, "weights_index.txt")
gene_manifest <- file.path(stage1_dir, "2_data/gene_manifests", paste0(argv$celltype, ".genes_autosome.tsv"))
out_dir <- file.path(stage2_dir, "3_outputs/focus_db")
db_path <- file.path(out_dir, paste0(argv$celltype, ".db"))

if (!file.exists(weights_index)) stop("Missing weights index: ", weights_index)
if (!file.exists(gene_manifest)) stop("Missing gene manifest: ", gene_manifest)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

idx <- fread(weights_index, data.table = FALSE)
man <- fread(gene_manifest, data.table = FALSE)
idx <- idx[file.exists(file.path(weights_dir, idx$WGT)), ]
idx <- idx[idx$ID %in% man$gene, ]
if (nrow(idx) == 0) stop("No valid weights in index for ", argv$celltype)
man_match <- man[match(idx$ID, man$gene), ]

if (file.exists(db_path)) file.remove(db_path)
con <- dbConnect(SQLite(), db_path)
on.exit(dbDisconnect(con), add = TRUE)

dbExecute(con, "CREATE TABLE refpanel (ref_name VARCHAR(128) NOT NULL, tissue VARCHAR(128) NOT NULL, assay VARCHAR(128), id INTEGER NOT NULL, PRIMARY KEY (id))")
dbExecute(con, "CREATE TABLE molecularfeature (ens_gene_id VARCHAR(64) NOT NULL, ens_tx_id VARCHAR(64), mol_name VARCHAR(64), type VARCHAR(64), chrom VARCHAR(10) NOT NULL, tx_start INTEGER, tx_stop INTEGER, id INTEGER NOT NULL, PRIMARY KEY (id))")
dbExecute(con, "CREATE TABLE model (inference VARCHAR(128) NOT NULL, ref_id INTEGER, mol_id INTEGER, id INTEGER NOT NULL, PRIMARY KEY (id), FOREIGN KEY(ref_id) REFERENCES refpanel(id), FOREIGN KEY(mol_id) REFERENCES molecularfeature(id))")
dbExecute(con, "CREATE TABLE modelattribute (attr_name VARCHAR(128) NOT NULL, value FLOAT, model_id INTEGER, id INTEGER NOT NULL, PRIMARY KEY (id), FOREIGN KEY(model_id) REFERENCES model(id))")
dbExecute(con, "CREATE TABLE weight (snp VARCHAR(128) NOT NULL, chrom VARCHAR(2) NOT NULL, pos INTEGER NOT NULL, effect_allele VARCHAR(32) NOT NULL, alt_allele VARCHAR(32) NOT NULL, effect FLOAT NOT NULL, std_error FLOAT, model_id INTEGER, id INTEGER NOT NULL, PRIMARY KEY (id))")
dbWriteTable(con, "refpanel",
             data.frame(ref_name = paste0("kidney_", argv$celltype, "_250kb"),
                        tissue = argv$celltype, assay = "rnaseq", id = 1L),
             append = TRUE)

molecular_features <- vector("list", nrow(idx))
models <- vector("list", nrow(idx))
attributes <- vector("list", nrow(idx))
weights <- vector("list", nrow(idx))
weight_id <- 0L
attribute_id <- 0L
gene_id <- 0L

for (i in seq_len(nrow(idx))) {
  wgt_path <- file.path(weights_dir, idx$WGT[i])
  e <- new.env()
  ok <- tryCatch({ load(wgt_path, envir = e); TRUE }, error = function(z) FALSE)
  if (!ok) next
  if (is.null(e$wgt.matrix) || is.null(e$cv.performance) || is.null(e$snps)) next
  if (is.null(dim(e$wgt.matrix))) next

  keep_model <- apply(e$wgt.matrix, 2, function(col) !any(is.na(col)) && any(col != 0))
  if (!any(keep_model)) next
  wgt_matrix <- e$wgt.matrix[, keep_model, drop = FALSE]
  cv_perf <- e$cv.performance[, keep_model, drop = FALSE]
  rsq <- suppressWarnings(as.numeric(cv_perf["rsq", ]))
  pval <- suppressWarnings(as.numeric(cv_perf["pval", ]))
  if (all(is.na(rsq))) next
  best_idx <- which.max(rsq)
  model_weights <- wgt_matrix[, best_idx]
  nonzero <- which(model_weights != 0 & !is.na(model_weights))
  if (length(nonzero) == 0) next

  gene_id <- gene_id + 1L
  gene_name <- idx$ID[i]
  n_nonzero <- length(nonzero)
  snps <- e$snps

  molecular_features[[gene_id]] <- data.frame(
    ens_gene_id = gene_name,
    ens_tx_id = NA_character_,
    mol_name = gene_name,
    type = NA_character_,
    chrom = as.character(man_match$chr[i]),
    tx_start = man_match$gene_start[i],
    tx_stop = man_match$gene_end[i],
    id = gene_id,
    stringsAsFactors = FALSE
  )
  models[[gene_id]] <- data.frame(
    inference = "ElasticNet",
    ref_id = 1L,
    mol_id = gene_id,
    id = gene_id,
    stringsAsFactors = FALSE
  )
  attributes[[gene_id]] <- data.frame(
    attr_name = c("cv.R2", "cv.R2.pval"),
    value = c(rsq[best_idx], pval[best_idx]),
    model_id = gene_id,
    id = c(attribute_id + 1L, attribute_id + 2L),
    stringsAsFactors = FALSE
  )
  attribute_id <- attribute_id + 2L
  weights[[gene_id]] <- data.frame(
    snp = as.character(snps$V2[nonzero]),
    chrom = as.character(snps$V1[nonzero]),
    pos = as.integer(snps$V4[nonzero]),
    effect_allele = as.character(snps$V5[nonzero]),
    alt_allele = as.character(snps$V6[nonzero]),
    effect = as.numeric(model_weights[nonzero]),
    std_error = NA_real_,
    model_id = gene_id,
    id = (weight_id + 1L):(weight_id + n_nonzero),
    stringsAsFactors = FALSE
  )
  weight_id <- weight_id + n_nonzero
}

if (gene_id == 0L) stop("No valid models were written for ", argv$celltype)
dbWriteTable(con, "molecularfeature", rbindlist(molecular_features[seq_len(gene_id)]), append = TRUE)
dbWriteTable(con, "model", rbindlist(models[seq_len(gene_id)]), append = TRUE)
dbWriteTable(con, "modelattribute", rbindlist(attributes[seq_len(gene_id)]), append = TRUE)
dbWriteTable(con, "weight", rbindlist(weights[seq_len(gene_id)]), append = TRUE)

cat(sprintf("[%s] FOCUS DB genes=%d weight_rows=%d -> %s\n",
            argv$celltype, gene_id, weight_id, db_path))
