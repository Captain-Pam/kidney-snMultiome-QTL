library(data.table)
library(rtracklayer)

# Step 3 - Write GCT expression matrices and the gene-level annotation.
#
# Builds raw.gct / tpm.gct per cell type from the QC'd count CSVs, and a
# gene-level GTF for GTEx expression normalization. Non-unique gene names are
# dropped, gene_id is set to gene_name, and features are relabeled as
# 'transcript' so the GTEx collapse step treats each gene as one interval.
#
# Inputs:
#   prepare_data/sample_qc/<ct>/{count_raw.csv, count_tpm.csv, sample_lookup.txt}
#   data/reference/gencode.v46.basic.annotation.gtf.gz
# Outputs:
#   prepare_data/sample_qc/<ct>/{raw.gct, tpm.gct}
#   prepare_data/sample_qc/gencode_v46_rev.gtf
#
# Run from the eqtl/ directory.

qc_dir <- 'prepare_data/sample_qc'
cts <- gsub('.csv', '', list.files('prepare_data/anndata/counts'))

write_gct <- function(csv, path) {
    m <- read.csv(csv)
    m <- m[, c('X', 'X', sample_ids)]
    colnames(m)[1:2] <- c('Name', 'Description')
    con <- file(path, 'w')
    writeLines('#1.2', con)
    writeLines(paste(dim(m) - c(0, 2), collapse = ' '), con)
    close(con)
    write.table(m, path, append = TRUE, sep = '\t', quote = FALSE, row.names = FALSE)
}

for (ct in cts) {
    sample_ids <- read.table(sprintf('%s/%s/sample_lookup.txt', qc_dir, ct),
                             header = TRUE)$sample_ids
    write_gct(sprintf('%s/%s/count_raw.csv', qc_dir, ct), sprintf('%s/%s/raw.gct', qc_dir, ct))
    write_gct(sprintf('%s/%s/count_tpm.csv', qc_dir, ct), sprintf('%s/%s/tpm.gct', qc_dir, ct))
    print(ct)
}

# gene-level annotation
gtf <- import('data/reference/gencode.v46.basic.annotation.gtf.gz')
gtf <- gtf[gtf$type == 'gene']
non_unique <- unique(gtf$gene_name[duplicated(gtf$gene_name)])
gtf <- gtf[!gtf$gene_name %in% non_unique]
gtf$gene_id <- gtf$gene_name
gtf$type <- 'transcript'
export(gtf, sprintf('%s/gencode_v46_rev.gtf', qc_dir))
