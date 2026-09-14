library(rtracklayer)
library(data.table)
library(arrow)

var2pos <- function(var_id) {
    pos <- as.numeric(sapply(strsplit(var_id, ':'), function(x) {x[2]}))
    return(pos)
}

cts <- gsub('.parquet', '', list.files('data/eqtl/tensorqtl/all_parquets/'))

res_list <- list()
for (i in cts) {
    x <- read_parquet(sprintf('data/eqtl/tensorqtl/all_parquets/%s.parquet', i))
    x <- as.data.table(x)
    x <- x[x$egene_qval<0.1]
    x$chr <- gsub(':.*', '', x$variant_id)
    x$pos <- var2pos(x$variant_id)
    x_gene <- x[, .(start = min(pos), end = max(pos), chr = first(chr)), by = phenotype_id]
    gr <- GRanges(x_gene$chr, IRanges(x_gene$start-10, x_gene$end+10))
    gr$gene <- x_gene$phenotype_id
    names(gr) <- gr$gene            
    res_list[[i]] <- gr
    print(i)
}

saveRDS(res_list, 'gene_gr.rds')

