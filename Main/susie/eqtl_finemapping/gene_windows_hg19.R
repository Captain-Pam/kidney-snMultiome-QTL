library(rtracklayer)
library(data.table)

var2pos <- function(var_id) {
    pos <- as.numeric(sapply(strsplit(var_id, ':'), function(x) {x[2]}))
    return(pos)
}

cts <- gsub('.parquet', '', list.files('data/eqtl/tensorqtl/all_parquets'))

res_list <- list()
for (i in cts) {
    res_path <- list.files(path = i, pattern = "res_susie.csv", recursive = TRUE, full.names = TRUE)
    res_gr <- do.call(c, lapply(res_path, function(x) {
        g <- gsub('.*/', '', gsub('/res_susie.csv', '', x))
        a <- fread(x)
        pos <- var2pos(a$id_harmonized)
        left <- min(pos)
        right <- max(pos)
        chr <- gsub(':.*', '', a$id_harmonized[1])
        gr <- GRanges(chr, IRanges(left, right))
        gr$gene <- g        
        return(gr)
    }))
    res_list[[i]] <- res_gr
    print(i)
}
saveRDS(res_list, 'gene_gr_hg19.rds')
