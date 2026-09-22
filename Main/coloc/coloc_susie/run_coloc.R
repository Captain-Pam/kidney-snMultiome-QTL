library(coloc)
library(data.table)
library(rtracklayer)

var2pos <- function(var_id) {
    pos <- as.numeric(sapply(strsplit(var_id, ':'), function(x) {x[2]}))
    return(pos)
}

########
# eGFR #
########
egfr_path <- '../../susie/gwas_finemapping/outputs/'
egfr_res <- read.csv('../../susie/gwas_finemapping/susie_meta.csv')
egfr_res <- egfr_res[egfr_res$cs>0,] # credible sets
egfr_loci <- egfr_res$X
egfr_gr <- readRDS('../../susie/gwas_finemapping/egfr_gr.rds')
egfr_gr <- egfr_gr[egfr_loci]

genes_gr <- readRDS('../../susie/eqtl_finemapping/gene_gr_hg19.rds')

dir.create('outputs')
cts <- gsub('.parquet', '', list.files('data/eqtl/tensorqtl/all_parquets'))

output_list <- list()

for (i in cts) {
    eqtl_res <- read.csv(sprintf('../../susie/eqtl_finemapping/susie_%s.csv', i), row.names=1)
    eqtl_res <- eqtl_res[which(eqtl_res$cs>0), ]
    eqtl_genes <- rownames(eqtl_res)

    eqtl_gene_gr <- genes_gr[[i]]
    names(eqtl_gene_gr) <- eqtl_gene_gr$gene
    eqtl_gene_gr <- eqtl_gene_gr[eqtl_genes]
    
    gene_path <- sprintf('../../susie/eqtl_finemapping/%s', i)
    idx <- queryHits(findOverlaps(egfr_gr, eqtl_gene_gr))
    tmp_loci <- names(egfr_gr)[idx]
    tmp_loci <- unique(tmp_loci)

    res_list <- list()
    for (l in tmp_loci) {
        
        ovlp <- subjectHits(findOverlaps(egfr_gr[l], eqtl_gene_gr))        
        tmp_genes <- eqtl_genes[unique(ovlp)]
        
        rss1 <- readRDS(sprintf('%s/%s/rss.rds', egfr_path, l))
        t1 <- fread(sprintf('%s/%s/res_susie.csv', egfr_path, l))
        setkey(t1, id_harmonized)
    
        for (g in tmp_genes) {
            rss2 <- readRDS(sprintf('%s/%s/rss.rds', gene_path, g))
            t2 <- fread(sprintf('%s/%s/res_susie.csv', gene_path, g))
    
            setkey(t2, id_harmonized)
            ovlp <- intersect(t1$id_harmonized, t2$id_harmonized)

            if (length(ovlp)<50) {next}
            
            toplot <- data.frame(vid=ovlp, 
                                 egfr=t1[ovlp,c('variable_prob','cs')], 
                                 eqtl=t2[ovlp,c('variable_prob','cs')], 
                                 pos=var2pos(ovlp))
    
            susie.res <- coloc.susie(rss1,rss2)
            res <- susie.res$summary
            res$loci <- l
            res$gene <- g
            
            res_list[[length(res_list) + 1 ]] <- res
    
            if (all(is.na(res$PP.H4.abf))) {next}
            if (max(res$PP.H4.abf) < 0.7) {next}
    
            outdir <- sprintf('outputs/%s_%s_%s', i, l, g)
            dir.create(outdir)
            saveRDS(susie.res, sprintf('%s/coloc_res.rds', outdir))
            fwrite(res, sprintf('%s/coloc_res.csv', outdir))
    
        }
        
        print(l)
    }
    res <- do.call(rbind, res_list)
    res <- res[!is.na(res$PP.H4.abf),]
    res$celltype <- i
    output_list[[i]] <- res
}
outputs <- do.call(rbind, output_list)
fwrite(outputs, 'coloc_all.csv')
