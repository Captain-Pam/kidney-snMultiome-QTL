library(data.table)
library(rtracklayer)
library(arrow)

#############
# sum stats #
#############
# hg19
sumstats <- 'data/gwas/eGFR_sumstats.harmonized.txt.gz' # harmonized eGFR meta-analysis summary statistics (hg19)
egfr <- fread(sumstats)
egfr <- egfr[!is.na(egfr$POS),]
egfr$vid  <- paste0('chr', egfr$CHR, ':', egfr$POS, ':', egfr$ALT,':', egfr$REF)
egfr$vid2 <- paste0('chr', egfr$CHR, ':', egfr$POS, ':', egfr$REF,':', egfr$ALT)
egfr$chr <- paste0('chr', egfr$CHR)
gr_egfr <- GRanges(seqnames=paste0('chr', egfr$CHR), ranges=IRanges(egfr$POS, egfr$POS))
mcols(gr_egfr) <- egfr

#############
# 878 locus #
#############
indep <- readRDS('data/gwas/eGFR_loci.rds') # independent genome-wide significant lead SNPs
indep <- egfr[egfr$MarkerName %in% indep$rsid]
indep <- indep[!is.na(indep$chr)]
gr_indep <- GRanges(indep$chr, 
                   IRanges(indep$POS,
                           indep$POS))
gr_indep$pos <- start(gr_indep)
gr_indep$rsid <- indep$MarkerName
names(gr_indep) <- gr_indep$rsid
gr_indep_1m <- resize(gr_indep, width=1000000, fix='center')

egfr_list <- list()
for (locus in names(gr_indep_1m)) {
    idx <- queryHits(findOverlaps(gr_egfr, gr_indep_1m[locus]))
    tmp <- gr_egfr[idx]    
    egfr_list[[locus]] <- tmp    
    print(locus)
}
saveRDS(gr_indep, 'loci_gr.rds')
saveRDS(gr_indep_1m, 'egfr_gr.rds')
saveRDS(egfr_list, 'egfr_list_1m.rds')
