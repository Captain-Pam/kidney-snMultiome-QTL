library(coloc)
library(susieR)
library(data.table)
library(arrow)
library(ggplot2)
library(gridExtra)
library(GenomicRanges)

read_ld <- function(ld_file){
    a <- data.frame(fread(ld_file))
    ids <- a[,1]
    a <- a[,-1]
    rownames(a) <- colnames(a) <- ids
    return(a)
}

##################
# input argument #
##################
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("No input provided. Usage: Rscript script.R <x>")
}
locus <- args[1]

gr_indep <- readRDS('egfr_gr.rds')
egfr_list <- readRDS('egfr_list_1m.rds')

susie_meta <- rep(NA, 7)
names(susie_meta) <- c('total_snps','outlier_snps','finemap_snps','lambda','alignment','converge','cs')
outlier_thres <- 5
N <- 1.7e6
LD_gr <- readRDS('../ld_reference/LD_gr.rds') # from ld_reference/process_LD.R

outdir <- paste0('outputs/', locus, '/')
dir.create(outdir)
        
egfr_tmp <- egfr_list[[locus]] # eGFR sum stats
egfr_tmp$hg19_id1 <- egfr_tmp$vid
egfr_tmp$hg19_id2 <- egfr_tmp$vid2

snp_file <- paste0(outdir, 'tmp_snps.csv')
ld_file <-  paste0(outdir, 'ld_m_tmp.csv')

fwrite(data.frame(mcols(egfr_tmp)), snp_file)
ovlp_snps <- table(subjectHits(findOverlaps(egfr_tmp, LD_gr)))
LD_id <- as.numeric(names(ovlp_snps)[which.max(ovlp_snps)])
ld_prefix <- gsub('.gz', '', LD_gr[LD_id]$LD_gz)
cmd <- sprintf('../ld_reference/query_LD.py --ld_prefix %s --input_file %s --output_file %s', 
               ld_prefix, snp_file, ld_file)
system(cmd)

ld_m <- tryCatch({
    read_ld(ld_file)
}, error = function(e) {
    NULL
})

if (is.null(ld_m)) {
    system(paste0('rm -r ', outdir))
    susie_meta['total_snps'] <- 0
    write.csv(susie_meta, paste0(outdir, 'meta.csv')) 
    stop("stopped due to snps not found in ref panel.")
}

colnames(ld_m) <- rownames(ld_m)

# harmonize effect allele with ref panel
egfr_dt <- data.table(data.frame(mcols(egfr_tmp)))
egfr_dt <- egfr_dt[egfr_dt$hg19_id1 %in% rownames(ld_m) | egfr_dt$hg19_id2 %in% rownames(ld_m)]
egfr_dt$id_harmonized <- egfr_dt$hg19_id2
toflip <- egfr_dt$hg19_id1 %in% rownames(ld_m)
if (sum(toflip)>0) {
    egfr_dt$id_harmonized[toflip] <- egfr_dt$hg19_id1[toflip]
    egfr_dt$BETA[toflip] <- (-egfr_dt$BETA[toflip])
}
egfr_dt <- egfr_dt[!duplicated(egfr_dt$id_harmonized)]
setkey(egfr_dt, id_harmonized)        
egfr_dt <- egfr_dt[rownames(ld_m)]
susie_meta['total_snps'] <- nrow(egfr_dt)

# remove outlier snps
# https://github.com/stephenslab/susieR/issues/182
z <- egfr_dt$BETA / egfr_dt$SE
condz_in <- kriging_rss(z, data.matrix(ld_m), n=N)
outlier <- which(abs(condz_in$conditional_dist$z_std_diff)>outlier_thres)
susie_meta['outlier_snps'] <- length(outlier)

if (length(outlier)>0) {
    egfr_dt <- egfr_dt[-outlier]
    ld_m <- ld_m[-outlier, -outlier]            
}
susie_meta['finemap_snps'] <- nrow(egfr_dt)

d2 <- list()
d2[['beta']] <- egfr_dt$BETA
d2[['varbeta']] <- egfr_dt$SE**2
d2[['snp']] <- egfr_dt$id_harmonized
d2[['position']] <- egfr_dt$POS
d2[['type']] <- 'quant'            
d2[['N']] <- N
d2[['MAF']] <- egfr_dt$MAF
d2[['LD']] <- data.matrix(ld_m)

susie_meta['lambda'] <- estimate_s_rss(d2[['beta']]/sqrt(d2[['varbeta']]), d2[['LD']], N)
susie_meta['alignment'] <- check_alignment(d2, do_plot=F)

fitted_rss <- tryCatch({
    runsusie(d2, max_iter=500, repeat_until_convergence=F)
}, error = function(e) {
    NULL
})
     
if (is.null(fitted_rss)) {
    system(paste0('rm -r ', outdir))
    susie_meta['converge'] <- 0
    write.csv(susie_meta, paste0(outdir, 'meta.csv'))
    stop("stopped due to susie not converge.")

}

susie_meta['converge'] <- 1
egfr_dt$variable <- 1:nrow(egfr_dt)
res_susie <- merge(x=egfr_dt, y=summary(fitted_rss)$vars, 
                   by="variable", all.x=TRUE)        
saveRDS(fitted_rss, paste0(outdir, 'rss.rds'))
write.csv(res_susie, paste0(outdir, 'res_susie.csv'))
tmp <- res_susie$cs
n_cs <- length(unique(tmp[tmp!=(-1)]))
susie_meta['cs'] <- n_cs
write.csv(susie_meta, paste0(outdir, 'meta.csv'))

res_susie$P.value <- as.numeric(res_susie$P.value)
if (n_cs>0) {
    # plot
    res_susie$cs <- as.factor(res_susie$cs)
    res_susie$P.value[res_susie$P.value==0] <- 1e-300
    pdf(paste(outdir, 'susie.pdf'), 10, 4)
    p1 <- ggplot(res_susie, aes(x=POS-gr_indep[locus]$pos, y=-log10(P.value))) + geom_point(alpha=0.8) + 
        geom_point(data=subset(res_susie, res_susie$cs!=-1), aes(color=cs), size=4) + theme_classic()

    p2 <- ggplot(res_susie, aes(x=POS-gr_indep[locus]$pos, y=variable_prob)) + geom_point() + 
        geom_point(data=subset(res_susie, res_susie$cs!=-1), aes(color=cs), size=4) + theme_classic()
    grid.arrange(p1, p2, nrow=2)
    dev.off()
}

