library(coloc)
library(data.table)
library(arrow)
library(ggplot2)
library(gridExtra)
library(rtracklayer)
library(susieR)


var2pos <- function(var_id) {
    pos <- as.numeric(sapply(strsplit(var_id, ':'), function(x) {x[2]}))
    return(pos)
}

cts <- gsub('.parquet', '', list.files('data/eqtl/tensorqtl/all_parquets'))

####################################
# take cell type as input argument #
####################################
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("No input provided. Usage: Rscript script.R <x>")
}
ct <- args[1]
dir.create(ct)

eqtl <- read_parquet(sprintf('data/eqtl/tensorqtl/all_parquets/%s.parquet', ct))
eqtl <- as.data.table(eqtl)
eqtl <- eqtl[eqtl$egene_qval<0.1]

# other info
query_ld_script <- '../ld_reference/query_LD.py'
chain <- import.chain('data/reference/hg38ToHg19.over.chain')
LD_gr <- readRDS('../ld_reference/LD_gr.rds') # from ld_reference/process_LD.R
outlier_thres <- 5 # susie snp outlier threshold

Y <- data.frame(fread(sprintf('data/eqtl/expression_bed/%s.expression.bed.gz', ct)))
rownames(Y) <- Y$gene_id
Y <- Y[, 5:ncol(Y)]
Y_sd <- apply(Y, 1, sd)
N <- ncol(Y)

susie_meta <- matrix(NA, nrow=length(unique(eqtl$phenotype_id)), ncol=7,
                     dimnames=list(sort(unique(eqtl$phenotype_id)), 
                                   c('total_snps','outlier_snps','finemap_snps','lambda','alignment','converge','cs')))

for (pheno_id in rownames(susie_meta)) {
        
    outdir <- paste0(ct, '/', pheno_id, '/')
    dir.create(outdir)
        
    tmp_eqtl <- eqtl[eqtl$phenotype_id == pheno_id]
    tmp_eqtl <- tmp_eqtl[!duplicated(tmp_eqtl$start_distance)] # exclude multi-allelic sites
    setkey(tmp_eqtl, variant_id)

    # convert to hg19
    tmp_chr <- gsub(':.*', '', tmp_eqtl$variant_id)
    tmp_pos <- var2pos(tmp_eqtl$variant_id)
    gr_hg38 <- GRanges(seqnames=tmp_chr,
                           ranges=IRanges(tmp_pos, tmp_pos))
    gr_hg38$id_hg38 <- tmp_eqtl$variant_id
    gr_hg19 <- unlist(liftOver(gr_hg38, chain))

    dup_snps <- gr_hg19$id_hg38[duplicated(gr_hg19$id_hg38)]
    gr_hg19 <- gr_hg19[!gr_hg19$id_hg38 %in% dup_snps]
    tmp_eqtl <- tmp_eqtl[gr_hg19$id_hg38]
    mcols(gr_hg19) <- tmp_eqtl
    gr_hg19$hg19_id1 <- paste(seqnames(gr_hg19), start(gr_hg19), gr_hg19$A1, gr_hg19$A2, sep=':')
    gr_hg19$hg19_id2 <- paste(seqnames(gr_hg19), start(gr_hg19), gr_hg19$A2, gr_hg19$A1, sep=':')

    if (length(gr_hg19)<50) {
        system(paste0('rm -r ', outdir))
        susie_meta[pheno_id, 'total_snps'] <- 0
        next
    }
    
    snp_file <- paste0(ct, '/tmp_snps.csv')
    fwrite(data.frame(mcols(gr_hg19)), snp_file)

    # pick which LD file region cover most SNPs
    ovlp_snps <- table(subjectHits(findOverlaps(gr_hg19, LD_gr)))
    LD_id <- as.numeric(names(ovlp_snps)[which.max(ovlp_snps)])
    ld_prefix <- gsub('.gz', '', LD_gr[LD_id]$LD_gz)

    ld_file <- paste0(ct, '/ld_m_tmp.csv')
    cmd <- sprintf('%s --ld_prefix %s --input_file %s --output_file %s', 
                   query_ld_script, ld_prefix, snp_file, ld_file)
    system(cmd)

    ld_m <- tryCatch({
        read.csv(ld_file, row.names=1)
    }, error = function(e) {
        NULL
    })

    if (is.null(ld_m)) {
        system(paste0('rm -r ', outdir))
        susie_meta[pheno_id, 'total_snps'] <- 0
        next
    }
    
    colnames(ld_m) <- rownames(ld_m)

    # harmonize effect allele with ref panel
    eqtl_hg19 <- data.table(data.frame(mcols(gr_hg19)))
    eqtl_hg19 <- eqtl_hg19[eqtl_hg19$hg19_id1 %in% rownames(ld_m) | eqtl_hg19$hg19_id2 %in% rownames(ld_m)]
    eqtl_hg19$id_harmonized <- eqtl_hg19$hg19_id1
    toflip <- eqtl_hg19$hg19_id2 %in% rownames(ld_m)
    eqtl_hg19$id_harmonized[toflip] <- eqtl_hg19$hg19_id2[toflip]
    eqtl_hg19$slope[toflip] <- (-eqtl_hg19$slope[toflip])
    setkey(eqtl_hg19, id_harmonized)        
    eqtl_hg19 <- eqtl_hg19[rownames(ld_m)]        
    susie_meta[pheno_id, 'total_snps'] <- nrow(eqtl_hg19)
        
    # remove outlier snps
    # https://github.com/stephenslab/susieR/issues/182
    z <- eqtl_hg19$slope / eqtl_hg19$slope_se
    condz_in <- kriging_rss(z, data.matrix(ld_m), n=N)
    outlier <- which(abs(condz_in$conditional_dist$z_std_diff)>outlier_thres)
    susie_meta[pheno_id, 'outlier_snps'] <- length(outlier)

    if (length(outlier)>0) {
        eqtl_hg19 <- eqtl_hg19[-outlier]
        ld_m <- ld_m[-outlier, -outlier]            
    }
    susie_meta[pheno_id, 'finemap_snps'] <- nrow(eqtl_hg19)

    d1 <- list()
    d1[['beta']] <- eqtl_hg19$slope
    d1[['varbeta']] <- eqtl_hg19$slope_se**2
    d1[['snp']] <- eqtl_hg19$id_harmonized
    d1[['position']] <- var2pos(eqtl_hg19$id_harmonized)
    d1[['type']] <- 'quant'
    d1[['sdY']] <- Y_sd[pheno_id]
    d1[['LD']] <- data.matrix(ld_m)
    d1[['N']] <- N
        
    susie_meta[pheno_id, 'lambda'] <- estimate_s_rss(d1[['beta']]/sqrt(d1[['varbeta']]), d1[['LD']], N)
    susie_meta[pheno_id, 'alignment'] <- check_alignment(d1, do_plot=F)

    fitted_rss <- tryCatch({
        runsusie(d1, max_iter=500, repeat_until_convergence=F)
    }, error = function(e) {
        NULL
    })
     
    if (is.null(fitted_rss)) {
        system(paste0('rm -r ', outdir))
        susie_meta[pheno_id, 'converge'] <- 0
        next
    }
    
    susie_meta[pheno_id, 'converge'] <- 1
    eqtl_hg19$variable <- 1:nrow(eqtl_hg19)
    res_susie <- merge(x=eqtl_hg19, y=summary(fitted_rss)$vars, 
                       by="variable", all.x=TRUE)
        
    saveRDS(fitted_rss, paste0(outdir, 'rss.rds'))
    write.csv(res_susie, paste0(outdir, 'res_susie.csv'))
    
    tmp <- res_susie$cs
    n_cs <- length(unique(tmp[tmp!=(-1)]))
    susie_meta[pheno_id, 'cs'] <- n_cs
    
    print(pheno_id)
    print(round(susie_meta[pheno_id, ], 3))
        
    if (n_cs==0) {
        system(paste0('rm -r ', outdir))
        next
    }
    
    # plot
    res_susie$cs <- as.factor(res_susie$cs)
    pdf(paste(outdir, 'susie.pdf'), 10, 4)
    p1 <- ggplot(res_susie, aes(x=start_distance, y=-log10(pval_nominal))) + geom_point(alpha=0.8) +
            geom_hline(yintercept = -log10(res_susie$pval_nominal_threshold[1]), linetype = "dashed", color = "red") +
            geom_point(data=subset(res_susie, res_susie$cs!=-1), aes(color=cs), size=2) + theme_classic()
    p2 <- ggplot(res_susie, aes(x=start_distance, y=variable_prob)) + geom_point() + 
            geom_point(data=subset(res_susie, res_susie$cs!=-1), aes(color=cs), size=2) + theme_classic()
    grid.arrange(p1, p2, nrow=2)
    dev.off()
}

write.csv(susie_meta, sprintf('susie_%s.csv', ct))
