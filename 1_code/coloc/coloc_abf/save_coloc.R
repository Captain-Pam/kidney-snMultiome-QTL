library(data.table)
library(rtracklayer)
library(coloc)
library(arrow)
library(Gviz)

var2pos <- function(var_id) {
    pos <- as.numeric(sapply(strsplit(var_id, ':'), function(x) {x[2]}))
    return(pos)
}

##############################
# 100k around each eGFR loci #
##############################
egfr_list <- readRDS('data/gwas/egfr_list.rds') # sumstats for each locus

cts <- c('CNT_CD_PC', 'DCT', 'Endothelial', 'Immune', 'Intercalated', 
         'PT', 'Stromal', 'Thick_limb', 'Thin_limb')

coloc_pairs <- fread('coloc_0.8.csv')
coloc_pairs[, REF := as.character(NA)]
coloc_pairs[, ALT := as.character(NA)]
coloc_pairs[, eQTL_z := as.numeric(NA)]
coloc_pairs[, eGFR_z := as.numeric(NA)]

dir.create('coloc_res')
N <- 1.7e6

for (i in cts) {
    
    tmp_pairs <- coloc_pairs[coloc_pairs$celltype==i]
    if (nrow(tmp_pairs)==0) {next}

    ######################
    # eQTL split by gene #
    ######################
    eqtl_all <- data.table(read_parquet(sprintf('data/eqtl/tensorqtl/all_parquets/%s.parquet', i)))
    eqtl_all_list <- split(eqtl_all, eqtl_all$phenotype_id)

    ###########
    # pheno m #
    ###########
    Y <- data.frame(fread(sprintf('data/eqtl/expression_bed/psd_DESeq2_%s_2.bed.gz', i)))
    rownames(Y) <- Y$feature
    Y <- Y[,5:ncol(Y)]
    N <- ncol(Y)
    Y_sd <-  apply(Y, 1, sd)
    
    ##################
    # loop over loci #
    ##################
    for (t in 1:nrow(tmp_pairs)) {
        locus <- tmp_pairs[t]$locus
        g <- tmp_pairs[t]$gene

        egfr_tmp <- as.data.table(mcols(egfr_list[[locus]]))
        egfr_tmp <- egfr_tmp[!duplicated(egfr_tmp$id_final)]
        setkey(egfr_tmp, id_final)

        eqtl_tmp <- eqtl_all_list[[g]] # eqtl sum stats
        setkey(eqtl_tmp, 'variant_id')
        
        ovlp_snp <- intersect(eqtl_tmp$variant_id, egfr_tmp$id_final)        
        if (length(ovlp_snp)<20) { next }
        egfr_tmp <- egfr_tmp[ovlp_snp]
        eqtl_tmp <- eqtl_tmp[ovlp_snp]

        # allele flip
        allele_flip <- which(eqtl_tmp$A1!=egfr_tmp$ALT)
        eqtl_tmp[allele_flip]$slope <- (-eqtl_tmp[allele_flip]$slope)
        
        d1 <- list()
        d1[['beta']] <- eqtl_tmp$slope
        d1[['varbeta']] <- eqtl_tmp$slope_se**2
        d1[['snp']] <- eqtl_tmp$variant_id
        d1[['position']] <- var2pos(eqtl_tmp$variant_id)
        d1[['type']] <- 'quant'
        d1[['sdY']] <- Y_sd[g]
            
        d2 <- list()
        d2[['beta']] <- egfr_tmp$BETA
        d2[['varbeta']] <- egfr_tmp$SE**2
        d2[['snp']] <- egfr_tmp$id_final
        d2[['position']] <- egfr_tmp$POS_hg38
        d2[['type']] <- 'quant'      
        d2[['N']] <- N
        d2[['MAF']] <- egfr_tmp$MAF
        
        res <- suppressWarnings(coloc.abf(dataset1=d1, dataset2=d2))
        top_idx <- which.max(res$result$SNP.PP.H4)[1]
        top_snp <- res$result$snp[top_idx]
        top_cond_pph4 <- res$result$SNP.PP.H4[top_idx]

        save_list <- list(eqtl=d1, egfr=d2, res=res)
        saveRDS(save_list, sprintf('coloc_res/%s.%s.%s.rds', i, g, locus))

        # additional info
        df <- res$results
        rownames(df) <- df$snp
        idx <- which((coloc_pairs$celltype==i) & (coloc_pairs$locus==locus) & (coloc_pairs$gene==g))
        coloc_pairs[idx, 'REF'] <- egfr_tmp[top_snp]$REF
        coloc_pairs[idx, 'ALT'] <- egfr_tmp[top_snp]$ALT
        coloc_pairs[idx, 'eQTL_z'] <- df[top_snp, 'z.df1']
        coloc_pairs[idx, 'eGFR_z'] <- df[top_snp, 'z.df2']

        print(t)
    }
}

fwrite(coloc_pairs, 'coloc_0.8_z.csv')
