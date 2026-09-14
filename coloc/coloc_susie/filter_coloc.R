library(data.table)
library(GenomicRanges)
library(arrow)

var2pos <- function(var_id) {
    pos <- as.numeric(sapply(strsplit(var_id, ':'), function(x) {x[2]}))
    return(pos)
}

# filter PIP>0.8
outputs <- fread('coloc_all.csv')
outputs <- unique(outputs) # duplicates because i computed some loci multiple times.
outputs <- outputs[which(outputs$PP.H4.abf>0.8),]
outputs[, top_snp := as.character(NA)]
outputs[, cond_pph4 := as.numeric(NA)]

outputs[, REF := as.character(NA)]
outputs[, ALT := as.character(NA)]
outputs[, eQTL_z := as.numeric(NA)]
outputs[, eGFR_z := as.numeric(NA)]

eqtl_path <- '../../susie/eqtl_finemapping/'
egfr_path <- '../../susie/gwas_finemapping/'

# load eQTL sum-stat
eqtl_sumstat <- list()
cts <- gsub('.parquet', '', list.files('data/eqtl/tensorqtl/all_parquets'))
for (ct in cts) {
    ct_sumstat <- as.data.table(read_parquet(sprintf('data/eqtl/tensorqtl/fdr/%s.parquet', ct)))
    dt_split <- split(ct_sumstat, by = "phenotype_id", keep.by = TRUE)
    names(dt_split) <- paste0(ct, '_', names(dt_split))
    eqtl_sumstat <- c(eqtl_sumstat, dt_split)
    print(ct)
}

# load eGFR sum-stat
egfr_list <- readRDS(sprintf('%s/egfr_list_1m.rds', egfr_path))

# add information of top snp
for (i in 1:nrow(outputs)) {
    locus <- outputs[i]$loci
    gene <- outputs[i]$gene
    celltype <- outputs[i]$celltype
    cs1 <- outputs[i]$idx1
    cs2 <- outputs[i]$idx2
    
    susie.res <- readRDS(sprintf('outputs/%s_%s_%s/coloc_res.rds', celltype, locus, gene))
    tmp <- susie.res$summary
    res_df <- data.frame(susie.res$results)

    if (nrow(tmp)==1) {
        col_name <- 'SNP.PP.H4.abf'
        id_max <- which.max(res_df[, col_name])
        top_snp <- res_df[id_max, 'snp']
        outputs[i, 'top_snp'] <- top_snp
        outputs[i, 'cond_pph4'] <- res_df[id_max, col_name]
    } else {
        rowid <-  which(tmp$idx1==cs1 & tmp$idx2==cs2)
        col_name <- paste0('SNP.PP.H4.row', rowid)
        id_max <- which.max(res_df[, col_name])
        top_snp <- res_df[id_max, 'snp']
        outputs[i, 'top_snp'] <- top_snp
        outputs[i, 'cond_pph4'] <- res_df[id_max, col_name]        
    }

    outputs[i, 'REF'] <- strsplit(top_snp, ':')[[1]][3]
    outputs[i, 'ALT'] <- strsplit(top_snp, ':')[[1]][4]

    # eQTL
    t1 <- fread(sprintf('%s/%s/%s/res_susie.csv', eqtl_path, celltype, gene))
    setkey(t1, id_harmonized)
    vid <- t1[top_snp]$variant_id
    tmp <- eqtl_sumstat[[paste0(celltype,'_',gene)]]
    tmp <- tmp[tmp$variant_id==vid]
    if (tmp$A1==outputs[i]$ALT) {
        outputs[i, 'eQTL_z'] <- tmp$slope/tmp$slope_se
    } else {
        outputs[i, 'eQTL_z'] <- (-tmp$slope)/tmp$slope_se        
    }
    
    # eGFR
    tmp <- egfr_list[[locus]]
    tmp <- tmp[tmp$POS==var2pos(top_snp)]
    if (tmp$ALT==outputs[i]$ALT) {
        outputs[i, 'eGFR_z'] <- tmp$Zscore
    } else {
        outputs[i, 'eGFR_z'] <- (-tmp$Zscore)
    }

    print(i)
}

fwrite(outputs, 'coloc_0.8.csv')
