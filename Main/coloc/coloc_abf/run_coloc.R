library(data.table)
library(rtracklayer)
library(coloc)
library(arrow)

var2pos <- function(var_id) {
    pos <- as.numeric(sapply(strsplit(var_id, ':'), function(x) {x[2]}))
    return(pos)
}

# Load shared data
indep <- readRDS('data/gwas/gr_indep.rds')
loci_gr <- resize(indep, width=100000, fix='center')
egfr_list <- readRDS('data/gwas/egfr_list.rds')
genes_list <- readRDS('gene_gr.rds')  # gene window

cts <- gsub('.parquet', '', list.files('data/eqtl/tensorqtl/all_parquets/'))

dir.create('outputs', showWarnings = FALSE)
dir.create('coloc_res', showWarnings = FALSE)

N_egfr <- 1.7e6
all_output_df <- list()

for (i in cts) {

    message("Processing: ", i)

    ######################
    # eQTL split by gene #
    ######################
    eqtl_all <- data.table(read_parquet(sprintf('data/eqtl/tensorqtl/all_parquets/%s.parquet', i)))
    eqtl_all$tmp <- sub(":[^:]+:[^:]+$", "", eqtl_all$variant_id)
    eqtl_all$id1 <- paste0(eqtl_all$tmp, ':', eqtl_all$A2, ':', eqtl_all$A1)
    eqtl_all$id2 <- paste0(eqtl_all$tmp, ':', eqtl_all$A1, ':', eqtl_all$A2)
    
    egenes <- unique(eqtl_all[eqtl_all$egene_qval < 0.05, ]$phenotype_id)
    eqtl_all_list <- split(eqtl_all, eqtl_all$phenotype_id)

    ###########
    # pheno m #
    ###########
    Y <- data.frame(fread(sprintf('data/eqtl/expression_bed/%s.expression.bed.gz', i)))
    rownames(Y) <- Y$gene_id
    Y <- Y[, 5:ncol(Y)]
    N_eqtl <- ncol(Y)
    Y_sd <- apply(Y, 1, sd)

    #######################
    # gene-locus matching #
    #######################
    gene_gr <- genes_list[[i]]
    ovlp_df <- data.frame(findOverlaps(gene_gr, loci_gr))
    output_df <- data.frame(locus = names(loci_gr)[ovlp_df[,2]],
                            gene = names(gene_gr)[ovlp_df[,1]],
                            nsnps = NA,
                            PP.H0.abf = NA,
                            PP.H1.abf = NA,
                            PP.H2.abf = NA,
                            PP.H3.abf = NA,
                            PP.H4.abf = NA,
                            top_snp = NA,
                            cond_pph4 = NA,
                            REF = NA,
                            ALT = NA,
                            eQTL_z = NA,
                            eGFR_z = NA
                            )

    for (j in 1:nrow(output_df)) {

        l <- output_df[j, 'locus']
        g <- output_df[j, 'gene']

        if (!g %in% names(eqtl_all_list)) next

        # egfr sum stats
        egfr_tmp <- as.data.table(mcols(egfr_list[[l]]))
        setkey(egfr_tmp, 'vid')
        
        # find overlapping snps, flip for eqtl
        eqtl_tmp <- eqtl_all_list[[g]]
        eqtl_tmp <- eqtl_tmp[(eqtl_tmp$id1 %in% egfr_tmp$vid) | (eqtl_tmp$id2 %in% egfr_tmp$vid)]
        eqtl_tmp[, id_final := ifelse(id1 %in% egfr_tmp$vid, id1, id2)]
        eqtl_tmp <- unique(eqtl_tmp, by = "id_final")
        eqtl_tmp[, slope := ifelse(id_final == id1, slope, -slope)]        
        setkey(eqtl_tmp, 'id_final')
        
        ovlp_snp <- eqtl_tmp$id_final
        output_df[j, 'nsnps'] = length(ovlp_snp)
        if (length(ovlp_snp) < 20) next

        egfr_tmp <- egfr_tmp[ovlp_snp]
        eqtl_tmp <- eqtl_tmp[ovlp_snp]

        # coloc input
        d1 <- list(beta = eqtl_tmp$slope,
                   varbeta = eqtl_tmp$slope_se^2,
                   snp = eqtl_tmp$id_final,
                   position = var2pos(eqtl_tmp$variant_id),
                   type = 'quant',
                   sdY = Y_sd[g])

        d2 <- list(beta = egfr_tmp$BETA,
                   varbeta = egfr_tmp$SE^2,
                   snp = egfr_tmp$vid,
                   position = egfr_tmp$POS_hg38,
                   type = 'quant',
                   N = N_egfr,
                   MAF = egfr_tmp$MAF)

        # coloc
        res <- suppressWarnings(coloc.abf(dataset1 = d1, dataset2 = d2))
        top_idx <- which.max(res$result$SNP.PP.H4)[1]
        top_snp <- res$result$snp[top_idx]
        top_cond_pph4 <- res$result$SNP.PP.H4[top_idx]

        output_df[j, 4:8] <- res$summary[2:6]
        output_df[j, 'top_snp'] <- top_snp
        output_df[j, 'cond_pph4'] <- top_cond_pph4

        # additional info
        df <- res$result
        rownames(df) <- df$snp
        output_df[j, 'REF'] <- egfr_tmp[top_snp]$REF
        output_df[j, 'ALT'] <- egfr_tmp[top_snp]$ALT
        output_df[j, 'eQTL_z'] <- df[top_snp, 'z.df1']
        output_df[j, 'eGFR_z'] <- df[top_snp, 'z.df2']

        # Save RDS when PP.H4.abf >= 0.7.
        if (output_df[j, 'PP.H4.abf'] >= 0.7) {
            save_list <- list(eqtl = d1, egfr = d2, res = res)
            saveRDS(save_list, sprintf('coloc_res/%s.%s.%s.rds', i, g, l))            
        }
    }

    output_df$celltype <- i
    write.csv(output_df, sprintf('outputs/%s.csv', i), row.names = FALSE)
    all_output_df[[i]] <- output_df
}

# Combine and filter for PP.H4.abf >= 0.7.
combined <- rbindlist(all_output_df)
combined <- combined[PP.H4.abf >= 0.7]
fwrite(combined, 'coloc_0.7_z.csv')
