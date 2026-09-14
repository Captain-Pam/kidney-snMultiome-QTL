library(data.table)
library(rtracklayer)
library(coloc)
library(arrow)
library(Gviz)
library(ggplot2)
library(grid)
library(gridExtra)
library(ggplotify)

var2pos <- function(var_id) {
    pos <- as.numeric(sapply(strsplit(var_id, ':'), function(x) {x[2]}))
    return(pos)
}

# gene track
exons_df <- readRDS('data/reference/gene_model_hg38.rds')
exons_df$transcript_support_level[exons_df$transcript_support_level=='NA'] <- '100'
exons_df$transcript_support_level[is.na(exons_df$transcript_support_level)] <- '100'
exons_df$transcript_support_level <- as.numeric(exons_df$transcript_support_level)
exons_dt <- data.table(exons_df)
exons_dt <- exons_dt[,.SD[transcript_support_level == min(transcript_support_level)],by = symbol]
exons_df <- data.frame(exons_dt)

top_transcripts <- exons_dt[, .(transcript = unique(transcript)[1:3]), by = gene]
exons_dt <- exons_dt[transcript %in% top_transcripts$transcript & gene %in% top_transcripts$gene]
exons_df <- data.frame(exons_dt)


# load eGFR sum-stat
gr_egfr <- readRDS('data/gwas/egfr_sumstat.rds')
gr_indep <- readRDS('data/gwas/gr_indep.rds')

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

coloc_res_path <- list.files('coloc_res', pattern='.rds', full.names=T)
dir.create('plots')

for (t in 1:length(coloc_res_path)) {
    input_str <- gsub('coloc_res/|.rds', '', coloc_res_path[[t]])
    
    tmp <- strsplit(input_str, '\\.')[[1]]
    i <- tmp[1]
    g <- tmp[2]
    locus <- tmp[3] 
 
    coloc_res <- readRDS(coloc_res_path[t])
    d1 <- coloc_res[[1]]
    d2 <- coloc_res[[2]]
    res <- coloc_res[[3]]
    snps <- d1$snp
    
    top_idx <- which.max(res$result$SNP.PP.H4)[1]
    top_snp <- res$result$snp[top_idx]
    top_cond_pph4 <- res$result$SNP.PP.H4[top_idx]
    
    chr <- gsub(':.*', '', top_snp)

    ########
    # plot #
    ########
    locus_window <- resize(gr_indep[locus], width=500000, fix='center')
    idx <- queryHits(findOverlaps(gr_egfr, locus_window))
    res_l <- data.table(data.frame(mcols(gr_egfr[idx])))
    res_l$color <- ifelse(res_l$id_final %in% snps,'coloc_tested','other')
    res_l$color[res_l$id_final==top_snp] <- 'coloc_lead'
    res_l$P.value <- as.numeric(res_l$P.value)

    res_g <- eqtl_sumstat[[paste0(i, '_', g)]]
    res_g$pos <- var2pos(res_g$variant_id)
    res_g$color <- ifelse(res_g$variant_id %in% snps,'coloc_tested','other')
    res_g$color[res_g$variant_id==top_snp] <- 'coloc_lead'
    
    winL <- min(c(res_l$POS_hg38, res_g$pos))
    winR <- max(c(res_l$POS_hg38, res_g$pos))

    if (g %in% exons_df$symbol) {
        gL <- min(exons_df[exons_df$symbol==g, 'start'])
        gR <- max(exons_df[exons_df$symbol==g, 'end'])
        winL <- min(winL,gL)
        winR <- max(winR,gR)
    }
    win_gr <- GRanges(chr, IRanges(winL, winR))

    ########
    # Gviz #
    ########
    # gene track
    gtrack <- GeneRegionTrack(exons_df, 
                              chromosome = chr, geneSymbol = TRUE, 
                              start=winL, end=winR, just.group='below',
                              fill= '#81D2C7', col= '#81D2C7', shape='arrow',
                              transcriptAnnotation = "symbol", name = "Gene")

    p_g <- as.grob(~plotTracks(list(gtrack), 
                               sizes = 1,                               
                               from = winL, to = winR, 
                               title.width = 0.85, 
                               margin = 27, frame = TRUE,
                               col.axis = "black", 
                               col.frame = "black", 
                               col.grid = "gray90", 
                               fontcolor.title="black", col.axis="black",
                               cex.axis = 0.6, cex.title = 0.4))
    dev.off()

    custom_theme <- theme(plot.margin = margin(1, 27, 5, 1), 
                          legend.position = 'none', legend.justification = c(0, 1),
                          axis.text = element_text(family = "mono"))

    # eGFR
    p1 <- ggplot(res_l, aes(x=POS_hg38, y=-log10(P.value), color=color)) + geom_point(size=1, alpha=0.8) + 
                scale_color_manual(values = c("coloc_tested" = "black", "other" = "gray", "coloc_lead"='red')) + 
                geom_point(data = res_l[color == 'coloc_lead'], color = "red") + 
                geom_point(data = res_l[color == 'coloc_lead'], shape = 1, size = 5, color = "red") + 
                labs(title='eGFR-2025', x=NULL, y='-log10(P)') + theme_classic() + 
                scale_y_continuous(labels = function(x) {formatC(x, width = 6)}) + 
                scale_x_continuous(expand=c(0,0), limits=c(winL,winR)) + custom_theme

    # sc-eQTL
    p2 <- ggplot(res_g, aes(x=pos, y=-log10(pval_nominal), color=color)) + geom_point(size=1, alpha=0.8) + 
                scale_color_manual(values = c("coloc_tested" = "black", "other" = "gray", "coloc_lead"='red')) + 
                geom_point(data = res_g[color == 'coloc_lead'], color = "red") + 
                geom_point(data = res_g[color == 'coloc_lead'], shape = 1, size = 5, color = "red") + 
                labs(title='sc-eQTL', x='hg38_coord', y='-log10(p_nominal)') + theme_classic() + 
                scale_y_continuous(labels = function(x) {formatC(x, width = 6)}) + 
                scale_x_continuous(expand=c(0,0), limits=c(winL,winR)) + custom_theme

    pdf(sprintf('plots/%s.%s.%s.pdf', i, locus, g), 7, 6)
    grid.arrange(p1, p2, p_g, nrow=3)
    dev.off()

    print(t)
}
