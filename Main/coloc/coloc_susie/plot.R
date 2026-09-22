library(data.table)
library(rtracklayer)
library(ggplot2)
library(grid)
library(gridExtra)
library(ggplotify)
library(Gviz)

var2pos <- function(var_id) {
    pos <- as.numeric(sapply(strsplit(var_id, ':'), function(x) {x[2]}))
    return(pos)
}

coloc_res <- fread('coloc_0.7.csv')

# make exons_df shorter
# use all transcripts of coloc genes
# for rest of genes, keep only protein coding ones
exons_df <- readRDS('data/reference/gene_model_hg19.rds')
exons_df$transcript_level <- 100
exons_df[which(exons_df$tag=='appris_principal'),]$transcript_level <- 1
exons_df[which(exons_df$tag=='basic'),]$transcript_level <- 2
exons_df[which(exons_df$tag=='appris_candidate_longest'),]$transcript_level <- 3
exons_df[which(exons_df$tag%in%c('appris_candidate', 'exp_conf')),]$transcript_level <- 4
exons_df[which(exons_df$tag%in%c('pseudo_consens', 'readthrough_transcript', 'cds_end_NF', 'cds_start_NF', 'mRNA_end_NF')),]$transcript_level <- 5
exons_dt <- data.table(exons_df)
exons_dt <- exons_dt[,.SD[transcript_level == min(transcript_level)],by = symbol]
exons_df <- data.frame(exons_dt)
a <- exons_df[(exons_df$symbol %in% coloc_res$gene), ]
exons_df <- exons_df[!(exons_df$symbol %in% coloc_res$gene), ]
exons_df <- exons_df[(exons_df$feature=='protein_coding'),]
exons_df <- rbind(a, exons_df)

top_transcripts <- exons_dt[, .(transcript = unique(transcript)[1:3]), by = gene]
exons_dt <- exons_dt[transcript %in% top_transcripts$transcript & gene %in% top_transcripts$gene]
exons_df <- data.frame(exons_dt)


dir.create('plots')

for (i in 1:nrow(coloc_res)) {
    ct <- coloc_res[i]$celltype
    g <- coloc_res[i]$gene
    l <- coloc_res[i]$loci
    top_snp <- coloc_res[i]$top_snp

    res_l <- fread(sprintf('../../susie/gwas_finemapping/outputs/%s/res_susie.csv', l))    
    res_g <- fread(sprintf('../../susie/eqtl_finemapping/%s/%s/res_susie.csv', ct, g))

    res_l$P.value <- as.numeric(res_l$P.value)
    res_l$cs <- as.factor(res_l$cs)
    res_l$P.value[res_l$P.value==0] <- 1e-300

    res_g$pval_nominal <- as.numeric(res_g$pval_nominal)
    res_g$cs <- as.factor(res_g$cs)
    res_g$pos <- var2pos(res_g$id_harmonized)

    # plot gene track
    chr <- paste0('chr', res_l$CHR[1])
    winL <- min(c(res_l$POS, res_g$pos)) - 1000
    winR <- max(c(res_l$POS, res_g$pos)) + 1000

    if (g %in% exons_df$symbol) {
        winL <- min(winL, min(exons_df[exons_df$symbol==g,'start'])) - 1000
        winR <- max(winR, max(exons_df[exons_df$symbol==g,'end'])) + 1000
    }

    gtrack <- GeneRegionTrack(exons_df, 
                              chromosome = chr, geneSymbol = TRUE, 
                              start=winL, end=winR, just.group='below',
                              fill= '#81D2C7', col= '#81D2C7', shape='arrow',
                              transcriptAnnotation = "symbol", name = "Gene")
    p_g <- as.grob(~plotTracks(list(gtrack),
                              from = winL, to = winR, 
                              title.width = 1, 
                              margin = 27, frame = TRUE,
                              col.axis = "black", 
                              col.frame = "black", 
                              col.grid = "gray90", 
                              fontcolor.title="black", col.axis="black",
                              cex.axis = 0.6, cex.axis = 0.6))
    dev.off()

    custom_theme <- theme(plot.margin = margin(1, 27, 5, 1), 
                          legend.position = 'none', legend.justification = c(0, 1),
                          axis.text = element_text(family = "mono"))

    p1 <- ggplot(res_l, aes(x=POS, y=-log10(P.value))) + geom_point(size=1, alpha=0.8) + 
                geom_point(data=subset(res_l, res_l$cs!=-1), aes(color=cs)) + theme_classic() + 
                geom_point(data = res_l[res_l$id_harmonized==top_snp,], shape = 1, size = 5, color = "red") + 
                labs(title='eGFR-2025', x=NULL, y='-log10(P)') +
                scale_y_continuous(labels = function(x) {formatC(x, width = 6)}) + 
                scale_x_continuous(expand=c(0,0), limits=c(winL,winR)) + custom_theme

    p2 <- ggplot(res_l, aes(x=POS, y=variable_prob)) + geom_point(size=1, alpha=0.8) + 
                geom_point(data=subset(res_l, res_l$cs!=-1), aes(color=cs)) + theme_classic() + 
                geom_point(data = res_l[res_l$id_harmonized==top_snp,], shape = 1, size = 5, color = "red") + 
                labs(x=NULL, y='PIP') +
                scale_y_continuous(labels = function(x) {formatC(x, width = 6)}) + 
                scale_x_continuous(expand=c(0,0), limits=c(winL,winR)) + custom_theme

    p3 <- ggplot(res_g, aes(x=pos, y=-log10(pval_nominal))) + geom_point(size=1, alpha=0.8) + 
                geom_point(data=subset(res_g, res_g$cs!=-1), aes(color=cs)) + theme_classic() + 
                geom_point(data = res_g[res_g$id_harmonized==top_snp,], shape = 1, size = 5, color = "red") + 
                labs(title='sc-eQTL', x=NULL, y='-log10(Q)') +
                scale_y_continuous(labels = function(x) {formatC(x, width = 6)}) + 
                scale_x_continuous(expand=c(0,0), limits=c(winL,winR)) + custom_theme

    p4 <- ggplot(res_g, aes(x=pos, y=variable_prob)) + geom_point(size=1, alpha=0.8) + 
                geom_point(data=subset(res_g, res_g$cs!=-1), aes(color=cs)) + theme_classic() + 
                geom_point(data = res_g[res_g$id_harmonized==top_snp,], shape = 1, size = 5, color = "red") + 
                labs(x='hg19 coord', y='PIP') +
                scale_y_continuous(labels = function(x) {formatC(x, width = 6)}) + 
                scale_x_continuous(expand=c(0,0), limits=c(winL,winR)) + custom_theme

    pdf(sprintf('plots/%s.%s.%s.pdf', ct, l, g), 7, 6)
    grid.arrange(p1, p2, p3, p4, p_g, nrow=5)
    dev.off()

    print(i)
}
