library(rtracklayer)

# make a rds file of genomic ranges of the LD blocks

UKBB_path <- 'data/reference/ukbb_ld_scores' # UKBB EUR LD blocks (Broad/alkesgroup), hg19: chr_start_end.gz + .npz
LD_gz <- list.files(UKBB_path, pattern=".gz", full.names=T)
pos_df <- t(data.frame(strsplit(gsub('.*/|.gz', '', LD_gz), '_')))
LD_gr <- GRanges(seqnames=pos_df[,1],
                 ranges=IRanges(as.numeric(pos_df[,2]), as.numeric(pos_df[,3])))
LD_gr$LD_gz <- LD_gz
LD_gr$LD_npz <- gsub('.gz', '.npz', LD_gz)
LD_gr <- sort(LD_gr)
saveRDS(LD_gr, 'LD_gr.rds')