library(data.table)
library(ggplot2)
library(gridExtra)

# Step 2 - Sample and gene QC on the pseudobulk counts.
#
# Per cell type: drop low-depth samples (log10 total count < mean - 1.5 SD),
# keep genes expressed in > 20% of samples, drop expression PCA outliers
# (|z(PC1)| or |z(PC2)| > 2.5), and keep only samples with WGS genotypes.
#
# Inputs:
#   prepare_data/anndata/counts/<ct>.csv
#   data/genotype/final_QCed.fam
# Outputs (per cell type, under prepare_data/sample_qc/<ct>/):
#   count_raw.csv, count_tpm.csv, sample_lookup.txt, samples.csv, QC pdfs
#   prepare_data/sample_qc/sample_qc.csv   summary table
#
# Run from the eqtl/ directory.

counts_dir <- 'prepare_data/anndata/counts'
out_dir <- 'prepare_data/sample_qc'
plink_fam <- 'data/genotype/final_QCed.fam'

cts <- gsub('.csv', '', list.files(counts_dir))
geno_samples <- fread(plink_fam)$V1

results <- data.frame(
    cell_type = character(),
    genes_before = integer(),
    samples_before = integer(),
    depth_filter_sample = integer(),
    pca_filter_sample = integer(),
    genes_after = integer(),
    samples_after = integer(),
    samples_ovlp_WGS = integer(),
    stringsAsFactors = FALSE
)

for (ct in cts) {

    dir.create(file.path(out_dir, ct), recursive = TRUE, showWarnings = FALSE)
    m <- read.csv(sprintf('%s/%s.csv', counts_dir, ct), row.names = 1)
    gene_before <- nrow(m)
    sample_before <- ncol(m)

    # filter samples on sequencing depth
    umis <- log10(colSums(m))
    umi_thres <- mean(umis) - 1.5 * sd(umis)
    samples <- names(umis[umis > umi_thres])

    pdf(sprintf('%s/%s/filter_depth.pdf', out_dir, ct), 8, 4)
    par(mfrow = c(1, 2))
    hist(umis, breaks = 100, xlab = 'log10(count)',
         main = sprintf('%d total samples', length(umis)))
    abline(v = umi_thres, col = 'red', lty = 2)
    hist(umis[samples], breaks = 100, xlab = 'log10(count)',
         main = sprintf('%d filtered samples', length(samples)))
    dev.off()

    depth_filter <- ncol(m) - length(samples)
    m <- m[, samples]

    # keep genes expressed in > 20% of samples
    thres <- ncol(m) / 5
    count_thres <- rownames(m)[rowSums(m > 0) > thres]
    m_filter <- m[count_thres, ]
    tpm_filter <- t(t(m_filter) / colSums(m_filter)) * 1e6

    # drop expression PCA outliers
    pca <- prcomp(t(log2(tpm_filter + 1)))
    toplot <- data.frame(pca$x[, 1:2])
    toplot$UMI <- log10(colSums(m_filter))
    toplot$outlier <- ifelse(abs(scale(toplot$PC1)) > 2.5 | abs(scale(toplot$PC2)) > 2.5,
                             'outlier', 'in-dist')

    p1 <- ggplot(toplot, aes(x = PC1, y = PC2, color = UMI)) +
        geom_point(size = 2) + theme_classic()
    p2 <- ggplot(toplot, aes(x = PC1, y = PC2, color = outlier)) +
        geom_point(size = 2) + theme_classic()
    pdf(sprintf('%s/%s/filter_pca.pdf', out_dir, ct), 8, 3)
    grid.arrange(grobs = list(p1, p2), ncol = 2)
    dev.off()

    samples <- rownames(toplot)[toplot$outlier == 'in-dist']
    pca_filter <- ncol(m_filter) - length(samples)

    # keep only samples with WGS genotypes
    final_samples <- samples[samples %in% geno_samples]
    df <- data.frame(sample_ids = final_samples, participant_id = final_samples)
    write.table(df, sprintf('%s/%s/sample_lookup.txt', out_dir, ct),
                sep = '\t', quote = FALSE, row.names = FALSE)

    m_filter <- m_filter[, final_samples]
    tpm_filter <- tpm_filter[, final_samples]

    write.csv(m_filter, sprintf('%s/%s/count_raw.csv', out_dir, ct))
    write.csv(tpm_filter, sprintf('%s/%s/count_tpm.csv', out_dir, ct))
    write.csv(toplot[final_samples, ], sprintf('%s/%s/samples.csv', out_dir, ct))

    results <- rbind(results, data.frame(
        cell_type = ct,
        genes_before = gene_before,
        samples_before = sample_before,
        depth_filter_sample = depth_filter,
        pca_filter_sample = pca_filter,
        genes_after = nrow(m_filter),
        samples_after = length(samples),
        samples_ovlp_WGS = length(final_samples),
        stringsAsFactors = FALSE
    ))
    print(ct)
}

fwrite(results, sprintf('%s/sample_qc.csv', out_dir))
