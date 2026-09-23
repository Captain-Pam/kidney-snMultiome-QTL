#!/usr/bin/env Rscript
# Purpose: Apply caQTL sample/peak QC and prepare ordered chr1 PEER input.
# Inputs:  pseudobulk matrices from step 5 and data/genotype/final_QCed.fam.
# Outputs: prepare_data/output/sample_qc/<cell_type>/ ordered QC files.
# Run:     Rscript prepare_data/6_sample_qc.R

suppressPackageStartupMessages({
    library(data.table)
    library(Matrix)
})

cell_types <- c(
    "CNT_CD_PC", "DCT", "DTL_ATL", "EC", "IC", "Immune",
    "PEC", "PTS", "Podocyte", "Stromal", "TAL", "injPT"
)
pseudobulk_root <- "prepare_data/output/pseudobulk"
output_root <- "prepare_data/output/sample_qc"
fam_file <- "data/genotype/final_QCed.fam"
chromosome_levels <- paste0("chr", 1:22)

dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
genotype_samples <- fread(fam_file, header = FALSE, select = 1)[[1]]

parse_peaks <- function(peak_ids) {
    data.table(
        peak_id = peak_ids,
        chr = sub(":.*", "", peak_ids),
        start = as.integer(sub("-.*", "", sub(".*:", "", peak_ids))),
        end = as.integer(sub(".*-", "", peak_ids))
    )
}

quantile_normalize <- function(x) {
    reference <- rowMeans(apply(x, 2, sort))
    apply(x, 2, function(column) {
        ranks <- rank(column, ties.method = "average")
        approx(seq_along(reference), reference, xout = ranks,
               ties = "ordered")$y
    })
}

inverse_normal_rows <- function(x) {
    t(apply(x, 1, function(values) {
        qnorm((rank(values, ties.method = "average") - 0.5) / length(values))
    }))
}

summary_rows <- list()

for (cell_type in cell_types) {
    input_dir <- file.path(pseudobulk_root, cell_type)
    output_dir <- file.path(output_root, cell_type)
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    counts <- readMM(gzfile(file.path(input_dir, "counts.mtx.gz")))
    peaks <- fread(file.path(input_dir, "peaks.tsv"))
    samples <- fread(file.path(input_dir, "samples.tsv"), header = FALSE)[[1]]
    qc <- fread(file.path(input_dir, "qc_summary.tsv"))
    colnames(counts) <- samples

    peak_table <- parse_peaks(peaks$peak_id)
    autosome_rows <- peak_table$chr %in% chromosome_levels
    peak_table <- peak_table[autosome_rows]
    counts <- counts[autosome_rows, , drop = FALSE]
    chr1_rows <- peak_table$chr == "chr1"
    chr1_counts <- counts[chr1_rows, , drop = FALSE]

    # The study used chr1 pseudobulk depth and chr1 PCA for sample QC.
    depth <- log10(Matrix::colSums(chr1_counts))
    depth_threshold <- mean(depth) - 2 * sd(depth)
    depth_samples <- names(depth)[depth > depth_threshold]
    chr1_counts <- chr1_counts[, depth_samples, drop = FALSE]

    detected <- Matrix::rowSums(chr1_counts > 0) > 0.2 * ncol(chr1_counts)
    pca_counts <- as.matrix(chr1_counts[detected, , drop = FALSE])
    scaled_counts <- t(t(pca_counts) / colSums(pca_counts)) * 1e6
    pca <- prcomp(t(log2(scaled_counts + 1)))
    pc1_z <- as.numeric(scale(pca$x[, 1]))
    pc2_z <- as.numeric(scale(pca$x[, 2]))
    pca_samples <- rownames(pca$x)[abs(pc1_z) <= 2.5 & abs(pc2_z) <= 2.5]
    final_samples <- pca_samples[pca_samples %in% genotype_samples]

    writeLines(final_samples, file.path(output_dir, "samples.qc.txt"))

    counts_qc <- counts[, final_samples, drop = FALSE]
    minimum_samples <- max(5L, ceiling(0.05 * length(final_samples)))
    keep_peak <- Matrix::rowSums(counts_qc >= 3) >= minimum_samples
    counts_qc <- counts_qc[keep_peak, , drop = FALSE]
    peak_qc <- peak_table[keep_peak]

    peak_qc[, chromosome_rank := match(chr, chromosome_levels)]
    ordering <- order(peak_qc$chromosome_rank, peak_qc$start, peak_qc$end)
    peak_qc <- peak_qc[ordering]
    counts_qc <- counts_qc[ordering, , drop = FALSE]
    peak_qc[, rasqual_index := seq_len(.N)]
    peak_qc[, chromosome_rank := NULL]

    fwrite(peak_qc, file.path(output_dir, "master_peaks.tsv"), sep = "\t")
    writeMM(counts_qc, file.path(output_dir, "counts.qc.mtx"))

    for (chromosome in chromosome_levels) {
        fwrite(
            peak_qc[chr == chromosome, .(chr, start, end, peak_id, rasqual_index)],
            file.path(output_dir, paste0(chromosome, ".peaks.tsv")),
            sep = "\t"
        )
    }

    # Prepare the normalized chr1 peak BED consumed by GTEx run_PEER.R.
    peer_rows <- peak_qc$chr == "chr1"
    peer_counts <- as.matrix(counts_qc[peer_rows, , drop = FALSE])
    peer_tpm <- t(t(peer_counts) / colSums(peer_counts)) * 1e6
    peer_keep <- rowSums(peer_tpm >= 0.1) >= 0.2 * ncol(peer_tpm) &
        rowSums(peer_counts >= 6) >= 0.2 * ncol(peer_counts)
    peer_qn <- quantile_normalize(peer_tpm[peer_keep, , drop = FALSE])
    peer_norm <- inverse_normal_rows(peer_qn)
    peer_peaks <- peak_qc[peer_rows][peer_keep]
    peer_bed <- data.table(
        "#chr" = peer_peaks$chr,
        start = peer_peaks$start,
        end = peer_peaks$end,
        phenotype_id = peer_peaks$peak_id
    )
    peer_bed <- cbind(peer_bed, as.data.table(peer_norm))
    setnames(peer_bed, 5:ncol(peer_bed), final_samples)
    fwrite(peer_bed, file.path(output_dir, paste0(cell_type, ".chr1_peak.bed.gz")),
           sep = "\t", compress = "gzip")

    qc_ordered <- qc[match(final_samples, qc$sample_id)]
    fwrite(qc_ordered, file.path(output_dir, "qc_summary.tsv"), sep = "\t")
    summary_rows[[cell_type]] <- data.table(
        cell_type = cell_type,
        samples_before = length(samples),
        samples_after = length(final_samples),
        peaks_after = nrow(peak_qc)
    )
    message(cell_type, ": ", length(final_samples), " donors; ", nrow(peak_qc), " peaks")
}

fwrite(rbindlist(summary_rows), file.path(output_root, "sample_qc_summary.tsv"), sep = "\t")
