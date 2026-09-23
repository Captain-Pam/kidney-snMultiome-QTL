#!/usr/bin/env Rscript
# Purpose: Write ordered RASQUAL count, size-factor and covariate binaries.
# Inputs:  QC matrices/peaks from step 6 and covariate tables from step 8.
# Outputs: prepare_data/output/rasqual/<cell_type>/{data,peaks}/.
# Run:     RASQUAL_DIR=software/rasqual Rscript prepare_data/9_prepare_rasqual_inputs.R

suppressPackageStartupMessages({
    library(data.table)
    library(Matrix)
    library(rasqualTools)
})

cell_types <- c(
    "CNT_CD_PC", "DCT", "DTL_ATL", "EC", "IC", "Immune",
    "PEC", "PTS", "Podocyte", "Stromal", "TAL", "injPT"
)
sample_qc_root <- "prepare_data/output/sample_qc"
covariate_root <- "prepare_data/output/covariates"
output_root <- "prepare_data/output/rasqual"
rasqual_dir <- Sys.getenv("RASQUAL_DIR", "software/rasqual")
txt2bin <- file.path(rasqual_dir, "R", "txt2bin.R")

for (cell_type in cell_types) {
    qc_dir <- file.path(sample_qc_root, cell_type)
    out_dir <- file.path(output_root, cell_type)
    data_dir <- file.path(out_dir, "data")
    peak_dir <- file.path(out_dir, "peaks")
    dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(peak_dir, recursive = TRUE, showWarnings = FALSE)

    samples <- readLines(file.path(qc_dir, "samples.qc.txt"))
    peaks <- fread(file.path(qc_dir, "master_peaks.tsv"))
    counts <- readMM(file.path(qc_dir, "counts.qc.mtx"))
    colnames(counts) <- samples
    rownames(counts) <- peaks$peak_id

    count_text <- file.path(data_dir, "counts.mtx.txt")
    size_text <- file.path(data_dir, "size_factors.mtx.txt")
    write.table(
        data.frame(peak_id = peaks$peak_id, as.matrix(counts), check.names = FALSE),
        count_text, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
    )
    size_factors <- rasqualCalculateSampleOffsets(counts, gc_correct = FALSE)
    rownames(size_factors) <- peaks$peak_id
    write.table(
        data.frame(peak_id = peaks$peak_id, as.matrix(size_factors), check.names = FALSE),
        size_text, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
    )

    covariate_files <- list.files(
        file.path(covariate_root, cell_type), pattern = "\\.tsv$", full.names = TRUE
    )
    for (covariate_file in covariate_files) {
        covariates <- fread(covariate_file)
        stopifnot(identical(as.character(covariates$sample_id), samples))
        covariates[, sample_id := NULL]
        covariate_text <- file.path(
            data_dir, paste0("covar_", sub("\\.tsv$", ".txt", basename(covariate_file)))
        )
        write.table(covariates, covariate_text, sep = "\t", quote = FALSE,
                    row.names = FALSE, col.names = FALSE)
        command <- paste(
            "R --vanilla --quiet --args",
            shQuote(count_text), shQuote(size_text), shQuote(covariate_text),
            "<", shQuote(txt2bin)
        )
        status <- system(command)
        if (status != 0) stop("txt2bin failed for ", covariate_file)
    }

    fwrite(peaks, file.path(peak_dir, "master_peaks.tsv"), sep = "\t")
    for (chromosome in paste0("chr", 1:22)) {
        fwrite(peaks[chr == chromosome, .(chr, start, end, peak_id, rasqual_index)],
               file.path(peak_dir, paste0(chromosome, ".peaks.tsv")), sep = "\t")
    }
    file.copy(file.path(qc_dir, "samples.qc.txt"), out_dir, overwrite = TRUE)
    message(cell_type, ": RASQUAL matrices completed")
}
