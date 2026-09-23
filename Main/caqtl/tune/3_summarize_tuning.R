#!/usr/bin/env Rscript
# Purpose: Summarize converged lead-mode RASQUAL fits across tuning models.
# Inputs:  tune/output/<cell_type>/<model>/rasqual_*.txt.
# Outputs: tune/tuning_summary.tsv and tune/final_covariate_choice.tsv.
# Run:     Rscript tune/3_summarize_tuning.R

suppressPackageStartupMessages(library(data.table))

cell_types <- c(
    "CNT_CD_PC", "DCT", "DTL_ATL", "EC", "IC", "Immune",
    "PEC", "PTS", "Podocyte", "Stromal", "TAL", "injPT"
)
models <- c(
    paste0("peer", c(0, 1, 2, 3, 4, 5, 10), "_pc0"),
    paste0("peer2_pc", c(1, 2, 3, 4, 5, 10, 15))
)
output_root <- "tune/output"

summaries <- list()
row_index <- 1L
for (cell_type in cell_types) {
    for (model in models) {
        files <- file.path(
            output_root, cell_type, model,
            paste0("rasqual_", cell_type, "_chr", 1:22, "_", model, ".txt")
        )
        missing_files <- files[!file.exists(files)]
        if (length(missing_files) > 0) {
            stop("Missing tuning result: ", missing_files[[1]])
        }
        result <- rbindlist(lapply(files, function(path) {
            fread(path, header = FALSE, select = c(10, 11, 23), showProgress = FALSE)
        }))
        result <- result[V23 == 0]
        rasqual_q <- 10^as.numeric(result$V10)
        nominal_p <- pchisq(as.numeric(result$V11), df = 1, lower.tail = FALSE)
        global_q <- p.adjust(nominal_p, method = "BH")

        summaries[[row_index]] <- data.table(
            cell_type = cell_type,
            model = model,
            converged_fits = nrow(result),
            rasqual_bh_0.1 = sum(p.adjust(rasqual_q, "BH") < 0.1, na.rm = TRUE),
            global_bh_0.01 = sum(global_q < 0.01, na.rm = TRUE),
            global_bh_0.05 = sum(global_q < 0.05, na.rm = TRUE),
            global_bh_0.1 = sum(global_q < 0.1, na.rm = TRUE)
        )
        row_index <- row_index + 1L
    }
}

fwrite(rbindlist(summaries), "tune/tuning_summary.tsv", sep = "\t")
fwrite(
    data.table(cell_type = cell_types, n_peer = 2L, n_genotype_pc = 0L),
    "tune/final_covariate_choice.tsv", sep = "\t"
)
