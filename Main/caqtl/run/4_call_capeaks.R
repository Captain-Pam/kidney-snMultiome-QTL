#!/usr/bin/env Rscript
# Purpose: Call empirical-FDR caPeaks and retain their full-mode associations.
# Inputs:  merged observed, two permutation and full-mode RASQUAL results.
# Outputs: run/capeaks/<cell_type> caPeak and full-association tables plus summary.
# Run:     Rscript run/4_call_capeaks.R

suppressPackageStartupMessages(library(data.table))
source("src/empirical_fdr.R")

cell_types <- c(
    "CNT_CD_PC", "DCT", "DTL_ATL", "EC", "IC", "Immune",
    "PEC", "PTS", "Podocyte", "Stromal", "TAL", "injPT"
)
rasqual_header <- c(
    "Feature_ID", "rs_ID", "Chromosome", "SNP_position", "Ref_allele",
    "Alt_allele", "Allele_frequency", "HWE_Chi_square_statistic",
    "Imputation_quality_score", "Log_10_BH_Qvalue", "Chi_square_statistic",
    "Effect_size", "Sequencing_mapping_error_rate",
    "Reference_allele_mapping_bias", "Overdispersion", "SNP_ID_within_region",
    "Num_feature_SNPs", "Num_tested_SNPs", "Num_iterations_null",
    "Num_iterations_alternative", "tie_lead_SNP", "Log_likelihood_null",
    "Convergence_status", "Squared_correlation_genotypes_fSNPs",
    "Squared_correlation_genotypes_rSNP"
)
input_root <- "run/merged"
output_root <- "run/capeaks"

read_rasqual <- function(path) {
    result <- fread(path, header = FALSE)
    setnames(result, rasqual_header)
    result <- result[Convergence_status == 0]
    result[, q_value := 10^as.numeric(Log_10_BH_Qvalue)]
    result[!is.finite(q_value) | q_value <= 0 | q_value > 1, q_value := NA_real_]
    result
}

summary_rows <- list()
for (cell_type in cell_types) {
    input_dir <- file.path(input_root, cell_type)
    output_dir <- file.path(output_root, cell_type)
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    observed <- read_rasqual(file.path(input_dir, paste0(cell_type, "_lead_observed.txt")))
    perm1 <- read_rasqual(file.path(input_dir, paste0(cell_type, "_lead_perm1.txt")))[, .(Feature_ID, q_perm1 = q_value)]
    perm2 <- read_rasqual(file.path(input_dir, paste0(cell_type, "_lead_perm2.txt")))[, .(Feature_ID, q_perm2 = q_value)]
    joined <- merge(observed, perm1, by = "Feature_ID")
    joined <- merge(joined, perm2, by = "Feature_ID")

    threshold <- empirical_fdr_threshold(
        joined$q_value, c(joined$q_perm1, joined$q_perm2), alpha = 0.1
    )
    capeaks <- joined[q_value < threshold]
    capeaks[, empirical_fdr_threshold := threshold]
    fwrite(capeaks, file.path(output_dir, paste0(cell_type, "_caPeaks_FDR0.1.tsv")), sep = "\t")

    full <- read_rasqual(file.path(input_dir, paste0(cell_type, "_full_observed.txt")))
    full_capeaks <- full[Feature_ID %in% capeaks$Feature_ID]
    fwrite(
        full_capeaks,
        file.path(output_dir, paste0(cell_type, "_full_caPeak_associations_FDR0.1.tsv")),
        sep = "\t"
    )
    summary_rows[[cell_type]] <- data.table(
        cell_type = cell_type,
        matched_peaks = nrow(joined),
        q_value_threshold = threshold,
        significant_caPeaks = uniqueN(capeaks$Feature_ID),
        full_associations = nrow(full_capeaks)
    )
}

fwrite(rbindlist(summary_rows), file.path(output_root, "caPeak_summary_FDR0.1.tsv"), sep = "\t")
