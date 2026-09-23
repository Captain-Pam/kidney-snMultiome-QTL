# Purpose: Reproduce the pooled-permutation empirical-FDR threshold used for caPeaks.
# Inputs:  observed and null RASQUAL q-value vectors plus an alpha threshold.
# Output:  largest refined q-value threshold with estimated FDR below alpha.
# Use:     source("src/empirical_fdr.R")

empirical_fdr_threshold <- function(q_observed, q_null, alpha = 0.1) {
    q_observed <- q_observed[!is.na(q_observed)]
    q_null <- q_null[!is.na(q_null)]
    if (length(q_observed) == 0 || length(q_null) == 0) {
        stop("Observed and null q-value vectors must contain finite values")
    }

    threshold <- 0
    for (iteration in seq_len(10)) {
        candidates <- rev(threshold + 0:100 / (100^iteration))
        estimates <- vapply(candidates, function(candidate) {
            mean(q_null < candidate) / mean(q_observed < candidate)
        }, numeric(1))
        accepted <- candidates[estimates < alpha]
        threshold <- max(c(0, accepted), na.rm = TRUE)
    }
    threshold
}
