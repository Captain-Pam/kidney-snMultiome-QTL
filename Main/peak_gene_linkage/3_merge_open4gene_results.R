# Purpose: Merge Open4Gene chunk results and calculate cell-type-specific BH FDR values.
# Inputs: Uncompressed res_chunk_NNN.tsv files in the directory supplied by --result-dir.
# Outputs: combine_peak_gene.txt.gz and hit_fdr005_peak_gene.txt.gz in --result-dir.
# Example:
# Rscript Main/peak_gene_linkage/3_merge_open4gene_results.R --result-dir results/open4gene

parser <- argparser::arg_parser("Merge Open4Gene chunk results")
parser <- argparser::add_argument(
  parser,
  "--result-dir",
  help = "Directory containing res_chunk_NNN.tsv files",
  type = "character"
)
args <- argparser::parse_args(parser)

if (is.na(args$result_dir)) {
  stop("Missing required option: --result-dir", call. = FALSE)
}

result_files <- list.files(
  args$result_dir,
  pattern = "^res_chunk_[0-9]{3}\\.tsv$",
  full.names = TRUE
)
if (length(result_files) == 0L) {
  stop("No res_chunk_NNN.tsv files found in --result-dir", call. = FALSE)
}

required_columns <- c(
  "Celltype",
  "hurdle.Res.zero.p",
  "hurdle.Res.count.p"
)
pvalue_columns <- c("hurdle.Res.zero.p", "hurdle.Res.count.p")

# Validate each chunk before combining and retain its basename for provenance.
chunk_results <- lapply(result_files, function(result_file) {
  chunk_result <- data.table::fread(result_file)
  source_name <- basename(result_file)

  missing_columns <- setdiff(required_columns, names(chunk_result))
  if (length(missing_columns) > 0L) {
    stop(
      sprintf(
        "Missing required columns in %s: %s",
        source_name,
        paste(missing_columns, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  invalid_celltype_rows <- which(
    is.na(chunk_result$Celltype) |
      trimws(as.character(chunk_result$Celltype)) == ""
  )
  if (length(invalid_celltype_rows) > 0L) {
    stop(
      sprintf(
        "Missing or blank Celltype in %s at row(s): %s",
        source_name,
        paste(invalid_celltype_rows, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  for (pvalue_column in pvalue_columns) {
    p_values <- chunk_result[[pvalue_column]]
    if (is.logical(p_values) && all(is.na(p_values))) {
      p_values <- as.numeric(p_values)
      data.table::set(chunk_result, j = pvalue_column, value = p_values)
    }
    if (!is.numeric(p_values)) {
      stop(
        sprintf(
          "Non-numeric %s in %s at row(s): %s",
          pvalue_column,
          source_name,
          paste(seq_len(nrow(chunk_result)), collapse = ", ")
        ),
        call. = FALSE
      )
    }

    invalid_pvalue_rows <- which(
      !is.na(p_values) &
        (!is.finite(p_values) | p_values < 0 | p_values > 1)
    )
    if (length(invalid_pvalue_rows) > 0L) {
      stop(
        sprintf(
          "Invalid %s in %s at row(s): %s; expected finite values in [0, 1] or NA",
          pvalue_column,
          source_name,
          paste(invalid_pvalue_rows, collapse = ", ")
        ),
        call. = FALSE
      )
    }
  }

  chunk_result[, source_file := basename(result_file)]
  chunk_result
})
combined <- data.table::rbindlist(chunk_results, use.names = TRUE, fill = TRUE)

# Preserve missing p-values while adjusting observed values within each cell type.
adjust_bh <- function(p_values) {
  adjusted <- rep(NA_real_, length(p_values))
  observed <- !is.na(p_values)
  adjusted[observed] <- stats::p.adjust(p_values[observed], method = "BH")
  adjusted
}

combined[, zero_FDR := adjust_bh(hurdle.Res.zero.p), by = Celltype]
combined[, count_FDR := adjust_bh(hurdle.Res.count.p), by = Celltype]

hits <- combined[
  (!is.na(zero_FDR) & zero_FDR < 0.05) |
    (!is.na(count_FDR) & count_FDR < 0.05)
]

# Write the complete table and the fixed-FDR hit subset.
data.table::fwrite(
  combined,
  file.path(args$result_dir, "combine_peak_gene.txt.gz"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
data.table::fwrite(
  hits,
  file.path(args$result_dir, "hit_fdr005_peak_gene.txt.gz"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
