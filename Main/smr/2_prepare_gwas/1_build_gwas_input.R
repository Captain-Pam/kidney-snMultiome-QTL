#!/usr/bin/env Rscript

# Purpose:
#   Convert one unified hg38 GWAS summary-statistics file into SMR format.
# Inputs:
#   --input        Gzipped or plain TSV with unified GWAS columns.
#   --output       Output eight-column SMR GWAS file.
#   --freq_impute  Optional value used only when FREQ_A1 is missing.
# Output:
#   A resolved SMR file with SNP, A1, A2, freq, b, se, p, and n columns.
# Example:
#   Rscript 1_build_gwas_input.R \
#     --input trait.unified.with_hg38.tsv.gz \
#     --output trait.resolved.smr.tsv

suppressPackageStartupMessages({
  library(argparser)
  library(data.table)
})

required_columns <- c(
  "SNP", "A1", "A2", "FREQ_A1", "BETA", "SE", "P", "N",
  "hg38_chr", "hg38_bp"
)

parse_freq_impute <- function(x) {
  value <- trimws(as.character(x))
  if (value == "" || toupper(value) %chin% c("NA", "NONE", "NULL")) {
    return(NA_real_)
  }
  value <- suppressWarnings(as.numeric(value))
  if (!is.finite(value) || value < 0 || value > 1) {
    stop("--freq_impute must be NA or a number between 0 and 1.")
  }
  value
}

normalize_chr <- function(x) {
  value <- trimws(as.character(x))
  value[toupper(value) %chin% c("", "NA", "NAN", ".")] <- NA_character_
  value <- sub("^chr", "", value, ignore.case = TRUE)
  fifelse(is.na(value), NA_character_, paste0("chr", value))
}

parser <- arg_parser("Convert unified hg38 GWAS statistics to SMR format.")
parser <- add_argument(parser, "--input", help = "Input unified GWAS TSV or TSV.GZ.", default = "")
parser <- add_argument(parser, "--output", help = "Output SMR TSV.", default = "")
parser <- add_argument(parser, "--freq_impute", help = "Optional replacement for missing FREQ_A1.", default = "NA")
args <- parse_args(parser, commandArgs(trailingOnly = TRUE))

if (args$input == "" || !file.exists(args$input)) {
  stop("Input GWAS file does not exist: ", args$input)
}
if (args$output == "") {
  stop("--output is required.")
}

# Read only the fields required by SMR.
if (grepl("\\.gz$", args$input, ignore.case = TRUE)) {
  input_cmd <- paste("gzip -cd", shQuote(args$input))
  input_columns <- names(fread(cmd = input_cmd, nrows = 0))
  missing_columns <- setdiff(required_columns, input_columns)
  if (length(missing_columns) > 0) {
    stop("Missing GWAS columns: ", paste(missing_columns, collapse = ", "))
  }
  gwas <- fread(cmd = input_cmd, select = required_columns)
} else {
  input_columns <- names(fread(args$input, nrows = 0))
  missing_columns <- setdiff(required_columns, input_columns)
  if (length(missing_columns) > 0) {
    stop("Missing GWAS columns: ", paste(missing_columns, collapse = ", "))
  }
  gwas <- fread(args$input, select = required_columns)
}

# Standardize coordinates, alleles, and numeric fields.
gwas[, source_row := .I]
setnames(gwas, "SNP", "rsID")
gwas[, rsID := trimws(as.character(rsID))]
gwas[, c("A1", "A2") := lapply(.SD, function(x) toupper(trimws(as.character(x)))), .SDcols = c("A1", "A2")]
gwas[, hg38_chr := normalize_chr(hg38_chr)]
gwas[, hg38_bp := suppressWarnings(as.integer(hg38_bp))]
gwas[, freq := suppressWarnings(as.numeric(FREQ_A1))]
gwas[, b := suppressWarnings(as.numeric(BETA))]
gwas[, se := suppressWarnings(as.numeric(SE))]
gwas[, p := suppressWarnings(as.numeric(P))]
gwas[, n := suppressWarnings(as.numeric(N))]

freq_impute <- parse_freq_impute(args$freq_impute)
if (is.finite(freq_impute)) {
  gwas[is.na(freq), freq := freq_impute]
}

# Apply the required SMR field rules.
gwas <- gwas[
  !is.na(hg38_chr) & !is.na(hg38_bp) & hg38_bp > 0L &
    is.finite(freq) & freq >= 0 & freq <= 1 &
    !is.na(A1) & A1 != "" & !(A1 %chin% c("NA", "NAN", ".")) &
    !is.na(A2) & A2 != "" & !(A2 %chin% c("NA", "NAN", ".")) &
    is.finite(b) & is.finite(se) & se > 0 &
    is.finite(p) & p >= 0 & p <= 1 &
    is.finite(n) & n > 0
]
gwas[p == 0, p := .Machine$double.xmin]
gwas[, SNP := paste0(hg38_chr, ":", hg38_bp)]

# Resolve duplicated chr:bp IDs deterministically.
base_alleles <- c("A", "C", "G", "T")
gwas[, snv_priority := fifelse(A1 %chin% base_alleles & A2 %chin% base_alleles, 0L, 1L)]
gwas[, rsid_missing_priority := fifelse(is.na(rsID) | rsID == "", 0L, 1L)]
setorder(gwas, SNP, snv_priority, p, -n, rsid_missing_priority, rsID, source_row)
gwas <- gwas[!duplicated(SNP)]

smr <- gwas[, .(
  SNP,
  A1,
  A2,
  freq,
  b,
  se,
  p,
  n
)]

dir.create(dirname(args$output), recursive = TRUE, showWarnings = FALSE)
fwrite(smr, args$output, sep = "\t", quote = FALSE, na = "NA")
message("Saved ", format(nrow(smr), big.mark = ","), " rows to ", args$output)
