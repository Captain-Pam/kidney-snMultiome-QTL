#!/usr/bin/env Rscript
# Purpose: Standardize GWAS summary statistics for scTWAS Stage 2 and FOCUS.
# Inputs:
#   - 2_stage2_focus/2_data/gwas_manifest.tsv with trait-specific column indices.
# Outputs:
#   - 2_stage2_focus/2_data/sumstats_by_trait/<trait>/<trait>.chr<chr>.sumstats.tsv.
#   - 2_stage2_focus/2_data/sumstats_by_trait/<trait>/<trait>.focus_input.tsv.
#   - 2_stage2_focus/3_outputs/summaries/gwas/<trait>.gwas_prepare_summary.tsv.

suppressMessages({
  library(argparser)
  library(data.table)
})

p <- arg_parser("Prepare GWAS summary statistics for scTWAS Stage 2 and FOCUS.")
p <- add_argument(p, "--trait", default = "all",
                  help = "Trait name from gwas_manifest.tsv, comma-separated traits, or 'all'.")
p <- add_argument(p, "--workflow-dir", default = "",
                  help = "Path to the public Main/sctwas directory.")
argv <- parse_args(p)

if (argv$workflow_dir == "") stop("--workflow-dir is required.")

stage_dir <- file.path(argv$workflow_dir, "2_stage2_focus")
manifest_path <- file.path(stage_dir, "2_data/gwas_manifest.tsv")
out_root <- file.path(stage_dir, "2_data/sumstats_by_trait")
summary_dir <- file.path(stage_dir, "3_outputs/summaries/gwas")

if (!file.exists(manifest_path)) stop("Missing GWAS manifest: ", manifest_path)
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

is_missing <- function(x) length(x) == 0 || is.na(x) || identical(x, "") || trimws(as.character(x)) == ""
as_bool <- function(x) tolower(as.character(x)) %in% c("true", "t", "1", "yes", "y")
as_col <- function(x) {
  if (is_missing(x)) return(NA_integer_)
  as.integer(x)
}
get_col <- function(df, idx, label, required = TRUE) {
  idx <- as_col(idx)
  if (is.na(idx)) {
    if (required) stop("Missing required column index for ", label)
    return(NULL)
  }
  if (idx < 1 || idx > ncol(df)) stop("Column index for ", label, " is out of range: ", idx)
  df[[idx]]
}

manifest <- fread(manifest_path, sep = "\t", header = TRUE, data.table = FALSE,
                  colClasses = "character", na.strings = c("", "NA"))
if (argv$trait == "all") {
  traits <- manifest$trait
} else {
  traits <- unlist(strsplit(argv$trait, ",", fixed = TRUE))
}

for (trait_name in traits) {
  row <- manifest[manifest$trait == trait_name, , drop = FALSE]
  if (nrow(row) != 1) stop("Trait must match exactly one manifest row: ", trait_name)
  if (!file.exists(row$input_file)) stop("Missing GWAS input for ", trait_name, ": ", row$input_file)

  header <- as_bool(row$has_header)
  raw <- fread(row$input_file, sep = "\t", header = header, data.table = FALSE, showProgress = FALSE)
  raw_n <- nrow(raw)

  chr_raw <- as.character(get_col(raw, row$chr_col, "CHR"))
  pos <- suppressWarnings(as.integer(get_col(raw, row$pos_col, "BP")))
  a1 <- toupper(as.character(get_col(raw, row$a1_col, "A1")))
  a2 <- toupper(as.character(get_col(raw, row$a2_col, "A2")))

  chr_num <- suppressWarnings(as.integer(sub("^chr", "", chr_raw, ignore.case = TRUE)))
  z_col <- get_col(raw, row$z_col, "Z", required = FALSE)
  if (is.null(z_col)) {
    beta <- suppressWarnings(as.numeric(get_col(raw, row$beta_col, "BETA")))
    se <- suppressWarnings(as.numeric(get_col(raw, row$se_col, "SE")))
    z <- beta / se
  } else {
    z <- suppressWarnings(as.numeric(z_col))
  }

  p_col <- get_col(raw, row$p_col, "P", required = FALSE)
  if (is.null(p_col)) {
    pval <- 2 * stats::pnorm(abs(z), lower.tail = FALSE)
    pval <- pmax(pval, .Machine$double.xmin)
  } else {
    pval <- suppressWarnings(as.numeric(p_col))
  }

  n_col <- get_col(raw, row$n_col, "N", required = FALSE)
  n_value <- if (is.null(n_col)) rep(NA_real_, length(z)) else suppressWarnings(as.numeric(n_col))

  df <- data.frame(
    SNP = paste0("chr", chr_num, ":", pos),
    CHR = chr_num,
    BP = pos,
    A1 = a1,
    A2 = a2,
    Z = z,
    P = pval,
    N = n_value,
    stringsAsFactors = FALSE
  )
  keep <- df$CHR %in% 1:22 & !is.na(df$BP) & !is.na(df$Z) & is.finite(df$Z) &
    !is.na(df$P) & is.finite(df$P) & df$P > 0 &
    !is.na(df$A1) & !is.na(df$A2) & df$A1 != "" & df$A2 != ""
  df <- df[keep, ]
  df <- df[order(df$SNP, df$P, -abs(df$Z), na.last = TRUE), ]
  dup_n <- sum(duplicated(df$SNP))
  df <- df[!duplicated(df$SNP), ]
  df <- df[order(df$CHR, df$BP), ]

  trait_dir <- file.path(out_root, trait_name)
  dir.create(trait_dir, recursive = TRUE, showWarnings = FALSE)
  for (cc in 1:22) {
    out_chr <- file.path(trait_dir, sprintf("%s.chr%d.sumstats.tsv", trait_name, cc))
    stage2_df <- df[df$CHR == cc, c("SNP", "A1", "A2", "Z")]
    write.table(stage2_df, out_chr, sep = "\t", quote = FALSE, row.names = FALSE)
  }

  focus_input <- file.path(trait_dir, paste0(trait_name, ".focus_input.tsv"))
  write.table(df[, c("SNP", "CHR", "BP", "A1", "A2", "Z", "P", "N")],
              focus_input, sep = "\t", quote = FALSE, row.names = FALSE)

  summary <- data.frame(
    trait = trait_name,
    input_file = row$input_file,
    rows_raw = raw_n,
    rows_kept = nrow(df),
    duplicate_snps_removed = dup_n,
    chromosomes_with_rows = length(unique(df$CHR)),
    n_missing_N = sum(is.na(df$N)),
    stage2_dir = trait_dir,
    focus_input = focus_input,
    stringsAsFactors = FALSE
  )
  summary_path <- file.path(summary_dir, paste0(trait_name, ".gwas_prepare_summary.tsv"))
  write.table(summary, summary_path, sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("[%s] raw=%d kept=%d duplicates_removed=%d -> %s\n",
              trait_name, raw_n, nrow(df), dup_n, trait_dir))
}
