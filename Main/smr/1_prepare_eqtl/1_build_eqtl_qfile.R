#!/usr/bin/env Rscript

# Purpose:
#   Convert one tensorQTL eQTL parquet file into a cleaned SMR qfile.
# Inputs:
#   --input       tensorQTL parquet file for one cell type.
#   --annotation  Gene annotation table with seqnames, TSS, strand, gene_id,
#                 and gene_symbol columns.
#   --output      Output SMR qfile in tab-separated format.
# Output:
#   A 14-column SMR qfile with probes lacking valid chromosome or TSS
#   annotations removed.
# Example:
#   Rscript 1_build_eqtl_qfile.R \
#     --input celltype.parquet \
#     --annotation gene_annotation.tsv.gz \
#     --output celltype.qfile.tsv

suppressPackageStartupMessages({
  library(argparser)
  library(arrow)
  library(data.table)
})

qtl_columns <- c(
  "phenotype_id", "variant_id", "af", "A1", "A2",
  "slope", "slope_se", "pval_nominal"
)

annotation_columns <- c("seqnames", "TSS", "strand", "gene_id", "gene_symbol")

qfile_columns <- c(
  "SNP", "Chr", "BP", "A1", "A2", "Freq", "Probe", "Probe_Chr",
  "Probe_bp", "Gene", "Orientation", "b", "se", "p"
)

read_annotation <- function(path) {
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    fread(cmd = paste("gzip -cd", shQuote(path)))
  } else {
    fread(path)
  }
}

to_smr_chr <- function(x) {
  cleaned <- toupper(sub("^chr", "", trimws(as.character(x)), ignore.case = TRUE))
  valid <- c(as.character(1:22), "X", "Y", "XY", "M", "MT")
  code <- suppressWarnings(as.integer(cleaned))
  code[cleaned == "X"] <- 23L
  code[cleaned == "Y"] <- 24L
  code[cleaned == "XY"] <- 25L
  code[cleaned %chin% c("M", "MT")] <- 26L
  code[!(cleaned %chin% valid)] <- NA_integer_
  code
}

parser <- arg_parser("Convert one tensorQTL parquet file to a cleaned SMR qfile.")
parser <- add_argument(parser, "--input", help = "Input tensorQTL parquet file.", default = "")
parser <- add_argument(parser, "--annotation", help = "Gene annotation TSV or TSV.GZ file.", default = "")
parser <- add_argument(parser, "--output", help = "Output SMR qfile TSV.", default = "")
args <- parse_args(parser, commandArgs(trailingOnly = TRUE))

if (args$input == "" || !file.exists(args$input)) {
  stop("Input parquet file does not exist: ", args$input)
}
if (args$annotation == "" || !file.exists(args$annotation)) {
  stop("Annotation file does not exist: ", args$annotation)
}
if (args$output == "") {
  stop("--output is required.")
}

# Read the required eQTL and annotation fields.
qtl <- as.data.table(read_parquet(
  args$input,
  col_select = qtl_columns,
  as_data_frame = TRUE
))
missing_qtl <- setdiff(qtl_columns, names(qtl))
if (length(missing_qtl) > 0) {
  stop("Missing tensorQTL columns: ", paste(missing_qtl, collapse = ", "))
}

annotation <- read_annotation(args$annotation)
missing_annotation <- setdiff(annotation_columns, names(annotation))
if (length(missing_annotation) > 0) {
  stop("Missing annotation columns: ", paste(missing_annotation, collapse = ", "))
}

# Build gene-ID and unique-symbol annotation maps.
annotation[, gene_id_key := sub("\\..*$", "", as.character(gene_id))]
annotation[, `:=`(
  Probe_Chr = to_smr_chr(seqnames),
  Probe_bp = suppressWarnings(as.integer(TSS)),
  Gene = as.character(gene_symbol),
  Orientation = as.character(strand)
)]

id_map <- annotation[!duplicated(gene_id_key), .(
  Probe_join = gene_id_key,
  Probe_Chr,
  Probe_bp,
  Gene,
  Orientation
)]

symbol_counts <- annotation[
  !is.na(gene_symbol) & trimws(as.character(gene_symbol)) != "",
  .N,
  by = gene_symbol
]
unique_symbols <- symbol_counts[N == 1, as.character(gene_symbol)]
symbol_map <- annotation[as.character(gene_symbol) %chin% unique_symbols, .(
  Probe_join = as.character(gene_symbol),
  Probe_Chr,
  Probe_bp,
  Gene,
  Orientation
)]

# Parse hg38 chr:bp variant IDs while retaining the original ID as SNP.
qtl[, variant_id := as.character(variant_id)]
variant_parts <- tstrsplit(qtl$variant_id, ":", fixed = TRUE)
if (length(variant_parts) != 2L) {
  stop("variant_id must use the hg38 chr:bp format.")
}
variant_chr <- to_smr_chr(variant_parts[[1]])
variant_bp <- suppressWarnings(as.integer(variant_parts[[2]]))
if (any(is.na(variant_chr) | is.na(variant_bp) | variant_bp <= 0L)) {
  stop("variant_id contains an invalid chromosome or position.")
}

# Attach probe annotations, using a unique gene symbol only as a fallback.
qtl[, `:=`(
  SNP = variant_id,
  Chr = variant_chr,
  BP = variant_bp,
  Freq = af,
  Probe = as.character(phenotype_id),
  Probe_join = sub("\\..*$", "", as.character(phenotype_id)),
  b = slope,
  se = slope_se,
  p = pval_nominal
)]

id_index <- match(qtl$Probe_join, id_map$Probe_join)
qtl[, `:=`(
  Probe_Chr = id_map$Probe_Chr[id_index],
  Probe_bp = id_map$Probe_bp[id_index],
  Gene = id_map$Gene[id_index],
  Orientation = id_map$Orientation[id_index]
)]

symbol_index <- match(qtl$Probe_join, symbol_map$Probe_join)
use_symbol <- is.na(id_index) & !is.na(symbol_index)
qtl[use_symbol, `:=`(
  Probe_Chr = symbol_map$Probe_Chr[symbol_index[use_symbol]],
  Probe_bp = symbol_map$Probe_bp[symbol_index[use_symbol]],
  Gene = symbol_map$Gene[symbol_index[use_symbol]],
  Orientation = symbol_map$Orientation[symbol_index[use_symbol]]
)]

# Remove a complete probe if its chromosome or TSS annotation is invalid.
bad_probes <- unique(qtl[
  is.na(Probe_Chr) | is.na(Probe_bp) | Probe_bp <= 0L,
  Probe
])
if (length(bad_probes) > 0) {
  qtl <- qtl[!(Probe %chin% bad_probes)]
}

qtl[is.na(Gene) | trimws(Gene) == "", Gene := "."]
qtl[is.na(Orientation) | trimws(Orientation) == "", Orientation := "."]

qfile <- qtl[, .(
  SNP,
  Chr,
  BP,
  A1 = as.character(A1),
  A2 = as.character(A2),
  Freq,
  Probe,
  Probe_Chr,
  Probe_bp,
  Gene,
  Orientation,
  b,
  se,
  p
)]
setcolorder(qfile, qfile_columns)

dir.create(dirname(args$output), recursive = TRUE, showWarnings = FALSE)
fwrite(qfile, args$output, sep = "\t", quote = FALSE, na = ".")
message("Saved ", format(nrow(qfile), big.mark = ","), " rows to ", args$output)
