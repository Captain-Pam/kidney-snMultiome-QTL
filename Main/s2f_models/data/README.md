# Pseudobulk coverage and model targets

Turns aligned multiome reads into the per-cell-type coverage both models train
on, and builds the targets table that tells the Cerberus model which tracks to
predict.

**The per-cell-type coverage files will be released with the paper.** Steps 1
and 2 document how they were produced; you will not normally rerun them. Step 3
is the one you do run, because the targets file has to point at wherever you put
the downloaded coverage.

The same coverage feeds both models: the merged **fragment** files are the input
to ChromBPNet peak calling and training, and the **`.w5`** coverage is what the
Cerberus targets file references.

## What each script does

| Script | In | Out |
|--------|----|-----|
| `1_split_sample_bam.sh` | one sample's Cell Ranger ARC GEX BAM and ATAC fragments, plus per-cell-type barcode lists | per-cell-type BAM and fragment file for that sample |
| `2_make_coverage.sh` | those per-sample outputs for one cell type | merged fragments (`.tsv.gz` + index) and coverage (`.bw` + `.w5`) for ATAC and both RNA strands |
| `3_make_targets.py` | `manifest.tsv` mapping each dataset to its `.w5` directory | `kidney_targets_w5_{human,mouse}.txt`, optionally with `_local` variants |
| `subset_fragments.py` | fragment file, barcode list | barcode-subset fragment file (helper for step 1) |

Run step 2 once per cell type for the pseudobulk tracks used in training, and
once per cell type × donor for the per-donor tracks the genotype-stratified
coverage figure needs.

## Things that will bite you

**RNA strand naming.** deeptools labels strands by read orientation, which for
this library chemistry is the reverse of the transcript strand:
`--filterRNAstrand forward` produces the **`rna-`** track and `reverse`
produces **`rna+`**. Output filenames follow the transcript strand, which is
what the targets file and every downstream track index assume.

**The `xf:i:25` filter** in step 1 keeps only confidently mapped, deduplicated
transcriptomic reads. Without it the RNA coverage is dominated by PCR
duplicates.

**`--local` in step 3** writes a second targets file carrying a `window` column.
That column is what makes `hound_snp --local_window` actually restrict its
statistic to a window around the variant; passing the flag with the plain
targets file is **silently ignored** and produces whole-sequence scores that are
not comparable to the published ones. See `../cerberus/README.md`.

**`description` format.** The targets `description` column is
`<ATAC|RNA+|RNA->:<cell_type>:<dataset_label>`. The evaluation and figure code
parses those three fields to group tracks, so the format matters.

## Datasets

The released model has a human head and a mouse head, each built from its own
manifest:

- **human** — Susztak CKD multiome (this study); Ledru control multiome
  (PRJNA909297); Muto PKD paired (PRJNA771425); Muto AKI paired (PRJNA1057812).
- **mouse** — White CKD multiome (Calico data, to be released); Chen time-series
  multiome (PRJNA1062423); Muto PKD multiome (PRJNA1117459); Muto IRI snATAC
  (PRJNA1057807); Kirita IRI snRNA (PRJNA578393).

These scripts are the ones used for the human data from this study. The mouse
atlas and the public BioProjects were processed to the same `.w5` form from
their own sources.

## Software

`subset-bam`, samtools, bgzip/tabix, `scatac_fragment_tools`, deeptools
`bamCoverage`, `bw_w5.py`, and Python with pandas.
