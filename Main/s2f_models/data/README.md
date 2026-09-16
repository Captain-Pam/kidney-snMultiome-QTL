# Pseudobulk coverage and model targets

**The per-cell-type coverage files are released with the paper.** Steps 1 and 2
here document how they were produced from raw Cell Ranger ARC output; you do not
need to rerun them. Step 3 is the one you will run, because the targets file
must point at wherever you put the downloaded coverage.

The same coverage feeds both models:

- the merged **fragment** files are the input to ChromBPNet peak calling and
  training (`../chrombpnet`),
- the **`.w5`** coverage is what the Cerberus targets file points at
  (`../cerberus`).

## Layout

| Script | In | Out |
|--------|----|-----|
| `1_split_sample_bam.sh` | one sample's Cell Ranger ARC `gex_possorted_bam.bam` + `atac_fragments.tsv.gz`; per-cell-type barcode CSVs | `<sample>/split_bams/<ct>_filter_final.bam`, `<sample>/split_fragments/<ct>_atac_fragments.tsv.gz` |
| `2_make_coverage.sh` | the step-1 outputs for one cell type, across all samples | `bam_frag/<ct>_merged_fragments.tsv.gz` (+ `.tbi`), `bw_w5/<ct>_{atac,rna+,rna-}.{bw,w5}` |
| `3_make_targets.py` | `manifest.tsv` listing dataset label → `.w5` directory | `targets/kidney_targets_w5_{human,mouse}.txt` (+ `_local` variants) |
| `subset_fragments.py` | fragment file + barcode list | barcode-subset fragment file (helper for step 1) |

## 1. Splitting a sample

Cells are assigned to a cell type upstream; this stage only needs the resulting
barcode lists, one headerless CSV per cell type.

**RNA.** `subset-bam` selects the cell type's barcodes, then reads are filtered
to `xf:i:25` — Cell Ranger's flag for confidently mapped, deduplicated
transcriptomic reads. Skipping that filter leaves coverage dominated by PCR
duplicates.

**ATAC.** `subset_fragments.py` filters the fragment file on the barcode column.

## 2. Merging to pseudobulk

Fragments are concatenated across samples, coordinate-sorted, bgzipped and
tabix-indexed; the RNA BAMs are merged and sorted.

- **ATAC** → bigWig with `scatac_fragment_tools bigwig -x` (Tn5 insertion sites).
- **RNA** → two bigWigs with `bamCoverage --binSize 1`, split by strand.

  Note the strand naming. deeptools labels strands by read orientation, which
  for this library chemistry is the reverse of the transcript strand, so
  `--filterRNAstrand forward` produces the **`rna-`** track and
  `--filterRNAstrand reverse` produces **`rna+`**. The output filenames follow
  the transcript strand, which is what the model targets assume.

Both are then converted to `.w5` (an HDF5 coverage format read directly by
`hound_data`). `bw_w5.py` ships with the **TensorFlow** baskerville repository
(`calico/baskerville`) — a different package from the PyTorch `baskerville` used
for training and scoring. Both are needed.

Run this per cell type for the pseudobulk tracks used in training, and per
cell type × donor for the per-donor tracks used in genotype-stratified coverage
plots (`../figures/chr16_caqtl_locus.py`).

## 3. Targets files

`3_make_targets.py` walks each dataset's `.w5` directory and emits the
baskerville targets table:

| Column | Value |
|--------|-------|
| `identifier` | `<prefix>_<FILENAME>`, uppercased with underscores removed |
| `file` | absolute path to the `.w5` |
| `clip` | 200 (ATAC), 300 (RNA) |
| `scale` | 1 |
| `sum_stat` | `sum_sqrt` — coverage is square-root-summed within each bin |
| `strand_pair` | index of the opposite-strand track; ATAC points at itself |
| `description` | `<ATAC\|RNA+\|RNA->:<cell_type>:<dataset_label>` |

The `description` column is what the evaluation and figure code parses to group
tracks by assay, cell type and dataset, so its three-field format matters.

`--local` additionally writes a `_local` copy carrying a `window` column set to
`local`. That column is what makes `hound_snp --local_window` restrict its
statistic to a window around the variant — **the flag alone does nothing without
it**, and caQTL scores computed with the plain targets file are not comparable
to the published ones. See `../cerberus/score/caqtl/`.

The released model has two heads, built from two manifests:

- **human**, 182 tracks (66 ATAC, 58 RNA+, 58 RNA−): Susztak CKD multiome (this
  study), Ledru control multiome (PRJNA909297), Muto PKD paired (PRJNA771425),
  Muto AKI paired (PRJNA1057812).
- **mouse**, 373 tracks (127 ATAC, 123 RNA+, 123 RNA−): White CKD multiome (this
  study), Chen time-series multiome (PRJNA1062423), Muto PKD multiome
  (PRJNA1117459), Muto IRI snATAC (PRJNA1057807), Kirita IRI snRNA
  (PRJNA578393).

Only the two "this study" datasets were produced by the scripts here; the public
BioProjects were processed to the same `.w5` form from their own accessions.

## Software

`subset-bam` · samtools 1.18 · bgzip/tabix (htslib 1.18) ·
`scatac_fragment_tools` · deeptools `bamCoverage` · baskerville (TensorFlow, for
`bw_w5.py`) · Python 3.11 with pandas
