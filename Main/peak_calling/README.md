# Cell-type-resolved peak calling

This directory documents peak calling and peak-matrix construction for the
kidney single-nucleus multiome dataset. The workflow starts from a
quality-controlled and annotated SnapATAC2 `AnnDataSet`.

Peaks are called for the following 12 kidney cell types and states:
`CNT_CD_PC`, `DCT`, `DTL_ATL`, `EC`, `IC`, `Immune`, `PEC`, `PTS`, `Podocyte`,
`Stromal`, `TAL` and `injPT`.

## Workflow

### 1. Cell-type-resolved peak calling

`1_call_consensus_peaks.py` runs MACS3 through SnapATAC2 for all 12 final cell
types and states.

### 2. Consensus peaks and peak-matrix filtering

`2_build_and_filter_peak_matrix.py`:

1. merges the cell-type-specific MACS3 peaks using the hg38 chromosome sizes
   and a peak half-width of 250 bp;
2. retains peaks on chromosomes 1–22;
3. removes peaks overlapping the hg38 blacklist;
4. constructs a nucleus-by-peak count matrix using paired insertions; and
5. retains the union of peaks passing the cell-type-specific nucleus and
   sample coverage criteria.

Within each cell type or state, a peak is retained when both conditions are met:

- it is detected in at least 0.5% of nuclei, with a minimum of 50 nuclei; and
- it is detected in at least 5% of represented samples, with a minimum of five
  samples.

A sample supports a peak when at least one nucleus from that sample has a
non-zero count at the peak. The final peak set is the union of retained peaks
across the 12 cell types and states.

## Inputs

- A writable SnapATAC2 `AnnDataSet` containing fragment information and the final cell-type/state annotation.
- The sample identifier column in the dataset `obs`.
- The ENCODE hg38 blacklist BED file.

## Outputs

- Cell-type-specific MACS3 peak calls stored in the SnapATAC2 dataset.
- `consensus_peaks.hg38.autosomes.blacklist_filtered.bed`: the final consensus peak coordinates before coverage filtering.
- `peak_filter_summary.tsv`: cell-type-specific filtering thresholds and peak counts.
- `retained_peaks.tsv`: the final union of retained peaks.
- `peak_matrix.filtered.h5ad`: the filtered nucleus-by-peak count matrix.
