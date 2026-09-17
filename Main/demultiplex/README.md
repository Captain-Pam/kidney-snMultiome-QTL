# Genotype demultiplexing

Each snMultiome library in this study was a **pool of two donors** (`HKxxxx_HKyyyy`).
This directory contains the code used to assign every RNA and ATAC barcode back
to its donor of origin (or flag it as a doublet), using **WGS genotypes** of the
pooled donors.

Robustness comes from running **five callers** — three on the RNA BAM and two on
the ATAC BAM — and taking a per-cell majority vote:

| Caller | Modality | Method |
|--------|----------|--------|
| `demuxlet` | RNA | genotype-based assignment ([popscle](https://github.com/statgen/popscle)) |
| `souporcell` (`souporcell_gt`) | RNA | genotype-aware clustering ([souporcell](https://github.com/wheaton5/souporcell)) |
| `vireo` (`vireo_gt`) | RNA | cellSNP-lite pileup + [vireo](https://github.com/single-cell-genetics/vireo) |
| `souporcell` (`souporcell_atac`) | ATAC | souporcell on the ATAC BAM |
| `vireo` (`vireo_atac`) | ATAC | cellSNP-lite pileup + vireo on the ATAC BAM |

All callers are seeded with the **known donor genotypes** (`k = 2`), so the
returned clusters/labels are already donor-named rather than anonymous.

## Layout

```
demuxlet/          demuxlet on the RNA BAM (+ VCF reordering helper)
souporcell_gt/     souporcell on the RNA BAM
souporcell_atac/   souporcell on the ATAC BAM
vireo_gt/          cellSNP-lite + vireo on the RNA BAM (+ VCF chr-renaming)
vireo_atac/        cellSNP-lite + vireo on the ATAC BAM
final_assignment/  merge the five callers into a consensus per-barcode call
```

Each caller folder holds an `example_command.sh` — the concrete commands run
for a single library (`HK2385_HK2814`, donors `HK2385 + HK2814`). In practice
one such Slurm job is submitted per library; results land in
`outputs/<library>/` (not committed).

## Inputs

- **Cell Ranger ARC BAMs** — `gex_possorted_bam.bam` (RNA) and
  `atac_possorted_bam.bam` (ATAC) per library (see `../cellranger_arc`).
- **Cell barcodes** — CellBender-called RNA barcodes and ATAC barcodes.
- **Donor WGS VCF** — joint-called multi-sample WGS (`wgs_donors.vcf.gz`, hg38).

## Pipeline

### 1. `demuxlet/` (RNA)

`example_command.sh`, per library:
1. Subset the WGS VCF to the two pooled donors and filter uninformative /
   low-quality sites (`bcftools`; drop hom-ref/hom-alt-in-both, `DP < 5`,
   `AF < 0.1`).
2. Reorder the VCF contigs to match the BAM header with
   `sort_vcf_same_as_bam.sh` (required by popscle), then `tabix`.
3. Run `demuxlet` with `--field GT`, restricted to the called barcodes, →
   `results.best`.

`sort_vcf_same_as_bam.sh` is the upstream popscle helper
([aertslab/popscle_helper_tools](https://github.com/aertslab/popscle_helper_tools)).
The per-donor VCF (steps 1–2) is prepared once and reused by the souporcell
steps.

### 2. `souporcell_gt/` (RNA) and `souporcell_atac/` (ATAC)

`example_command.sh`, per library: unzip the donor VCF prepared by demuxlet and
run `souporcell_pipeline.py` against the RNA or ATAC BAM with
`--known_genotypes`, `--known_genotypes_sample_names`, `-k 2`. The ATAC run adds
`--no_umi True`. Output: `output/clusters.tsv`. Runs in the `souporcell` conda env.

### 3. `vireo_gt/` (RNA) and `vireo_atac/` (ATAC)

Preprocessing (`vireo_gt/preprocess_vcf/rename.sh`): cellSNP-lite/vireo expect
`1,2,...` contig names, so rename `chr1→1` in the WGS VCF with
`bcftools annotate --rename-chrs chr_rename.txt`.

`example_command.sh`, per library:
1. Subset the (renamed) VCF to the two donors (`vireo_gt`).
2. Pileup informative sites with `cellsnp-lite` over the called barcodes
   (RNA: `--minMAF 0.1 --minCOUNT 20`; ATAC: `--minMAF 0 --minCOUNT 2 --UMItag None`).
3. Restrict the donor VCF to the pileup sites and run `vireo -c ... -d ...` →
   `vireo_out/donor_ids.tsv`, `vireo_out/prob_singlet.tsv.gz`.

Runs in the `cellSNP` conda env. The ATAC run reuses the per-donor VCF produced
by the RNA run.

### 4. `final_assignment/`

- `combine.py` — harmonizes every caller into `Sing / Doub / Sing_LLK / Doub_LLK`
  columns (prefixes `DMXT_`, `SPrna_`, `VRrna_`, `SPatac_`, `VRatac_`), writes
  per-library `rna.csv` / `atac.csv` / `all.csv`, then takes the **per-cell
  majority vote** across the five callers → `final.csv`. ATAC-only barcodes fall
  back to the vireo-ATAC call. Finally concatenates all libraries into
  `all_genotype.csv` (`<barcode>-<library> → singlet, doublet`).
- `compare.py` — QC: pairwise caller-agreement box plots.

Two libraries (`HK2524_HK3043`, `HK3039_HK3007`) contained only a single donor
and are assigned that donor directly.

## Sample fixes

Two pools had a mislabeled donor, applied in every caller command and in
`combine.py`:

- `HK2518_HK2763` → donors `HK2518, HK2762`
- `HK3106_HK2598` → donors `HK3106, HK2648`

## Software

- popscle / `demuxlet`
- souporcell (conda env `souporcell`)
- cellSNP-lite + vireo (conda env `cellSNP`)
- bcftools, samtools, tabix
- Slurm; Python 3 with pandas / numpy / matplotlib

## Notes on paths

Real input/output locations have been replaced with `path/to/...` placeholders.
Edit them (BAMs, barcodes, donor VCF, genome FASTA, per-caller output dirs) for
your environment before running.
