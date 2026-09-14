# SuSiE fine-mapping

Statistical fine-mapping of GWAS and single-cell eQTL summary statistics with
[SuSiE-RSS](https://stephenslab.github.io/susieR/) (`runsusie`). Fine-mapping is
run **independently** on each dataset; the fitted objects are consumed later by
`../coloc/coloc_susie` for colocalization.

Reference LD comes from the precomputed **UKBB European** LD blocks
(Broad / alkesgroup, genome build **hg19**). Because the eQTLs are called in
hg38, they are lifted over to hg19 so both datasets share the LD panel.

## Layout

```
ld_reference/      Build/query the UKBB hg19 LD panel
gwas_finemapping/  Fine-map the eGFR GWAS (per independent locus)
eqtl_finemapping/  Fine-map the sc-eQTLs (per cell type / gene)
```

## Paths

| Placeholder | Contents |
|-------------|----------|
| `data/eqtl/tensorqtl/all_parquets/<ct>.parquet` | tensorQTL eQTL sumstats per cell type (hg38) |
| `data/eqtl/expression_bed/<ct>.expression.bed.gz` | expression matrices (for `sdY`, `N`) |
| `data/gwas/eGFR_sumstats.harmonized.txt.gz` | harmonized eGFR meta-analysis sumstats (hg19) |
| `data/gwas/eGFR_loci.rds` | independent GWAS lead SNPs |
| `data/reference/ukbb_ld_scores/` | UKBB EUR LD blocks (`chr_start_end.gz` + `.npz`, hg19) |
| `data/reference/hg38ToHg19.over.chain` | liftOver chain |
| `../ld_reference/LD_gr.rds`, `query_LD.py` | produced/kept in `ld_reference/` |

## `ld_reference/`

| File | Inputs | Outputs |
|------|--------|---------|
| `process_LD.R` | UKBB LD block files `chr_start_end.gz` + `.npz` in `ld_scores/` | `LD_gr.rds` — a `GRanges` index of every LD block region (with paths to its `.gz`/`.npz`), used to pick the block covering a locus |
| `query_LD.py` | `--ld_prefix` (an LD block), `--input_file` (SNP table with `hg19_id1`/`hg19_id2` columns), `--output_file` | Dense LD (correlation) matrix CSV, subset to the query SNPs, `chr:pos:A1:A2` row/col IDs, multi-allelic sites dropped |

## `gwas_finemapping/` (eGFR GWAS)

| File | Inputs | Outputs |
|------|--------|---------|
| `process_egfr.R` | Harmonized eGFR meta-analysis sumstats (hg19); independent lead-SNP list (`eGFR_loci.rds`) | `egfr_gr.rds` (1 Mb window `GRanges` per lead locus), `loci_gr.rds`, `egfr_list_1m.rds` (per-locus sumstats, both allele orientations `vid`/`vid2`) |
| `run_susie.R` | `<locus>` arg; `egfr_list_1m.rds`, `egfr_gr.rds`, `LD_gr.rds`; calls `query_LD.py` | Per locus: `rss.rds` (fitted SuSiE), `res_susie.csv` (per-SNP PIP = `variable_prob`, credible set = `cs`), `meta.csv` (QC), `susie.pdf` |
| `make_submission.R` | `egfr_list_1m.rds` | One Slurm job per locus (`outputs/<locus>/submission.sh`), submits them, then aggregates all `meta.csv` → `susie_meta.csv` |

## `eqtl_finemapping/` (sc-eQTL)

| File | Inputs | Outputs |
|------|--------|---------|
| `run_susie.R` | `<cell_type>` arg; tensorQTL parquet (hg38, `egene_qval < 0.1`); expression BED (for `sdY`, `N`); hg38→hg19 chain; `LD_gr.rds`; calls `query_LD.py` | Per gene: `<ct>/<gene>/rss.rds`, `res_susie.csv`, `susie.pdf`; per cell type `susie_<ct>.csv` (QC table over all genes) |
| `gene_windows_hg19.R` | The `res_susie.csv` files above | `gene_gr_hg19.rds` — hg19 span (`GRanges`) of each fine-mapped gene, used for GWAS↔gene overlap in coloc |

## Per-locus pipeline (both GWAS and eQTL)

1. (eQTL only) **liftOver hg38→hg19**; drop SNPs lifting to duplicate positions.
2. Pick the LD block covering the most SNPs; query the LD matrix (`query_LD.py`).
3. **Harmonize effect alleles** to the LD panel (flip `BETA`/`slope` when the
   alternate orientation matches).
4. Detect and drop **outlier SNPs** via `kriging_rss` (|`z_std_diff`| > 5),
   see susieR issue #182.
5. QC: `estimate_s_rss` (`lambda`, sumstat–LD consistency) and `check_alignment`.
6. `runsusie(..., max_iter = 500)`; save fitted object + per-SNP table.

## Key `meta.csv` / `susie_<ct>.csv` columns

`total_snps`, `outlier_snps`, `finemap_snps`, `lambda`, `alignment`,
`converge` (0/1), `cs` (number of credible sets).
