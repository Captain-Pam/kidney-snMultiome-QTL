#!/usr/bin/env python
"""Build the scored caQTL table: one row per variant x cell type.

Carries the measured RASQUAL effect alongside both models' predicted effects,
so a reader can reproduce the caQTL concordance analyses without rerunning any
model.

| Column | Meaning |
|---|---|
| `variant_id`, `chrom`, `pos`, `ref`, `alt`, `is_snp` | variant identity |
| `peak`, `peak_start`, `peak_end`, `dist_to_peak_center` | the RASQUAL peak |
| `celltype`, `celltype_label` | cell type, internal and display naming |
| `rasqual_chisq`, `rasqual_log10_bh_q` | association strength |
| `rasqual_effect_size`, `rasqual_effect_centered` | measured allelic effect |
| `chrombpnet_logfc`, `chrombpnet_abs_logfc`, `chrombpnet_logfc_q` | ChromBPNet |
| `borzoi_atac_logSUM` | Cerberus ATAC logSUM, from the fold ensemble |

Usage:
    python make_caqtl_scored.py [--out <path>]
"""
import argparse
import os

import pandas as pd

from caqtl_data import CELL_TYPES, load_merged
from figlib import CT_LABELS, cfg


def build(merged):
    records = []
    for celltype in CELL_TYPES:
        subset = merged[merged["cell_type"] == celltype].copy()
        if not len(subset):
            print(f"  warning: no rows for {celltype}")
            continue

        # Feature_ID encodes the peak as chr-start-end.
        peak = subset["Feature_ID"].str.rsplit("-", n=2, expand=True)
        peak_start = peak[1].astype(int)
        peak_end = peak[2].astype(int)

        record = pd.DataFrame(
            {
                "variant_id": (
                    subset["Chromosome"] + ":" + subset["SNP_position"].astype(str)
                    + ":" + subset["Ref_allele"] + ":" + subset["Alt_allele"]
                ),
                "chrom": subset["Chromosome"],
                "pos": subset["SNP_position"],
                "ref": subset["Ref_allele"],
                "alt": subset["Alt_allele"],
                "is_snp": (
                    (subset["Ref_allele"].str.len() == 1)
                    & (subset["Alt_allele"].str.len() == 1)
                ),
                "peak": subset["Feature_ID"],
                "peak_start": peak_start,
                "peak_end": peak_end,
                "celltype": celltype,
                "celltype_label": CT_LABELS[celltype],
                "dist_to_peak_center": (
                    subset["SNP_position"] - (peak_start + peak_end) // 2
                ),
                "rasqual_chisq": subset["Chi-square_statistic"],
                "rasqual_log10_bh_q": subset["Log_10_BH_Qvalue"],
                "rasqual_effect_size": subset["Effect_size"],
                # RASQUAL effect size lives on [0, 1] with no imbalance at 0.5;
                # centring gives a signed quantity comparable to an eQTL beta.
                "rasqual_effect_centered": subset["Effect_size"].astype(float) - 0.5,
                "chrombpnet_logfc": subset[f"logfc_{celltype}"],
                "chrombpnet_abs_logfc": subset[f"abs_logfc_{celltype}"],
                "chrombpnet_logfc_q": subset[f"logfc_q_{celltype}"],
                "borzoi_atac_logSUM": subset[f"logSUM_{celltype}"],
            }
        )

        # rs_ID in the RASQUAL tables is itself chrom:pos:ref:alt. Check rather
        # than carrying a duplicate column: a mismatch means the join key is
        # wrong and every model score in this block belongs to another variant.
        mismatched = (record["variant_id"].values != subset["rs_ID"].values).sum()
        if mismatched:
            raise ValueError(
                f"{celltype}: {mismatched} rows where the constructed variant_id "
                "disagrees with rs_ID"
            )

        records.append(record)
        print(
            f"  {celltype:12s} {len(record):7,d} rows, "
            f"{record['peak'].nunique():6,d} peaks"
        )

    return pd.concat(records, ignore_index=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default=None)
    args = parser.parse_args()

    out_path = args.out or os.path.join(cfg("FIGURE_DIR"), "caqtl_scored.csv")
    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)

    print("Loading caQTL scores...")
    table = build(load_merged())

    table.to_csv(out_path, index=False)
    print(
        f"\n{len(table):,} rows, {table['variant_id'].nunique():,} variants, "
        f"{table['peak'].nunique():,} peaks -> {out_path}"
    )


if __name__ == "__main__":
    main()
