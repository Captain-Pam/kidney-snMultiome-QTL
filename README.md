# Kidney single-nucleus multiome project

Most kidney disease-associated genetic variants lie in non-coding regions, yet
linking them to regulatory elements, effector genes and relevant cell types
remains challenging. We integrated paired single-nucleus gene expression and
chromatin accessibility profiles from 191 human kidney samples with
whole-genome sequencing data from 94 donors. This study establishes a multiomic
resource spanning 12 kidney cell types and states. Cell-type-resolved molecular
QTL mapping identified regulatory variants associated with chromatin
accessibility (caQTLs) and gene expression (eQTLs). These results provide a
foundation for translating complex genetic associations into testable,
cell-type-specific mechanistic hypotheses for kidney disease.

This repository provides analysis code supporting data processing, quality
control and cell-type-resolved molecular QTL mapping. It also includes workflows
for statistical fine-mapping, characterization of molecular QTL genetic
architecture, sequence-to-function models, GWAS–molecular QTL colocalization,
CASCADE-based regulatory cascade analysis and gene prioritization.

## Analysis code

- [`cellbender`](Main/cellbender): GEX matrix preparation and CellBender background removal.
- [`eqtl`](Main/eqtl): cell-type-resolved cis-eQTL mapping.
- [`susie`](Main/susie): GWAS and eQTL statistical fine-mapping.
- [`coloc`](Main/coloc): GWAS-eQTL colocalization.
