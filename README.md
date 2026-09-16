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

![Study overview of the kidney single-nucleus multiome project](Figures/Figure1_Study%20overview.png)

This repository provides analysis code supporting data processing, quality
control and cell-type-resolved molecular QTL mapping. It also includes workflows
for statistical fine-mapping, characterization of molecular QTL genetic
architecture, sequence-to-function models, GWAS–molecular QTL colocalization,
CASCADE classification and gene prioritization.

## Analysis code

- [`cellranger_arc`](Main/cellranger_arc): alignment and quantification of paired Gene Expression and Chromatin Accessibility reads.
- [`cellbender`](Main/cellbender): GEX matrix preparation and CellBender background removal.
- [`eqtl`](Main/eqtl): cell-type-resolved cis-eQTL mapping.
- [`susie`](Main/susie): GWAS and eQTL statistical fine-mapping.
- [`coloc`](Main/coloc): GWAS-eQTL colocalization.
- [`s2f_models`](Main/s2f_models): sequence-to-function models (ChromBPNet and a fine-tuned Cerberus) — training, variant effect scoring, attribution and figures.

## Data availability

The principal summary-level results from this study will be made publicly
available upon publication. These resources will include single-nucleus eQTL
and caQTL results, molecular QTL fine-mapping results, GWAS–molecular QTL
colocalization results, sequence-to-function model outputs and
gene-prioritization results for five major kidney-function GWAS loci.

Key results can be explored through the [interactive Kidney single-nucleus
multiome portal](https://35.255.159.98.sslip.io/kidney-snmultiome/). This
temporary URL will be replaced with a permanent domain upon publication.

## Citation

Until the associated manuscript becomes publicly available, please cite this
repository using its URL and the specific commit or release used. Once
available, the manuscript should be cited as the primary reference.

## Contact

For enquiries about this study, please contact the corresponding authors:
Katalin Susztak (ksusztak@pennmedicine.upenn.edu) and David Kelley
(drk@calicolabs.com).
