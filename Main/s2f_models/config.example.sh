# Configuration for the sequence-to-function pipelines.
#
#   cp config.example.sh config.sh    # then edit the paths below
#
# Shell scripts do:   source ../config.sh
# Python scripts do:  from figlib import load_config   (parses this same file)
#
# Every path below is an input you supply; nothing in this repository is a
# data file. See README.md for what each one is and where to obtain it.

# ─── Reference genome ─────────────────────────────────────────────────────────
# The two models were run against separate hg38 copies on our cluster. They are
# the same assembly; keep them pointed at one file unless you have a reason not
# to. CERBERUS_FASTA must be the copy the model was trained against.
export GENOME_FASTA=/path/to/hg38.fa
export GENOME_CHROM_SIZES=/path/to/hg38.chrom.sizes
export GENOME_BLACKLIST=/path/to/hg38.blacklist.bed.gz
export CERBERUS_FASTA="${GENOME_FASTA}"
export CHROMBPNET_FASTA="${GENOME_FASTA}"

export MM10_FASTA=/path/to/mm10.fa

# ─── Raw sequencing output (only needed to rerun data/, which you should not) ─
# One directory per sample, each holding Cell Ranger ARC's outs/.
export CELLRANGER_DIR=/path/to/cellranger_arc

# ─── Released coverage data (see data/README.md) ──────────────────────────────
# Per-cell-type pseudobulk fragments and coverage. These are released with the
# paper; data/ shows how they were produced but you do not need to rerun it.
export FRAG_DIR=/path/to/pseudobulk/bam_frag                 # <ct>_merged_fragments.tsv.gz
export FRAG_DIR_INDIVIDUAL=/path/to/pseudobulk/bam_frag_individual  # <ct>_<donor>_fragments.tsv.gz
export W5_DIR=/path/to/pseudobulk/bw_w5                      # <ct>_{atac,rna+,rna-}.w5

# ─── Cerberus training data and models ────────────────────────────────────────
export CERBERUS_DATA_HG38=/path/to/cerberus_data/hg38        # hound_data output
export CERBERUS_DATA_MM10=/path/to/cerberus_data/mm10
export CERBERUS_MODEL_DIR=/path/to/cerberus_folds            # f{0..7}c0/train/{params.json,model_best.pth}
export CERBERUS_N_FOLDS=8

# Pretrained Cerberus foundation checkpoints, one per fold, used to initialise
# fine-tuning: ${CERBERUS_PRETRAINED_DIR}/f{I}c0/train/model_best.pth
export CERBERUS_PRETRAINED_DIR=/path/to/cerberus_foundation

# Targets files written by data/3_make_targets.py.
export TARGETS_HUMAN="${PWD}/targets/kidney_targets_w5_human.txt"
export TARGETS_HUMAN_LOCAL="${PWD}/targets/kidney_targets_w5_human_local.txt"
export TARGETS_MOUSE="${PWD}/targets/kidney_targets_w5_mouse.txt"
# 11-track ATAC subset used for the example-locus ISM runs.
export TARGETS_ATAC_SUBSET="${PWD}/targets/targets_atac_subset.txt"

# ─── ChromBPNet models ────────────────────────────────────────────────────────
export CHROMBPNET_MODEL_DIR=/path/to/chrombpnet_models       # <ct>_fold_<k>/models/chrombpnet_nobias.h5
export CHROMBPNET_PEAK_DIR=/path/to/chrombpnet/peak_beds     # <ct>_peaks.bed, <ct>_nonpeaks_negatives.bed
export CHROMBPNET_MACS2_DIR=/path/to/chrombpnet/macs2_peaks  # <ct>_peaks.narrowPeak
export CHROMBPNET_N_FOLDS=5

# Chromosome-split definitions, one JSON per fold: ${SPLITS_DIR}/fold_<k>.json
export SPLITS_DIR=/path/to/chrombpnet/splits
# Pre-trained Tn5 bias models (ENCODE ENCSR291GJU, HepG2), one per fold:
#   ${BIAS_MODEL_DIR}/fold_<k>/model.bias.fold_<k>.ENCSR291GJU.h5
export BIAS_MODEL_DIR=/path/to/chrombpnet/bias_models
# Blacklist extended by 1057 bp either side to cover the model receptive field.
export BLACKLIST_EXT=/path/to/blacklist_ext1057.bed

# ─── External tools ───────────────────────────────────────────────────────────
# Local clone of the kundajelab variant-scorer repository.
export VARIANT_SCORER_SRC=/path/to/variant-scorer/src
# baskerville (PyTorch) and chrombpnet are expected on PATH via their conda
# environments; only their console entry points are used.
export CONDA_ENV_BASKERVILLE=baskerville
export CONDA_ENV_CHROMBPNET=chrombpnet

# ─── Annotation and motifs ────────────────────────────────────────────────────
export GENCODE_GTF=/path/to/gencode.v48.basic.annotation.gtf
# Subset of the above restricted to eQTL target genes, written by
# cerberus/score/eqtl/1_prepare_gtf.py.
export EQTL_GTF="${PWD}/annotation/eqtl_genes_gencode48.gtf"
export EQTL_GENE_ID_MAP="${PWD}/annotation/eqtl_gene_id_map.tsv"
# JASPAR2022 CORE vertebrates non-redundant, MEME format. Used for the motif
# alignment shown alongside attribution logos.
export MOTIF_DB=/path/to/JASPAR2022_CORE_vertebrates_non-redundant_v2.meme
# Database `modisco report` matches discovered motifs against. We used HOCOMOCO
# v12 core, as distributed with gReLU.
export MODISCO_REPORT_MEME=/path/to/H12CORE_meme_format.meme

# ─── Molecular QTL inputs (controlled access — see README.md) ─────────────────
export RASQUAL_CAQTL_DIR=/path/to/rasqual/caqtl              # <ct>_full_caPeaks_FDR0.1.tsv
export CAQTL_SUSIE_TSV=/path/to/susie/caqtl_merged_all_celltypes.tsv.gz
export EQTL_SUSIE_TSV=/path/to/susie/eqtl_merged_all_celltypes.tsv.gz
export GWAS_SUSIE_DIR=/path/to/susie/eGFR                    # one directory per locus
export GENOTYPE_VCF=/path/to/donor_genotypes.vcf.gz          # phased, indexed

# ─── Working directories ──────────────────────────────────────────────────────
export WORK_DIR=/path/to/work                                # scores, ISM, predictions
export FIGURE_DIR="${WORK_DIR}/figures"

# ─── Cluster ──────────────────────────────────────────────────────────────────
# Adjust to your scheduler, or ignore the sbatch headers and run the commands
# directly. Scripts read these when building their SLURM directives.
export SLURM_GPU_PARTITION=gpu
export SLURM_GPU_GRES=gpu:1
export SLURM_CPU_PARTITION=standard
