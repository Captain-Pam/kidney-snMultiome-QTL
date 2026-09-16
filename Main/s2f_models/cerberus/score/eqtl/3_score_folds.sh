#!/bin/bash
# Score fine-mapped eQTL variants with every Cerberus fold.
#
# The statistic is covgene/logSED: predicted RNA coverage summed over the
# annotated exons of the target gene, reported as the log ratio between alleles.
# Unlike the caQTL scoring there is no --local_window — the gene's exons define
# the region, so the score can pick up effects far from the variant.
#
#   conda activate baskerville && bash 3_score_folds.sh

set -euo pipefail

source "$(dirname "$0")/../../../config.sh"

VCF="${WORK_DIR}/eqtl/eqtl_susie_cs.vcf"
CHUNK_SIZE=1000

N_VARIANTS=$(grep -vc "^#" "${VCF}")
N_CHUNKS=$(( (N_VARIANTS + CHUNK_SIZE - 1) / CHUNK_SIZE ))
N_TASKS=$(( CERBERUS_N_FOLDS * N_CHUNKS ))

echo "${N_VARIANTS} variants -> ${N_CHUNKS} chunks x ${CERBERUS_N_FOLDS} folds = ${N_TASKS} tasks"

OUT_BASE="${WORK_DIR}/eqtl/scores/logSED"
mkdir -p "${OUT_BASE}"

sbatch --job-name=eqtl_score \
       --partition="${SLURM_GPU_PARTITION}" \
       --gres="${SLURM_GPU_GRES}" \
       --mem=24G \
       --time=12:00:00 \
       --array="0-$(( N_TASKS - 1 ))%64" \
       --output="${OUT_BASE}/job_%A_%a.out" \
       --error="${OUT_BASE}/job_%A_%a.err" \
       --export=ALL,VCF="${VCF}",N_CHUNKS="${N_CHUNKS}",CHUNK_SIZE="${CHUNK_SIZE}",OUT_BASE="${OUT_BASE}" \
       "$(cd "$(dirname "$0")" && pwd)/_score_one.sh"
