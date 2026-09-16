#!/bin/bash
# Score caQTL variants with every Cerberus fold.
#
# Submits an array over (fold x variant chunk). The statistic is logSUM over a
# 1,024 bp window centred on the variant, taken from the cell-type-matched ATAC
# track. Restricting to a local window is what makes the score comparable to
# ChromBPNet's, which only sees 2,114 bp of context.
#
#   conda activate baskerville && bash 2_score_folds.sh

set -euo pipefail

source "$(dirname "$0")/../../../config.sh"

VCF="${WORK_DIR}/caqtl/caqtl_full_fdr0.1_inPeaks.vcf"
CHUNK_SIZE=1000

N_VARIANTS=$(grep -vc "^#" "${VCF}")
N_CHUNKS=$(( (N_VARIANTS + CHUNK_SIZE - 1) / CHUNK_SIZE ))
N_TASKS=$(( CERBERUS_N_FOLDS * N_CHUNKS ))

echo "${N_VARIANTS} variants -> ${N_CHUNKS} chunks x ${CERBERUS_N_FOLDS} folds = ${N_TASKS} tasks"

OUT_BASE="${WORK_DIR}/caqtl/scores/logSUM"
mkdir -p "${OUT_BASE}"

sbatch --job-name=caqtl_score \
       --partition="${SLURM_GPU_PARTITION}" \
       --gres="${SLURM_GPU_GRES}" \
       --mem=24G \
       --time=12:00:00 \
       --array="0-$(( N_TASKS - 1 ))%64" \
       --output="${OUT_BASE}/job_%A_%a.out" \
       --error="${OUT_BASE}/job_%A_%a.err" \
       --export=ALL,VCF="${VCF}",N_CHUNKS="${N_CHUNKS}",CHUNK_SIZE="${CHUNK_SIZE}",OUT_BASE="${OUT_BASE}" \
       "$(cd "$(dirname "$0")" && pwd)/_score_one.sh"
