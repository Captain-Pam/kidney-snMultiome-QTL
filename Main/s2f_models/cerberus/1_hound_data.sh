#!/bin/bash
# Extract binned coverage into the per-fold training examples.
#
# Reads each target track's .w5 coverage over the sequence intervals defined by
# the cross-fold split, bins it, and writes the zarr examples the trainer
# consumes. Run once per genome.
#
# The split files (sequences.bed, contigs.bed, contig_components.txt, nets*.bed,
# mseqs_unmap.npy) are NOT generated here — see README.md. They come with the
# pretrained Cerberus dataset and must be placed in the output directory first,
# so that our folds line up with the ones the foundation model was trained on.
#
#   conda activate baskerville && bash 1_hound_data.sh hg38

set -euo pipefail

source "$(dirname "$0")/../config.sh"

if [[ $# -lt 1 ]]; then
    echo "usage: 1_hound_data.sh {hg38|mm10}" >&2
    exit 1
fi
GENOME=$1

case "${GENOME}" in
    hg38)
        FASTA="${CERBERUS_FASTA}"
        TARGETS="${TARGETS_HUMAN}"
        OUT_DIR="${CERBERUS_DATA_HG38}"
        ;;
    mm10)
        FASTA="${MM10_FASTA}"
        TARGETS="${TARGETS_MOUSE}"
        OUT_DIR="${CERBERUS_DATA_MM10}"
        ;;
    *)
        echo "unknown genome: ${GENOME}" >&2
        exit 1
        ;;
esac

if [[ ! -f "${OUT_DIR}/sequences.bed" ]]; then
    echo "Missing ${OUT_DIR}/sequences.bed - see README.md on obtaining the splits" >&2
    exit 1
fi

# -l 786432     input sequence length
# -w 32         output bin width -> 24576 bins per sequence
# -f 16         cross-validation folds
# -d 3          down-weight/limit on duplicated sequence
# -r 64         sequences per write job
# --umap_clip   clip coverage in unmappable regions to this fraction
# -z 2          zarr compression level
# --restart     resume a partially completed run instead of starting over
hound_data \
    --restart \
    -l 786432 \
    -w 32 \
    -f 16 \
    -d 3 \
    -r 64 \
    --umap_clip 0.5 \
    -z 2 \
    -p 64 \
    -b "${GENOME_BLACKLIST}" \
    -u "${OUT_DIR}/umap.bed" \
    -o "${OUT_DIR}" \
    "${FASTA}" \
    "${TARGETS}"

echo "Wrote ${OUT_DIR}/examples/fold{0..15}.zarr and statistics.json"
