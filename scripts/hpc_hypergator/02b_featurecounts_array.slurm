#!/bin/bash
#SBATCH --job-name=fc_batch
#SBATCH --account=qsong1
#SBATCH --output=/orange/qsong1/Yansheng/logs/fc_batch_%A_%a.out
#SBATCH --error=/orange/qsong1/Yansheng/logs/fc_batch_%A_%a.err
#SBATCH --time=16:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --mail-user=yanshengluo@ufl.edu
#SBATCH --mail-type=BEGIN,END,FAIL

set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "USAGE: sbatch --array=0-(N-1) 02b_featurecounts_array.slurm <DATASET>" >&2
  exit 1
fi

DATASET="$1"

PROJECT_ROOT="/orange/qsong1/Yansheng"
OUT_ROOT="${PROJECT_ROOT}/03_counts/${DATASET}"
FC_ROOT="${OUT_ROOT}/featurecounts"
META_DIR="${OUT_ROOT}/metadata"

DATASET_LOWER=$(echo "${DATASET}" | tr '[:upper:]' '[:lower:]')
REF_MOUSE="${PROJECT_ROOT}/00_metadata/reference/mm39_gencodeM38"
REF_HUMAN="${PROJECT_ROOT}/00_metadata/reference/GRCh38_gencode49"

if [[ "${DATASET_LOWER}" == *_human ]]; then
  REF_ROOT="${REF_HUMAN}"
else
  REF_ROOT="${REF_MOUSE}"
fi

GTF_GZ="${REF_ROOT}/genes.gtf.gz"
[[ -s "${GTF_GZ}" ]] || { echo "ERROR: missing GTF: ${GTF_GZ}" >&2; exit 1; }

BATCH_FILE=$(printf "%s/batches/batch_%02d.txt" "${FC_ROOT}" "${SLURM_ARRAY_TASK_ID}")
OUT_FILE=$(printf "%s/gene_counts.batch_%02d.tsv" "${FC_ROOT}" "${SLURM_ARRAY_TASK_ID}")

[[ -s "${BATCH_FILE}" ]] || { echo "ERROR: missing batch file: ${BATCH_FILE}" >&2; exit 1; }
[[ -s "${META_DIR}/library_layout.simple.txt" ]] || { echo "ERROR: missing layout file" >&2; exit 1; }

LIB_LAYOUT=$(cat "${META_DIR}/library_layout.simple.txt")

module load subread

GTF_TMP="${FC_ROOT}/genes.${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.gtf"
trap 'rm -f "${GTF_TMP}"' EXIT

zcat "${GTF_GZ}" > "${GTF_TMP}"
[[ -s "${GTF_TMP}" ]] || { echo "ERROR: failed to create temp GTF" >&2; exit 1; }

if ! awk '$3=="exon"{found=1; exit} END{exit(found?0:1)}' "${GTF_TMP}"; then
  echo "ERROR: Annotation has no exon features: ${GTF_TMP}" >&2
  exit 1
fi

mapfile -t BAMS < "${BATCH_FILE}"

FC_MODE_FLAGS=()
if [[ "${LIB_LAYOUT}" == "PAIRED_END" ]]; then
  FC_MODE_FLAGS=(-p --countReadPairs -B -C)
fi

featureCounts \
  -T "${SLURM_CPUS_PER_TASK}" \
  -a "${GTF_TMP}" \
  -o "${OUT_FILE}" \
  -g gene_id \
  -t exon \
  --primary \
  "${FC_MODE_FLAGS[@]}" \
  "${BAMS[@]}"

echo "DONE: ${OUT_FILE}"