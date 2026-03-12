#!/bin/bash
#SBATCH --job-name=fc_prepare
#SBATCH --account=qsong1
#SBATCH --output=/orange/qsong1/Yansheng/logs/fc_prepare_%j.out
#SBATCH --error=/orange/qsong1/Yansheng/logs/fc_prepare_%j.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --mail-user=yanshengluo@ufl.edu
#SBATCH --mail-type=BEGIN,END,FAIL

set -euo pipefail

if [[ $# -lt 1 || $# -gt 2 ]]; then
  echo "USAGE: sbatch 02a_prepare_fc_batches.slurm <FASTQ_DIR> [BATCH_SIZE]" >&2
  exit 1
fi

DATA_DIR=$(realpath "$1")
DATASET=$(basename "${DATA_DIR}")
BATCH_SIZE="${2:-16}"

PROJECT_ROOT="/orange/qsong1/Yansheng"
OUT_ROOT="${PROJECT_ROOT}/03_counts/${DATASET}"
MANIFEST="${OUT_ROOT}/manifest/fastq_manifest.tsv"
STAR_ROOT="${OUT_ROOT}/star"
FC_ROOT="${OUT_ROOT}/featurecounts"
META_DIR="${OUT_ROOT}/metadata"

mkdir -p "${FC_ROOT}/batches" "${META_DIR}"

module load samtools

[[ -s "${MANIFEST}" ]] || { echo "ERROR: missing manifest: ${MANIFEST}" >&2; exit 1; }

BAM_LIST="${FC_ROOT}/bam_list.txt"
: > "${BAM_LIST}"

pe=0
se=0

while IFS=$'\t' read -r sample fq1 fq2; do
  [[ "${sample}" == "sample_id" ]] && continue

  BAM="${STAR_ROOT}/${sample}/Aligned.sortedByCoord.out.bam"
  [[ -s "${BAM}" ]] || { echo "ERROR: missing BAM for ${sample}: ${BAM}" >&2; exit 1; }
  samtools quickcheck -v "${BAM}" || { echo "ERROR: BAM fails quickcheck: ${BAM}" >&2; exit 1; }

  echo "${BAM}" >> "${BAM_LIST}"

  npaired=$(samtools view -c -f 1 -@ 1 "${BAM}" 2>/dev/null || echo 0)
  if (( npaired > 0 )); then
    pe=$((pe+1))
  else
    se=$((se+1))
  fi
done < "${MANIFEST}"

if (( pe > 0 && se > 0 )); then
  echo "ERROR: mixed paired-end and single-end BAMs detected" >&2
  exit 1
fi

if (( pe > 0 )); then
  LIB_LAYOUT="PAIRED_END"
else
  LIB_LAYOUT="SINGLE_END"
fi

echo "${LIB_LAYOUT}" > "${META_DIR}/library_layout.simple.txt"

rm -f "${FC_ROOT}"/batches/batch_*.txt
split -d -l "${BATCH_SIZE}" "${BAM_LIST}" "${FC_ROOT}/batches/batch_"

for f in "${FC_ROOT}"/batches/batch_*; do
  mv "$f" "${f}.txt"
done

NBATCH=$(find "${FC_ROOT}/batches" -maxdepth 1 -name 'batch_*.txt' | wc -l | awk '{print $1}')
echo "${NBATCH}" > "${META_DIR}/num_batches.txt"

echo "Dataset       : ${DATASET}"
echo "BAM list      : ${BAM_LIST}"
echo "Library layout: ${LIB_LAYOUT}"
echo "Batch size    : ${BATCH_SIZE}"
echo "Num batches   : ${NBATCH}"