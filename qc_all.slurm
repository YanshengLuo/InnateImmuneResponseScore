#!/bin/bash
#SBATCH --job-name=bulkRNA_qc
#SBATCH --output=/orange/qsong1/Yansheng/logs/qc_%j.out
#SBATCH --error=/orange/qsong1/Yansheng/logs/qc_%j.err
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --mail-user=yanshengluo@ufl.edu
#SBATCH --mail-type=BEGIN,END,FAIL

set -euo pipefail

# =========================
# Argument check
# =========================
if [[ $# -ne 1 ]]; then
  echo "USAGE: sbatch qc_all.slurm <FASTQ_DIR>" >&2
  exit 1
fi

DATA_DIR=$(realpath "$1")

if [[ ! -d "${DATA_DIR}" ]]; then
  echo "ERROR: FASTQ directory not found: ${DATA_DIR}" >&2
  exit 1
fi

# =========================
# Dataset naming
# =========================
DATASET=$(basename "${DATA_DIR}")

QC_ROOT="/orange/qsong1/Yansheng/02_qc/${DATASET}"
FASTQC_OUT="${QC_ROOT}/fastqc"
MULTIQC_OUT="${QC_ROOT}/multiqc"
SUMMARY_DIR="${QC_ROOT}/summary"
LOG_DIR="/orange/qsong1/Yansheng/logs"

# =========================
# Ensure all directories exist
# =========================
for DIR in "${QC_ROOT}" "${FASTQC_OUT}" "${MULTIQC_OUT}" "${SUMMARY_DIR}" "${LOG_DIR}"; do
  if [[ ! -d "${DIR}" ]]; then
    mkdir -p "${DIR}" || {
      echo "ERROR: Failed to create directory ${DIR}" >&2
      exit 1
    }
  fi
done

echo "=== QC job started ==="
date
echo "Dataset       : ${DATASET}"
echo "Input FASTQs  : ${DATA_DIR}"
echo "QC output dir : ${QC_ROOT}"
echo

# =========================
# Environment (HPC modules)
# =========================
module load fastqc
module load multiqc

# =========================
# Sanity checks
# =========================
echo "=== Sanity checks ==="

mapfile -t FQS < <(find "${DATA_DIR}" -maxdepth 1 -type f -name "*.fastq.gz" | sort)
N=${#FQS[@]}

if [[ "${N}" -eq 0 ]]; then
  echo "ERROR: No .fastq.gz files found in ${DATA_DIR}" >&2
  exit 1
fi

echo "FASTQ file count: ${N}"

# Paired-end naming heuristic
N_R1=$(find "${DATA_DIR}" -maxdepth 1 -type f \( -name "*_1.fastq.gz" -o -name "*_R1*.fastq.gz" \) | wc -l || true)
N_R2=$(find "${DATA_DIR}" -maxdepth 1 -type f \( -name "*_2.fastq.gz" -o -name "*_R2*.fastq.gz" \) | wc -l || true)

if [[ "${N_R1}" -gt 0 || "${N_R2}" -gt 0 ]]; then
  echo "Paired-end naming detected (heuristic): R1=${N_R1}, R2=${N_R2}"
else
  echo "No paired-end naming detected → likely single-end"
fi

# Gzip integrity check (sampled)
CHECK_K=$(( N < 5 ? N : 5 ))
echo "Checking gzip integrity (${CHECK_K}/${N})..."
for i in $(seq 0 $((CHECK_K-1))); do
  gzip -t "${FQS[$i]}"
done
echo "Gzip integrity check passed"

# Read length inference
READLEN=$(zcat "${FQS[0]}" | awk 'NR==2{print length($0); exit}')
HEADER=$(zcat "${FQS[0]}" | awk 'NR==1{print $0; exit}')

echo "Example FASTQ header:"
echo "${HEADER}"
echo "Inferred read length: ${READLEN}"

echo -e "file\tread_length\n${FQS[0]}\t${READLEN}" \
  > "${SUMMARY_DIR}/read_length_first_file.tsv"

# =========================
# FastQC
# =========================
echo
echo "=== Running FastQC ==="

PARALLEL=$(( SLURM_CPUS_PER_TASK > 1 ? SLURM_CPUS_PER_TASK : 1 ))

printf "%s\n" "${FQS[@]}" > "${SUMMARY_DIR}/fastq_list.txt"

cat "${SUMMARY_DIR}/fastq_list.txt" | \
  xargs -I {} -P "${PARALLEL}" \
  fastqc --quiet -t 1 -o "${FASTQC_OUT}" {}

echo "FastQC finished"

# =========================
# MultiQC
# =========================
echo
echo "=== Running MultiQC ==="

multiqc "${FASTQC_OUT}" -o "${MULTIQC_OUT}"

if [[ ! -f "${MULTIQC_OUT}/multiqc_report.html" ]]; then
  echo "ERROR: MultiQC report not generated" >&2
  exit 1
fi

echo "MultiQC report:"
echo "${MULTIQC_OUT}/multiqc_report.html"

# =========================
# Final summary
# =========================
echo
echo "=== QC completed successfully ==="
echo "Dataset      : ${DATASET}"
echo "FASTQ files  : ${N}"
echo "Read length  : ${READLEN}"
echo "QC directory : ${QC_ROOT}"
date
