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
trap 'echo "ERROR at line $LINENO. Exit code $?"; exit 1' ERR

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
  mkdir -p "${DIR}" || { echo "ERROR: Failed to create ${DIR}" >&2; exit 1; }
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

N_R1=$(find "${DATA_DIR}" -maxdepth 1 -type f \( -name "*_1.fastq.gz" -o -name "*_R1*.fastq.gz" \) | wc -l || true)
N_R2=$(find "${DATA_DIR}" -maxdepth 1 -type f \( -name "*_2.fastq.gz" -o -name "*_R2*.fastq.gz" \) | wc -l || true)

if [[ "${N_R1}" -gt 0 || "${N_R2}" -gt 0 ]]; then
  echo "Paired-end naming detected (heuristic): R1=${N_R1}, R2=${N_R2}"
else
  echo "No paired-end naming detected → likely single-end"
fi

CHECK_K=$(( N < 5 ? N : 5 ))
echo "Checking gzip integrity (${CHECK_K}/${N})..."
for i in $(seq 0 $((CHECK_K-1))); do
  gzip -t "${FQS[$i]}"
done
echo "Gzip integrity check passed"
echo

# Write FASTQ list early (always useful)
printf "%s\n" "${FQS[@]}" > "${SUMMARY_DIR}/fastq_list.txt"

# =========================
# Read length inference (SAFE: no SIGPIPE)
# =========================
echo "Inferring read length from first file..."
FIRST="${FQS[0]}"

python3 - <<PY
import gzip
f = "${FIRST}"
out = "${SUMMARY_DIR}/read_length_first_file.tsv"
try:
    with gzip.open(f, "rt") as fh:
        header = fh.readline().rstrip("\n")
        seq = fh.readline().rstrip("\n")
    readlen = len(seq) if seq else "NA"
    with open(out, "w") as w:
        w.write("file\tread_length\n")
        w.write(f"{f}\t{readlen}\n")
    print("Example FASTQ header:")
    print(header)
    print(f"Inferred read length: {readlen}")
except Exception as e:
    with open(out, "w") as w:
        w.write("file\tread_length\n")
        w.write(f"{f}\tNA\n")
    print("WARNING: Could not infer read length:", e)
PY

# =========================
# FastQC
# =========================
echo
echo "=== Running FastQC on all files ==="
PARALLEL=$(( SLURM_CPUS_PER_TASK > 1 ? SLURM_CPUS_PER_TASK : 1 ))

cat "${SUMMARY_DIR}/fastq_list.txt" | \
  xargs -I {} -P "${PARALLEL}" \
  fastqc --quiet -t 1 -o "${FASTQC_OUT}" {}

echo "FastQC finished"
echo

# =========================
# MultiQC
# =========================
echo "=== Running MultiQC ==="
multiqc "${FASTQC_OUT}" -o "${MULTIQC_OUT}"

if [[ ! -f "${MULTIQC_OUT}/multiqc_report.html" ]]; then
  echo "ERROR: MultiQC report not generated" >&2
  exit 1
fi

echo "MultiQC report:"
echo "${MULTIQC_OUT}/multiqc_report.html"
echo

echo "=== QC completed successfully ==="
date
