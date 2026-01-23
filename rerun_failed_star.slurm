#!/bin/bash
set -euo pipefail

# ------------------------------------------------------------
# rerun_failed_star.sh
#
# Purpose:
#   Rerun STAR alignment only for samples that failed in the first run.
#
# Requirements:
#   - Manifest exists:
#       03_counts/<DATASET>/manifest/fastq_manifest.tsv
#   - Failure list exists:
#       03_counts/<DATASET>/qc_alignment/star_failed_samples.tsv
#     with header: sample  task_id  reason
#   - STAR array script exists (the one you run for alignment):
#       01_star_align_array.slurm
#
# Usage:
#   bash rerun_failed_star.sh <FASTQ_DIR> <STAR_ARRAY_SCRIPT> [MAX_CONCURRENT]
#
# Example:
#   bash rerun_failed_star.sh /orange/qsong1/Yansheng/01_raw/fastq/GSE264344 \
#       /orange/qsong1/Yansheng/Hypergator_scripts/01_star_align_array.slurm 1
# ------------------------------------------------------------

if [[ $# -lt 2 || $# -gt 3 ]]; then
  echo "USAGE: bash rerun_failed_star.sh <FASTQ_DIR> <STAR_ARRAY_SCRIPT> [MAX_CONCURRENT]" >&2
  exit 1
fi

FASTQ_DIR=$(realpath "$1")
STAR_SCRIPT=$(realpath "$2")
MAX_CONCURRENT="${3:-1}"   # default: rerun one at a time (polite + safer)

DATASET=$(basename "$FASTQ_DIR")

# Adjust if your project root differs; this matches your current layout
PROJECT_ROOT="/orange/qsong1/Yansheng"
OUT_ROOT="${PROJECT_ROOT}/03_counts/${DATASET}"
MANIFEST="${OUT_ROOT}/manifest/fastq_manifest.tsv"
FAILED="${OUT_ROOT}/qc_alignment/star_failed_samples.tsv"
RERUN_DIR="${OUT_ROOT}/qc_alignment/rerun"
RERUN_TASKS="${RERUN_DIR}/rerun_task_ids.txt"
RERUN_SAMPLES="${RERUN_DIR}/rerun_samples.txt"

mkdir -p "${RERUN_DIR}"

[[ -s "${STAR_SCRIPT}" ]] || { echo "ERROR: STAR script not found: ${STAR_SCRIPT}" >&2; exit 1; }
[[ -s "${MANIFEST}" ]] || { echo "ERROR: Manifest not found: ${MANIFEST}" >&2; exit 1; }
[[ -s "${FAILED}" ]] || { echo "ERROR: Failed list not found (no failures or file missing): ${FAILED}" >&2; exit 1; }

# Extract failed sample IDs (skip header)
awk -F'\t' 'NR>1 && $1!="" {print $1}' "${FAILED}" | sort -u > "${RERUN_SAMPLES}"

if [[ ! -s "${RERUN_SAMPLES}" ]]; then
  echo "No failed samples found in ${FAILED}. Nothing to rerun."
  exit 0
fi

# Map sample IDs to task IDs using manifest line numbers:
# manifest format:
#   NR==1 header
#   NR==2 corresponds to task 0
#   NR==3 corresponds to task 1
# so task_id = NR - 2
awk -F'\t' '
  NR==FNR { bad[$1]=1; next }
  NR>1 && ($1 in bad) { print NR-2 }
' "${RERUN_SAMPLES}" "${MANIFEST}" | sort -n > "${RERUN_TASKS}"

if [[ ! -s "${RERUN_TASKS}" ]]; then
  echo "ERROR: Could not map failed samples to task IDs. Check naming consistency." >&2
  echo "Failed samples file: ${RERUN_SAMPLES}" >&2
  echo "Manifest: ${MANIFEST}" >&2
  exit 1
fi

# Create comma-separated list for sbatch --array
ARRAY_SPEC=$(paste -sd, "${RERUN_TASKS}")

echo "Dataset          : ${DATASET}"
echo "FASTQ_DIR         : ${FASTQ_DIR}"
echo "STAR array script : ${STAR_SCRIPT}"
echo "Failed samples    : $(wc -l < "${RERUN_SAMPLES}")"
echo "Task IDs to rerun  : $(wc -l < "${RERUN_TASKS}")"
echo "Array spec        : ${ARRAY_SPEC}"
echo "Max concurrent     : ${MAX_CONCURRENT}"
echo

echo "Submitting rerun array..."
jid=$(sbatch --array="${ARRAY_SPEC}%${MAX_CONCURRENT}" "${STAR_SCRIPT}" "${FASTQ_DIR}" | awk '{print $4}')
echo "Submitted rerun jobid: ${jid}"

echo
echo "Saved:"
echo "  Samples: ${RERUN_SAMPLES}"
echo "  TaskIDs : ${RERUN_TASKS}"
echo
echo "Monitor:"
echo "  squeue -u \$USER"
echo "  sacct -j ${jid} --format=JobID,State,ExitCode,Elapsed,MaxRSS -P"
