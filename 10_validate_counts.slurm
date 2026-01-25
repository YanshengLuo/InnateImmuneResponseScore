#!/bin/bash
#SBATCH --job-name=validate_counts
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G

###############################################################################
# STEP 1 — GENE EXPRESSION TABLE VALIDATION & CLEANING
#
# PURPOSE
#   Validate and clean a gene-level RNA-seq count matrix BEFORE any normalization
#   or differential expression analysis.
#
# USAGE
#   sbatch 01_validate_counts.slurm <COUNTS_FILE> [OUTDIR] [SEP]
#
# ARGUMENTS
#   COUNTS_FILE : path to gene count matrix (required)
#   OUTDIR      : output directory (optional)
#                 default:
#                 /orange/qsong1/qsong1/qsong1/Yansheng/03_counts/GSE264344/featurecounts/validation
#   SEP         : delimiter (optional)
#                 default: tab
#
# OUTPUTS
#   - gene_counts_clean.tsv
#   - library_sizes.tsv
#   - qc_report.txt
#
# NEXT STEP
#   If PASS → proceed to Step 2 (design matrix alignment)
###############################################################################

set -euo pipefail

# =========================
# INPUT ARGUMENTS
# =========================
if [[ $# -lt 1 ]]; then
  echo "USAGE: sbatch 01_validate_counts.slurm <COUNTS_FILE> [OUTDIR] [SEP]"
  exit 1
fi

COUNTS="$1"

DEFAULT_OUTDIR="/orange/qsong1/qsong1/qsong1/Yansheng/03_counts/GSE264344/featurecounts/validation"
OUTDIR="${2:-$DEFAULT_OUTDIR}"

SEP="${3:-$'\t'}"

# =========================
# SETUP
# =========================
mkdir -p logs
mkdir -p "${OUTDIR}"

echo "--------------------------------------------------"
echo "Step 1: Gene expression table validation"
echo "Counts file : ${COUNTS}"
echo "Output dir  : ${OUTDIR}"
echo "Delimiter   : ${SEP}"
echo "--------------------------------------------------"

# Load Python if needed
# module load python/3.10

python3 validate_counts.py \
  --counts "${COUNTS}" \
  --outdir "${OUTDIR}" \
  --sep "${SEP}" \
  --drop-non-gene-rows \
  --collapse-duplicates \
  --outlier-factor 5
