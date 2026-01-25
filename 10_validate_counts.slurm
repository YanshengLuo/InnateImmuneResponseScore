#!/bin/bash
#SBATCH --job-name=validate_counts
#SBATCH --output=/orange/qsong1/Yansheng/logs/%x_%j.out
#SBATCH --error=/orange/qsong1/Yansheng/logs/%x_%j.err
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
# AUTOMATIC FEATURECOUNTS HANDLING
#   If input is raw featureCounts output (starts with "# Program:featureCounts"),
#   this script will first generate a clean matrix:
#     - keeps only Geneid + count columns (drops Chr/Start/End/Strand/Length)
#     - renames BAM path columns to SRR IDs
#     - removes meta rows starting with "__"
#   Then it runs validate_counts.py on the cleaned matrix.
#
# USAGE
#   sbatch 01_validate_counts.slurm <COUNTS_FILE> [OUTDIR] [SEP]
#
# ARGUMENTS
#   COUNTS_FILE : path to gene count table (required)
#   OUTDIR      : output directory (optional)
#                 default:
#                 /orange/qsong1/qsong1/qsong1/Yansheng/03_counts/GSE264344/featurecounts/validation
#   SEP         : delimiter (optional)
#                 default: tab
#
# OUTPUTS (OUTDIR)
#   - gene_counts.matrix.tsv     (only if featureCounts preprocessing is triggered)
#   - gene_counts_clean.tsv
#   - library_sizes.tsv
#   - qc_report.txt
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
mkdir -p /orange/qsong1/Yansheng/logs
mkdir -p "${OUTDIR}"

echo "--------------------------------------------------"
echo "Step 1: Gene expression table validation"
echo "Counts file : ${COUNTS}"
echo "Output dir  : ${OUTDIR}"
echo "Delimiter   : ${SEP}"
echo "--------------------------------------------------"

# =========================
# CONDA ENV
# =========================
source /blue/wenjunxie/yanshengluo/miniconda3/etc/profile.d/conda.sh
conda activate imrs

# =========================
# AUTO-DETECT FEATURECOUNTS RAW OUTPUT
# =========================
INPUT_FOR_VALIDATION="${COUNTS}"
MATRIX_OUT="${OUTDIR}/gene_counts.matrix.tsv"

if head -n 1 "${COUNTS}" | grep -q "^# Program:featureCounts"; then
  echo "[INFO] Detected raw featureCounts output. Converting to clean matrix..."
  echo "[INFO] Writing: ${MATRIX_OUT}"

  # Convert featureCounts output -> matrix with:
  #   Geneid + counts (fields 7..NF)
  #   BAM headers renamed to SRR IDs
  #   __meta rows removed
  awk -F'\t' -v OFS='\t' '
  BEGIN { header_done=0 }
  /^#/ { next }  # skip comment lines

  header_done==0 {
    # header line: keep Geneid (1) and counts columns (7..NF)
    out[1]=$1
    j=2
    for(i=7;i<=NF;i++){
      # convert "/.../star/SRRxxxx/Aligned...bam" -> "SRRxxxx"
      n=split($i, parts, "/")
      srr=parts[n-1]
      out[j]=srr
      j++
    }
    # print header
    for(k=1;k<j;k++){
      printf "%s%s", out[k], (k==j-1?ORS:OFS)
    }
    header_done=1
    next
  }

  # data lines
  {
    if($1 ~ /^__/) next
    printf "%s", $1
    for(i=7;i<=NF;i++){
      printf "%s%s", OFS, $i
    }
    printf "%s", ORS
  }
  ' "${COUNTS}" > "${MATRIX_OUT}"

  # quick sanity checks
  echo "[INFO] Clean matrix header columns:"
  head -n 1 "${MATRIX_OUT}" | awk -F'\t' '{print "NF="NF}'
  echo "[INFO] First data line columns:"
  sed -n '2p' "${MATRIX_OUT}" | awk -F'\t' '{print "NF="NF}'

  INPUT_FOR_VALIDATION="${MATRIX_OUT}"
else
  echo "[INFO] Input does not look like raw featureCounts output. Using it directly."
fi

# =========================
# RUN VALIDATION
# =========================
python3 validate_counts.py \
  --counts "${INPUT_FOR_VALIDATION}" \
  --outdir "${OUTDIR}" \
  --sep $'\t' \
  --drop-non-gene-rows \
  --collapse-duplicates \
  --outlier-factor 5

echo "[DONE] Step 1 completed."
echo "Outputs in: ${OUTDIR}"
