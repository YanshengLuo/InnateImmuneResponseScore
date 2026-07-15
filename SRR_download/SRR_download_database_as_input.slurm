#!/bin/bash
#SBATCH --job-name=sra_range
#SBATCH --output=/orange/qsong1/Yansheng/01_raw/sra/logs/%x_%j.out 
#SBATCH --error=/orange/qsong1/Yansheng/01_raw/sra/logs/%x_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=6
#SBATCH --mem=32G
#SBATCH --mail-user=yanshengluo@ufl.edu
#SBATCH --mail-type=BEGIN,END,FAIL
# NOTE: #SBATCH --output/--error cannot reliably use variables.
# We write per-dataset logs ourselves under ${SRADIR}/logs.

set -euo pipefail

# ============================================================
# USAGE
#   sbatch SRR_download.slurm <GSE_ID> [START END]
#
# BEHAVIOR
#   1) GSE_ID is REQUIRED. If missing -> prompt + exit.
#   2) If START/END are missing -> download ALL (1..N lines in SRR list).
# ============================================================

# ---------- Require dataset id ----------
if [[ $# -lt 1 ]]; then
  echo "ERROR: Missing dataset id (e.g., GSE264344)."
  echo "USAGE: sbatch SRR_download.slurm <GSE_ID> [START END]"
  exit 1
fi

DATASET="$1"

# ---------- Defaults: download all if no range provided ----------
LIST="/orange/qsong1/Yansheng/01_raw/sra/${DATASET}/SRR_Acc_List.txt"
if [[ ! -s "$LIST" ]]; then
  echo "ERROR: SRR list not found or empty:"
  echo "  $LIST"
  exit 2
fi

NSAMPLES="$(wc -l < "$LIST" | tr -d ' ')"

if [[ $# -eq 1 ]]; then
  START=1
  END="$NSAMPLES"
elif [[ $# -eq 3 ]]; then
  START="$2"
  END="$3"
else
  echo "ERROR: Invalid arguments."
  echo "USAGE: sbatch SRR_download.slurm <GSE_ID> [START END]"
  exit 3
fi

# ---------- Validate START/END ----------
if ! [[ "$START" =~ ^[0-9]+$ && "$END" =~ ^[0-9]+$ ]]; then
  echo "ERROR: START/END must be integers."
  exit 4
fi
if (( START < 1 )); then
  echo "ERROR: START must be >= 1"
  exit 5
fi
if (( END < START )); then
  echo "ERROR: END must be >= START"
  exit 6
fi
if (( END > NSAMPLES )); then
  echo "NOTE: END (${END}) > NSAMPLES (${NSAMPLES}); clamping END to ${NSAMPLES}"
  END="$NSAMPLES"
fi

# =========================
# PATHS (auto-recognized by dataset id)
# =========================
SRADIR="/orange/qsong1/Yansheng/01_raw/sra/${DATASET}"
FASTQDIR="/orange/qsong1/Yansheng/01_raw/fastq/${DATASET}"
LOGDIR="${SRADIR}/logs"

# Use node-local scratch if available (fast), fallback to /tmp
SCRATCH_BASE="${SLURM_TMPDIR:-${TMPDIR:-/tmp}}"
TMPBASE="${SCRATCH_BASE}/sratmp_${SLURM_JOB_ID}"

mkdir -p "$TMPBASE" "$LOGDIR" "$FASTQDIR"

# Optional: dataset-scoped job log that doesn't rely on #SBATCH --output
JOBLOG="${LOGDIR}/${DATASET}_job_${SLURM_JOB_ID}.log"
exec > >(tee -a "$JOBLOG") 2>&1

source /blue/wenjunxie/yanshengluo/miniconda3/etc/profile.d/conda.sh
conda activate sratools

echo "DATASET: ${DATASET}"
echo "SRR LIST: ${LIST}"
echo "NSAMPLES: ${NSAMPLES}"
echo "Processing lines ${START} to ${END}"
echo "Scratch TMPBASE: ${TMPBASE}"
echo "FASTQDIR: ${FASTQDIR}"
date

# =========================
# MAIN LOOP (bulk RNA-seq)
# =========================
for IDX in $(seq "$START" "$END"); do

  SRR=$(sed -n "${IDX}p" "$LIST" | tr -d '\r' | awk '{print $1}')

  if [[ -z "${SRR}" ]]; then
    echo "Line ${IDX}: empty or out of range — skipping"
    continue
  fi

  echo "=== [${IDX}/${NSAMPLES}] ${SRR} ==="

  # Skip if already downloaded (either SE or PE exists)
  if [[ -s "${FASTQDIR}/${SRR}_1.fastq.gz" || -s "${FASTQDIR}/${SRR}.fastq.gz" || -s "${FASTQDIR}/${SRR}_R1.fastq.gz" ]]; then
    echo "SKIP: ${SRR} appears already downloaded in ${FASTQDIR}"
    echo "${SRR}" >> "${LOGDIR}/skipped_srrs.txt"
    continue
  fi

  # Per-SRR scratch workspace
  TMPDIR="${TMPBASE}/${SRR}"
  mkdir -p "$TMPDIR"

  OUTTMP="${TMPDIR}/fastq_out"
  mkdir -p "$OUTTMP"

  # 1) Download .sra
  cd "$SRADIR"
  if ! prefetch "$SRR" --output-directory "$SRADIR"; then
    echo "prefetch failed for ${SRR}"
    echo "${SRR}" >> "${LOGDIR}/failed_srrs.txt"
    rm -rf "$TMPDIR"
    continue
  fi

  SRAFILE=$(find "$SRADIR" -maxdepth 3 -type f -name "${SRR}.sra" | head -n 1 || true)
  if [[ -z "${SRAFILE}" ]]; then
    echo "Missing .sra for ${SRR}"
    echo "${SRR}" >> "${LOGDIR}/failed_srrs.txt"
    rm -rf "$TMPDIR"
    continue
  fi

  # 2) Convert to FASTQ (safe default: split-files)
  if ! fasterq-dump "$SRAFILE" \
        -O "$OUTTMP" \
        -t "$TMPDIR" \
        --split-files \
        -e "${SLURM_CPUS_PER_TASK}"; then
    echo "fasterq-dump failed for ${SRR}"
    echo "${SRR}" >> "${LOGDIR}/failed_srrs.txt"
    rm -f "$SRAFILE"
    rm -rf "$TMPDIR"
    continue
  fi

  # 3) Compress and move
  if ls "${OUTTMP}/${SRR}"*.fastq 1>/dev/null 2>&1; then
    gzip -f "${OUTTMP}/${SRR}"*.fastq
  else
    echo "No FASTQ output for ${SRR}"
    echo "${SRR}" >> "${LOGDIR}/failed_srrs.txt"
    rm -f "$SRAFILE"
    rm -rf "$TMPDIR"
    continue
  fi

  if ls "${OUTTMP}/${SRR}"*.fastq.gz 1>/dev/null 2>&1; then
    mv -f "${OUTTMP}/${SRR}"*.fastq.gz "$FASTQDIR"/
  else
    echo "No gz FASTQ output for ${SRR}"
    echo "${SRR}" >> "${LOGDIR}/failed_srrs.txt"
    rm -f "$SRAFILE"
    rm -rf "$TMPDIR"
    continue
  fi

  # 4) Cleanup
  rm -f "$SRAFILE"
  rm -rf "$TMPDIR"

  echo "${SRR}" >> "${LOGDIR}/success_srrs.txt"
  echo "DONE ${SRR}"

done

date
echo "DATASET ${DATASET}: range ${START}-${END} completed"
echo "Job log: ${JOBLOG}"
