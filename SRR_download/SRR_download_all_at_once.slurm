#!/bin/bash
#SBATCH --job-name=sra_download
#SBATCH --output=/orange/qsong1/Yansheng/logs/%x_%j.out
#SBATCH --error=/orange/qsong1/Yansheng/logs/%x_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=6
#SBATCH --mem=32G
#SBATCH --mail-user=yanshengluo@ufl.edu
#SBATCH --mail-type=BEGIN,END,FAIL

set -euo pipefail

# =========================
# INPUT ARGUMENT
# =========================
if [[ $# -ne 1 ]]; then
  echo "USAGE: sbatch download_fastq_gse.sbatch GSEXXXX"
  exit 1
fi

GSE=$1

# =========================
# PATHS
# =========================
BASE="/orange/qsong1/Yansheng"

LIST="${BASE}/01_raw/sra/${GSE}/SRR_Acc_List.txt"
SRADIR="${BASE}/01_raw/sra/${GSE}"
FASTQDIR="${BASE}/01_raw/fastq/${GSE}"
LOGDIR="${SRADIR}/logs"

mkdir -p "$SRADIR" "$FASTQDIR" "$LOGDIR"

# =========================
# SCRATCH
# =========================
SCRATCH_BASE="${SLURM_TMPDIR:-${TMPDIR:-/tmp}}"
TMPBASE="${SCRATCH_BASE}/sratmp_${SLURM_JOB_ID}"

mkdir -p "$TMPBASE"

# =========================
# LOAD SRA TOOLS
# =========================
source /blue/wenjunxie/yanshengluo/miniconda3/etc/profile.d/conda.sh
conda activate sratools

echo "Processing dataset ${GSE}"
echo "SRR list: ${LIST}"
echo "Scratch TMPBASE: ${TMPBASE}"
date

# =========================
# MAIN LOOP
# =========================
while read -r SRR; do

  SRR=$(echo "$SRR" | tr -d '\r')

  if [[ -z "${SRR}" ]]; then
    continue
  fi

  echo "=== ${SRR} ==="

  TMPDIR="${TMPBASE}/${SRR}"
  mkdir -p "$TMPDIR"

  OUTTMP="${TMPDIR}/fastq_out"
  mkdir -p "$OUTTMP"

  cd "$SRADIR"

  # =========================
  # DOWNLOAD SRA
  # =========================
  if ! prefetch "$SRR" --output-directory "$SRADIR"; then
    echo "${SRR}" >> "${LOGDIR}/failed_srrs.txt"
    rm -rf "$TMPDIR"
    continue
  fi

  SRAFILE=$(find "$SRADIR" -maxdepth 3 -type f -name "${SRR}.sra" | head -n 1)

  if [[ -z "${SRAFILE}" ]]; then
    echo "Missing .sra for ${SRR}"
    echo "${SRR}" >> "${LOGDIR}/failed_srrs.txt"
    rm -rf "$TMPDIR"
    continue
  fi

  # =========================
  # CONVERT TO FASTQ
  # =========================
  if ! fasterq-dump "$SRAFILE" \
        -O "$OUTTMP" \
        -t "$TMPDIR" \
        --split-files \
        -e "$SLURM_CPUS_PER_TASK"; then
    echo "fasterq-dump failed for ${SRR}"
    echo "${SRR}" >> "${LOGDIR}/failed_srrs.txt"
    rm -f "$SRAFILE"
    rm -rf "$TMPDIR"
    continue
  fi

  # =========================
  # COMPRESS FASTQ
  # =========================
  if ls "${OUTTMP}/${SRR}"*.fastq 1>/dev/null 2>&1; then
    gzip -f "${OUTTMP}/${SRR}"*.fastq
  else
    echo "No FASTQ output for ${SRR}"
    echo "${SRR}" >> "${LOGDIR}/failed_srrs.txt"
    rm -f "$SRAFILE"
    rm -rf "$TMPDIR"
    continue
  fi

  # =========================
  # MOVE TO FINAL LOCATION
  # =========================
  if ls "${OUTTMP}/${SRR}"*.fastq.gz 1>/dev/null 2>&1; then
    mv -f "${OUTTMP}/${SRR}"*.fastq.gz "$FASTQDIR"/
  else
    echo "No gz FASTQ output for ${SRR}"
    echo "${SRR}" >> "${LOGDIR}/failed_srrs.txt"
    rm -f "$SRAFILE"
    rm -rf "$TMPDIR"
    continue
  fi

  # =========================
  # CLEANUP
  # =========================
  rm -f "$SRAFILE"
  rm -rf "$TMPDIR"

  echo "${SRR}" >> "${LOGDIR}/success_srrs.txt"

  echo "DONE ${SRR}"

done < "$LIST"

date
echo "Dataset ${GSE} completed"