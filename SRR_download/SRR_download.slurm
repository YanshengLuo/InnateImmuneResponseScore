#!/bin/bash
#SBATCH --job-name=GSE264344_sra_range
#SBATCH --output=/orange/qsong1/Yansheng/01_raw/sra/GSE264344/logs/%x_%j.out
#SBATCH --error=/orange/qsong1/Yansheng/01_raw/sra/GSE264344/logs/%x_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --mail-user=yanshengluo@ufl.edu
#SBATCH --mail-type=BEGIN,END,FAIL 

set -euo pipefail

# =========================
# INPUT ARGUMENTS
# =========================
if [[ $# -ne 2 ]]; then
  echo "USAGE: sbatch download_fastq_range.sbatch START END"
  exit 1
fi

START=$1
END=$2

# =========================
# PATHS
# =========================
LIST="/orange/qsong1/Yansheng/01_raw/sra/GSE264344/SRR_Acc_List.txt"
SRADIR="/orange/qsong1/Yansheng/01_raw/sra/GSE264344"
FASTQDIR="/orange/qsong1/Yansheng/01_raw/fastq/GSE264344"
LOGDIR="${SRADIR}/logs"
TMPBASE="${SLURM_TMPDIR}/sratmp_${SLURM_JOB_ID}"

mkdir -p "$TMPBASE"

source /blue/wenjunxie/yanshengluo/miniconda3/etc/profile.d/conda.sh
conda activate sratools

echo "Processing SRR list lines ${START} to ${END}"
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

  echo "=== [${IDX}] ${SRR} ==="

  TMPDIR="${TMPBASE}/${SRR}"
  mkdir -p "$TMPDIR"

  # 1) Download .sra
  cd "$SRADIR"
  if ! prefetch "$SRR" --output-directory "$SRADIR"; then
    echo "${SRR}" >> "${LOGDIR}/failed_srrs.txt"
    continue
  fi

  SRAFILE=$(find "$SRADIR" -maxdepth 3 -type f -name "${SRR}.sra" | head -n 1)
  if [[ -z "${SRAFILE}" ]]; then
    echo "Missing .sra for ${SRR}"
    echo "${SRR}" >> "${LOGDIR}/failed_srrs.txt"
    continue
  fi

  # 2) Convert to FASTQ (bulk RNA-seq safe)
  if ! fasterq-dump "$SRAFILE" \
        -O "$FASTQDIR" \
        -t "$TMPDIR" \
        --split-files \
        -e "$SLURM_CPUS_PER_TASK"; then
    echo "fasterq-dump failed for ${SRR}"
    echo "${SRR}" >> "${LOGDIR}/failed_srrs.txt"
    continue
  fi

  # 3) Compress FASTQs
  if ls "${FASTQDIR}/${SRR}"*.fastq 1>/dev/null 2>&1; then
    gzip -f "${FASTQDIR}/${SRR}"*.fastq
  else
    echo "No FASTQ output for ${SRR}"
    echo "${SRR}" >> "${LOGDIR}/failed_srrs.txt"
    continue
  fi

  # 4) Cleanup
  rm -f "$SRAFILE"
  rm -rf "$TMPDIR"

  echo "${SRR}" >> "${LOGDIR}/success_srrs.txt"
  echo "DONE ${SRR}"

done

date
echo "Range ${START}-${END} completed"
