#!/bin/bash
# ============================================================
# 02_featurecounts_merge.slurm
#
# PURPOSE
#   1) Verify STAR outputs are complete for every sample:
#        - BAM exists and passes samtools quickcheck
#        - Log.final.out exists and contains completion markers
#   2) Aggregate alignment QC into tables for review
#   3) Run featureCounts ONCE to produce a gene-level count matrix
#
# INPUT
#   sbatch 02_featurecounts_merge.slurm <FASTQ_DIR> [auto|paired|single]
#
#   Optional layout override:
#     auto   (default) : detect from BAM flags and enforce not-mixed
#     paired           : force paired-end counting and verify BAMs are paired
#     single           : force single-end counting and verify BAMs are single
#
# OUTPUTS (per dataset)
#   03_counts/<DATASET>/featurecounts/gene_counts.tsv
#   03_counts/<DATASET>/featurecounts/gene_counts.tsv.summary
#   03_counts/<DATASET>/qc_alignment/alignment_stats.tsv
#   03_counts/<DATASET>/qc_alignment/log_final_metrics.tsv
#   03_counts/<DATASET>/qc_alignment/alignment_outliers.tsv
#   03_counts/<DATASET>/qc_alignment/alignment_catastrophic.tsv
#   03_counts/<DATASET>/qc_alignment/alignment_summary.txt
#   03_counts/<DATASET>/metadata/library_layout.txt
#
# NOTES
#   - Species detection:
#       mouse default for GSExxxx
#       human for GSExxxx_HUMAN (case-insensitive)
#   - Library layout detection:
#       uses BAM flags (samtools view -c -f 1) to detect paired-end
#       hard-stop if mixed PE + SE within one dataset
#   - Paired-end counting:
#       uses -p + --countReadPairs (count fragments, not individual reads)
#       plus -B -C (recommended for standard paired RNA-seq)
# ============================================================

#SBATCH --job-name=featurecounts
#SBATCH --account=qsong1
#SBATCH --output=/orange/qsong1/Yansheng/logs/featurecounts_%j.out
#SBATCH --error=/orange/qsong1/Yansheng/logs/featurecounts_%j.err
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --mail-user=yanshengluo@ufl.edu
#SBATCH --mail-type=BEGIN,END,FAIL

set -euo pipefail

# -------------------------
# Args
# -------------------------
if [[ $# -lt 1 || $# -gt 2 ]]; then
  echo "USAGE: sbatch 02_featurecounts_merge.slurm <FASTQ_DIR> [auto|paired|single]" >&2
  echo "  Optional layout override:" >&2
  echo "    auto   (default) : detect from BAM flags and enforce not-mixed" >&2
  echo "    paired           : force paired-end counting (fragments) and verify BAMs are paired" >&2
  echo "    single           : force single-end counting and verify BAMs are single" >&2
  exit 1
fi

DATA_DIR=$(realpath "$1")
DATASET=$(basename "${DATA_DIR}")

FORCE_LAYOUT="${2:-auto}"
case "${FORCE_LAYOUT}" in
  auto|paired|single) ;;
  *)
    echo "ERROR: invalid layout override: ${FORCE_LAYOUT} (use auto|paired|single)" >&2
    exit 1
    ;;
esac

PROJECT_ROOT="/orange/qsong1/Yansheng"
OUT_ROOT="${PROJECT_ROOT}/03_counts/${DATASET}"
MANIFEST="${OUT_ROOT}/manifest/fastq_manifest.tsv"
STAR_ROOT="${OUT_ROOT}/star"
ALN_QC="${OUT_ROOT}/qc_alignment"
FC_ROOT="${OUT_ROOT}/featurecounts"
LOG_DIR="${PROJECT_ROOT}/logs"
META_DIR="${OUT_ROOT}/metadata"

# -------------------------
# Reference selection (mouse default, HUMAN override)
# -------------------------
DATASET_LOWER=$(echo "${DATASET}" | tr '[:upper:]' '[:lower:]')

REF_MOUSE="${PROJECT_ROOT}/00_metadata/reference/mm39_gencodeM38"
REF_HUMAN="${PROJECT_ROOT}/00_metadata/reference/GRCh38_gencode49"

if [[ "${DATASET_LOWER}" == *_human ]]; then
  REF_ROOT="${REF_HUMAN}"
  SPECIES="human (GRCh38 + GENCODE v49)"
else
  REF_ROOT="${REF_MOUSE}"
  SPECIES="mouse (mm39 + GENCODE M38)"
fi

GTF_GZ="${REF_ROOT}/genes.gtf.gz"

mkdir -p "${LOG_DIR}" "${ALN_QC}" "${FC_ROOT}" "${META_DIR}"

module load samtools
module load subread

echo "Dataset             : ${DATASET}"
echo "FASTQ dir (tag)     : ${DATA_DIR}"
echo "Detected species    : ${SPECIES}"
echo "Using reference     : ${REF_ROOT}"
echo "Layout override     : ${FORCE_LAYOUT}"
echo

[[ -s "${MANIFEST}" ]] || { echo "ERROR: missing manifest: ${MANIFEST}" >&2; exit 1; }
[[ -s "${GTF_GZ}" ]] || { echo "ERROR: missing GTF: ${GTF_GZ}" >&2; exit 1; }

NSAMPLES=$(( $(wc -l < "${MANIFEST}") - 1 ))
echo "Expected samples    : ${NSAMPLES}"
echo

# -------------------------
# Collect BAMs and alignment rows
# -------------------------
BAMS=()

echo -e "sample\tinput_fastq\tuniq_map_pct\tmulti_map_pct\ttoo_short_pct\tn_reads" > "${ALN_QC}/alignment_stats.tsv"

LOGFINAL_TSV="${ALN_QC}/log_final_metrics.tsv"
echo -e "sample\tuniq_map_pct\tmulti_map_pct\ttoo_short_pct\tn_reads" > "${LOGFINAL_TSV}"

while IFS=$'\t' read -r sample fq1 fq2; do
  [[ "${sample}" == "sample_id" ]] && continue

  BAM="${STAR_ROOT}/${sample}/Aligned.sortedByCoord.out.bam"
  ROW="${STAR_ROOT}/${sample}/alignment_row.tsv"
  LOGFINAL="${STAR_ROOT}/${sample}/Log.final.out"

  [[ -s "${BAM}" ]] || { echo "ERROR: missing BAM for ${sample}: ${BAM}" >&2; exit 1; }
  samtools quickcheck -v "${BAM}" || { echo "ERROR: BAM fails quickcheck: ${BAM}" >&2; exit 1; }

  [[ -s "${ROW}" ]] || { echo "ERROR: missing alignment row for ${sample}: ${ROW}" >&2; exit 1; }

  [[ -s "${LOGFINAL}" ]] || { echo "ERROR: missing Log.final.out for ${sample}: ${LOGFINAL}" >&2; exit 1; }
  grep -q "Uniquely mapped reads %" "${LOGFINAL}" \
    || { echo "ERROR: Log.final.out incomplete/corrupt for ${sample}: ${LOGFINAL}" >&2; exit 1; }

  cat "${ROW}" >> "${ALN_QC}/alignment_stats.tsv"

  UM=$(awk -F'\\|' '/Uniquely mapped reads %/ {gsub(/[%[:space:]]/,"",$2); print $2}' "${LOGFINAL}")
  MM=$(awk -F'\\|' '/% of reads mapped to multiple loci/ {gsub(/[%[:space:]]/,"",$2); print $2}' "${LOGFINAL}")
  TS=$(awk -F'\\|' '/% of reads unmapped: too short/ {gsub(/[%[:space:]]/,"",$2); print $2}' "${LOGFINAL}")
  RD=$(awk -F'\\|' '/Number of input reads/ {gsub(/[[:space:]]/,"",$2); print $2}' "${LOGFINAL}")

  [[ -n "${UM}" && -n "${MM}" && -n "${TS}" && -n "${RD}" ]] \
    || { echo "ERROR: failed to parse Log.final.out metrics for ${sample}" >&2; exit 1; }

  echo -e "${sample}\t${UM}\t${MM}\t${TS}\t${RD}" >> "${LOGFINAL_TSV}"

  BAMS+=("${BAM}")
done < "${MANIFEST}"

echo "BAMs found          : ${#BAMS[@]}"
echo

# -------------------------
# Detect library layout from BAM flags (robust)
# -------------------------
echo "[CHECK] Detecting library layout from BAM flags..."
pe=0
se=0
pe_samples=()
se_samples=()

for bam in "${BAMS[@]}"; do
  s=$(basename "$(dirname "${bam}")")

  # Count records with paired flag (0x1). If >0, paired-end reads exist in this BAM.
  npaired=$(samtools view -c -f 1 -@ 1 "${bam}" 2>/dev/null || echo 0)

  if (( npaired > 0 )); then
    pe=$((pe+1))
    pe_samples+=("${s}")
  else
    se=$((se+1))
    se_samples+=("${s}")
  fi
done

echo "[CHECK] paired-end BAMs     : ${pe}"
echo "[CHECK] single-end BAMs     : ${se}"

AUTO_LAYOUT="SINGLE_END"
if (( pe > 0 && se > 0 )); then
  AUTO_LAYOUT="MIXED"
elif (( pe > 0 )); then
  AUTO_LAYOUT="PAIRED_END"
else
  AUTO_LAYOUT="SINGLE_END"
fi

LIB_LAYOUT="${AUTO_LAYOUT}"
FC_MODE_FLAGS=()

if [[ "${FORCE_LAYOUT}" == "auto" ]]; then
  if [[ "${AUTO_LAYOUT}" == "MIXED" ]]; then
    echo "ERROR: Mixed paired-end and single-end BAMs detected in this dataset." >&2
    echo "       Fix upstream manifest/STAR step or split the dataset." >&2
    echo "       Example paired sample: ${pe_samples[0]:-NA}" >&2
    echo "       Example single sample: ${se_samples[0]:-NA}" >&2
    exit 1
  fi

  if [[ "${AUTO_LAYOUT}" == "PAIRED_END" ]]; then
    echo "[CHECK] Auto-detected paired-end: using -p --countReadPairs -B -C"
    FC_MODE_FLAGS=(-p --countReadPairs -B -C)
  else
    echo "[CHECK] Auto-detected single-end: using single-end counting"
  fi

elif [[ "${FORCE_LAYOUT}" == "paired" ]]; then
  LIB_LAYOUT="PAIRED_END"
  echo "[CHECK] Forced paired-end counting: using -p --countReadPairs -B -C"
  FC_MODE_FLAGS=(-p --countReadPairs -B -C)

  if (( pe == 0 )); then
    echo "ERROR: Forced 'paired' but no BAM shows paired-end flag (0x1)." >&2
    echo "       Example single sample: ${se_samples[0]:-NA}" >&2
    exit 1
  fi
  if (( se > 0 )); then
    echo "ERROR: Forced 'paired' but some BAMs look single-end (mixed dataset)." >&2
    echo "       Example paired sample: ${pe_samples[0]:-NA}" >&2
    echo "       Example single sample: ${se_samples[0]:-NA}" >&2
    exit 1
  fi

elif [[ "${FORCE_LAYOUT}" == "single" ]]; then
  LIB_LAYOUT="SINGLE_END"
  echo "[CHECK] Forced single-end counting"

  if (( se == 0 )); then
    echo "ERROR: Forced 'single' but all BAMs show paired-end flag (0x1)." >&2
    echo "       Example paired sample: ${pe_samples[0]:-NA}" >&2
    exit 1
  fi
  if (( pe > 0 )); then
    echo "ERROR: Forced 'single' but some BAMs look paired-end (mixed dataset)." >&2
    echo "       Example paired sample: ${pe_samples[0]:-NA}" >&2
    echo "       Example single sample: ${se_samples[0]:-NA}" >&2
    exit 1
  fi
fi

echo "[CHECK] Auto-detected layout : ${AUTO_LAYOUT}"
echo "[CHECK] Requested override   : ${FORCE_LAYOUT}"
echo "[CHECK] Effective layout     : ${LIB_LAYOUT}"
echo

# Extra guard: prevent accidental single-end counting on paired data
if (( pe > 0 )) && [[ "${LIB_LAYOUT}" != "PAIRED_END" ]]; then
  echo "ERROR: BAMs contain paired-end reads, but effective LIB_LAYOUT is not PAIRED_END." >&2
  echo "       This would count paired reads in single-end mode (WRONG)." >&2
  exit 1
fi

# -------------------------
# Log library layout metadata
# -------------------------
LIB_LAYOUT_FILE="${META_DIR}/library_layout.txt"

{
  echo "Dataset                    : ${DATASET}"
  echo "Species                    : ${SPECIES}"
  echo "Library layout (effective) : ${LIB_LAYOUT}"
  echo "Library layout (auto)      : ${AUTO_LAYOUT}"
  echo "Library layout (override)  : ${FORCE_LAYOUT}"
  echo "Total samples              : ${NSAMPLES}"
  echo "Paired-end count           : ${pe}"
  echo "Single-end count           : ${se}"
  echo "Reference used             : ${REF_ROOT}"
  echo "Gene
