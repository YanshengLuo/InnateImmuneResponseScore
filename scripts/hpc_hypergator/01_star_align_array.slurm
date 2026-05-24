#!/bin/bash
#SBATCH --job-name=star_align
#SBATCH --account=qsong1
#SBATCH --output=/orange/qsong1/Yansheng/logs/star_align_%A_%a.out
#SBATCH --error=/orange/qsong1/Yansheng/logs/star_align_%A_%a.err
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=34G
#SBATCH --mail-user=yanshengluo@ufl.edu
#SBATCH --mail-type=BEGIN,END,FAIL

set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "USAGE: sbatch --array=0-(N-1)%5 01_star_align_array.slurm <FASTQ_DIR>" >&2
  exit 1
fi

DATA_DIR=$(realpath "$1")
DATASET=$(basename "${DATA_DIR}")

PROJECT_ROOT="/orange/qsong1/Yansheng"
OUT_ROOT="${PROJECT_ROOT}/03_counts/${DATASET}"
MANIFEST_DIR="${OUT_ROOT}/manifest"
STAR_ROOT="${OUT_ROOT}/star"
ALN_QC="${OUT_ROOT}/qc_alignment"
LOG_DIR="${PROJECT_ROOT}/logs"

REF_ROOT="${PROJECT_ROOT}/00_metadata/reference/mm39_gencodeM38"
STAR_INDEX="${REF_ROOT}/STAR_index"

mkdir -p "${LOG_DIR}" "${OUT_ROOT}" "${MANIFEST_DIR}" "${STAR_ROOT}" "${ALN_QC}"

module load star
module load samtools

MANIFEST="${MANIFEST_DIR}/fastq_manifest.tsv"
LOCKFILE="${MANIFEST_DIR}/.manifest.lock"

# ------------------------------------------------------------
# QC pass/fail logging (array-safe)
# ------------------------------------------------------------
FAILED_LIST="${ALN_QC}/failed_samples.txt"
PASSED_LIST="${ALN_QC}/passed_samples.txt"
QC_LOCK="${ALN_QC}/.qc.lock"

append_qc_safe() {
  # $1 = "fail" or "pass"
  # $2 = reason (for fail)
  local mode="$1"
  local reason="${2:-}"
  (
    flock -x 9
    if [[ "${mode}" == "fail" ]]; then
      echo -e "${SAMPLE}\t${reason}\t${SDIR}" >> "${FAILED_LIST}"
    else
      echo -e "${SAMPLE}\tPASS\t${SDIR}" >> "${PASSED_LIST}"
    fi
  ) 9>"${QC_LOCK}"
}

fail_with_reason() {
  local reason="$1"
  echo "${reason}" > "${SDIR}/FAIL_REASON.txt" 2>/dev/null || true
  append_qc_safe "fail" "${reason}"
  echo "FAIL: ${SAMPLE} — ${reason}" >&2
  exit 1
}

# Log unexpected errors ONLY after SAMPLE/SDIR exist
on_error() {
  local exit_code=$?
  if [[ -n "${SAMPLE:-}" && -n "${SDIR:-}" ]]; then
    echo "script_error_exit_${exit_code}" > "${SDIR}/FAIL_REASON.txt" 2>/dev/null || true
    append_qc_safe "fail" "script_error_exit_${exit_code}"
  fi
  exit "${exit_code}"
}
trap on_error ERR

# ------------------------------------------------------------
# Build manifest once (task 0 does it; others wait)
# Manifest format:
# sample_id<TAB>fq1<TAB>fq2_or_dot
# ------------------------------------------------------------
build_manifest() {
  # Detect paired-end by presence of *_1 or *_R1 files
  local n_r1 n_r2
  n_r1=$(find "${DATA_DIR}" -maxdepth 1 -type f \( -name "*_1.fastq.gz" -o -name "*_R1*.fastq.gz" \) | wc -l || true)
  n_r2=$(find "${DATA_DIR}" -maxdepth 1 -type f \( -name "*_2.fastq.gz" -o -name "*_R2*.fastq.gz" \) | wc -l || true)

  echo -e "sample_id\tfq1\tfq2" > "${MANIFEST}"

  if [[ "${n_r1}" -gt 0 && "${n_r2}" -gt 0 ]]; then
    echo "Manifest: paired-end detected (heuristic). R1=${n_r1}, R2=${n_r2}"
    mapfile -t R1S < <(find "${DATA_DIR}" -maxdepth 1 -type f \( -name "*_1.fastq.gz" -o -name "*_R1*.fastq.gz" \) | sort)
    for R1 in "${R1S[@]}"; do
      R2="${R1/_1.fastq.gz/_2.fastq.gz}"
      R2="${R2/_R1/_R2}"
      [[ -f "${R2}" ]] || { echo "ERROR: missing R2 for ${R1}" >&2; exit 1; }

      base=$(basename "${R1}" .fastq.gz)
      base=${base%_1}
      base=${base%_R1}

      echo -e "${base}\t${R1}\t${R2}" >> "${MANIFEST}"
    done
  else
    echo "Manifest: single-end detected (heuristic)."
    mapfile -t FQS < <(find "${DATA_DIR}" -maxdepth 1 -type f -name "*.fastq.gz" | sort)
    ((${#FQS[@]} > 0)) || { echo "ERROR: no *.fastq.gz in ${DATA_DIR}" >&2; exit 1; }
    for FQ in "${FQS[@]}"; do
      base=$(basename "${FQ}" .fastq.gz)
      echo -e "${base}\t${FQ}\t." >> "${MANIFEST}"
    done
  fi
}

# Task 0 creates manifest under a lock; others wait for it.
if [[ "${SLURM_ARRAY_TASK_ID}" -eq 0 ]]; then
  (
    flock -n 9 || exit 0
    if [[ ! -s "${MANIFEST}" ]]; then
      echo "Building manifest at ${MANIFEST}"
      build_manifest
    fi
  ) 9>"${LOCKFILE}"
fi

# Wait until manifest exists (all tasks)
for i in $(seq 1 60); do
  [[ -s "${MANIFEST}" ]] && break
  sleep 2
done
[[ -s "${MANIFEST}" ]] || { echo "ERROR: manifest not created: ${MANIFEST}" >&2; exit 1; }

# How many samples?
NSAMPLES=$(( $(wc -l < "${MANIFEST}") - 1 ))
if [[ "${SLURM_ARRAY_TASK_ID}" -ge "${NSAMPLES}" ]]; then
  echo "Nothing to do: task_id ${SLURM_ARRAY_TASK_ID} >= NSAMPLES ${NSAMPLES}"
  exit 0
fi

# Fetch the row for this task (skip header)
ROW=$(awk -v n=$((SLURM_ARRAY_TASK_ID+2)) 'NR==n {print; exit}' "${MANIFEST}")
SAMPLE=$(echo "${ROW}" | cut -f1)
FQ1=$(echo "${ROW}" | cut -f2)
FQ2=$(echo "${ROW}" | cut -f3)

SDIR="${STAR_ROOT}/${SAMPLE}"
mkdir -p "${SDIR}"

echo "=== STAR align ==="
date
echo "Dataset : ${DATASET}"
echo "Sample  : ${SAMPLE}"
echo "FQ1     : ${FQ1}"
echo "FQ2     : ${FQ2}"
echo "Index   : ${STAR_INDEX}"
echo "Threads : ${SLURM_CPUS_PER_TASK}"
echo

[[ -s "${STAR_INDEX}/Genome" ]] || { echo "ERROR: STAR index missing at ${STAR_INDEX}. Run prep script." >&2; exit 1; }

# ------------------------------------------------------------
# STAR alignment (no QC checks happen until STAR returns)
# ------------------------------------------------------------
if [[ "${FQ2}" == "." ]]; then
  STAR --runThreadN ${SLURM_CPUS_PER_TASK} \
       --genomeDir "${STAR_INDEX}" \
       --readFilesIn "${FQ1}" \
       --readFilesCommand zcat \
       --outFileNamePrefix "${SDIR}/" \
       --outSAMtype BAM SortedByCoordinate
else
  STAR --runThreadN ${SLURM_CPUS_PER_TASK} \
       --genomeDir "${STAR_INDEX}" \
       --readFilesIn "${FQ1}" "${FQ2}" \
       --readFilesCommand zcat \
       --outFileNamePrefix "${SDIR}/" \
       --outSAMtype BAM SortedByCoordinate
fi

# ------------------------------------------------------------
# POST-STAR validation ONLY (runs after STAR is fully done)
# ------------------------------------------------------------
BAM="${SDIR}/Aligned.sortedByCoord.out.bam"
LOGF="${SDIR}/Log.final.out"

[[ -s "${BAM}" ]]  || fail_with_reason "missing_or_empty_bam"
[[ -s "${LOGF}" ]] || fail_with_reason "missing_or_empty_Log.final.out"

samtools quickcheck -v "${BAM}" >/dev/null 2>&1 || fail_with_reason "bam_quickcheck_failed"

grep -q "Uniquely mapped reads %" "${LOGF}" || fail_with_reason "log_missing_uniquely_mapped_marker"
grep -q "Finished on" "${LOGF}"            || fail_with_reason "log_missing_finished_on_marker"

# mark success (only after passing checks)
append_qc_safe "pass"

# proceed with indexing + stats
samtools index "${BAM}"

# Append alignment stats (safe-ish in arrays: we write per-sample then merge later)
# Write a per-sample stats file; merge happens in featureCounts stage.
NREADS=$(awk -F'\t' '$1 ~ /Number of input reads/ {gsub(/ /,"",$2); print $2}' "${LOGF}")
UMAP=$(awk -F'\t' '$1 ~ /Uniquely mapped reads %/ {gsub(/%| /,"",$2); print $2}' "${LOGF}")
MMAP=$(awk -F'\t' '$1 ~ /% of reads mapped to multiple loci/ {gsub(/%| /,"",$2); print $2}' "${LOGF}")
TSHO=$(awk -F'\t' '$1 ~ /% of reads unmapped: too short/ {gsub(/%| /,"",$2); print $2}' "${LOGF}")

echo -e "${SAMPLE}\t${FQ1}\t${UMAP}\t${MMAP}\t${TSHO}\t${NREADS}" > "${SDIR}/alignment_row.tsv"

echo "DONE: ${SAMPLE}"
date
