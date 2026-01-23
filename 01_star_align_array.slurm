#!/bin/bash
# ============================================================
# 01_star_align_array.slurm  (AUTO-RETRY + FAILURE LIST)
#
# PURPOSE
#   Align each sample (SE or PE) with STAR in a Slurm job array.
#   This version adds robustness:
#     1) Runs STAR (attempt 1)
#     2) Validates outputs (BAM + BAI + Log.final.out completion marker + samtools quickcheck)
#     3) If validation fails, retries ONCE (attempt 2) after cleaning the sample dir
#     4) If still failing, writes a per-sample failure status AND appends to a dataset-level failure list
#
# HOW TO SUBMIT (example)
#   sbatch --array=0-999%5 01_star_align_array.slurm <FASTQ_DIR>
#
# CONCURRENCY / RESOURCE CAP
#   - Total CPU cap (example): 20 CPUs
#   - Use cpus-per-task=4 and array concurrency %5 (4*5 = 20)
#   - Memory per task should be sized so concurrent tasks fit your group cap
#
# INPUT
#   FASTQ_DIR: directory containing *.fastq.gz for a dataset
#   Dataset name is inferred as: DATASET = basename(FASTQ_DIR)
#
# OUTPUTS (per dataset)
#   03_counts/<DATASET>/manifest/fastq_manifest.tsv
#   03_counts/<DATASET>/star/<SAMPLE>/
#     Aligned.sortedByCoord.out.bam
#     Aligned.sortedByCoord.out.bam.bai
#     Log.final.out
#     alignment_row.tsv
#     star_status.tsv              (attempt-level status + reason)
#   03_counts/<DATASET>/qc_alignment/
#     star_failed_samples.tsv      (aggregate list across array tasks)
#
# NOTES
#   - This script does not send emails (no personal info).
#   - It is safe to run repeatedly; successful samples will be skipped.
# ============================================================

#SBATCH --job-name=star_align
#SBATCH --account=qsong1
#SBATCH --output=/orange/qsong1/Yansheng/logs/star_align_%A_%a.out
#SBATCH --error=/orange/qsong1/Yansheng/logs/star_align_%A_%a.err
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=34G

set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "USAGE: sbatch --array=0-999%5 01_star_align_array.slurm <FASTQ_DIR>" >&2
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

FAILED_LIST="${ALN_QC}/star_failed_samples.tsv"
FAILED_LOCK="${ALN_QC}/.star_failed.lock"

# -------------------------
# Build manifest once (task 0 does it; others wait)
# Manifest columns: sample_id<TAB>fq1<TAB>fq2_or_dot
# -------------------------
build_manifest() {
  local n_r1 n_r2
  n_r1=$(find "${DATA_DIR}" -maxdepth 1 -type f \( -name "*_1.fastq.gz" -o -name "*_R1*.fastq.gz" \) | wc -l || true)
  n_r2=$(find "${DATA_DIR}" -maxdepth 1 -type f \( -name "*_2.fastq.gz" -o -name "*_R2*.fastq.gz" \) | wc -l || true)

  echo -e "sample_id\tfq1\tfq2" > "${MANIFEST}"

  if [[ "${n_r1}" -gt 0 && "${n_r2}" -gt 0 ]]; then
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
    mapfile -t FQS < <(find "${DATA_DIR}" -maxdepth 1 -type f -name "*.fastq.gz" | sort)
    ((${#FQS[@]} > 0)) || { echo "ERROR: no *.fastq.gz in ${DATA_DIR}" >&2; exit 1; }
    for FQ in "${FQS[@]}"; do
      base=$(basename "${FQ}" .fastq.gz)
      echo -e "${base}\t${FQ}\t." >> "${MANIFEST}"
    done
  fi
}

if [[ "${SLURM_ARRAY_TASK_ID}" -eq 0 ]]; then
  (
    flock -n 9 || exit 0
    if [[ ! -s "${MANIFEST}" ]]; then
      build_manifest
    fi
  ) 9>"${LOCKFILE}"
fi

for i in $(seq 1 90); do
  [[ -s "${MANIFEST}" ]] && break
  sleep 2
done
[[ -s "${MANIFEST}" ]] || { echo "ERROR: manifest not created: ${MANIFEST}" >&2; exit 1; }

NSAMPLES=$(( $(wc -l < "${MANIFEST}") - 1 ))
if [[ "${SLURM_ARRAY_TASK_ID}" -ge "${NSAMPLES}" ]]; then
  echo "Nothing to do: task_id ${SLURM_ARRAY_TASK_ID} >= NSAMPLES ${NSAMPLES}"
  exit 0
fi

ROW=$(awk -v n=$((SLURM_ARRAY_TASK_ID+2)) 'NR==n {print; exit}' "${MANIFEST}")
SAMPLE=$(echo "${ROW}" | cut -f1)
FQ1=$(echo "${ROW}" | cut -f2)
FQ2=$(echo "${ROW}" | cut -f3)

SDIR="${STAR_ROOT}/${SAMPLE}"
STATUS_TSV="${SDIR}/star_status.tsv"

mkdir -p "${SDIR}"

# If already successful, skip
if [[ -s "${SDIR}/Aligned.sortedByCoord.out.bam" && -s "${SDIR}/Aligned.sortedByCoord.out.bam.bai" && -s "${SDIR}/Log.final.out" ]]; then
  if grep -q "Uniquely mapped reads %" "${SDIR}/Log.final.out"; then
    if samtools quickcheck "${SDIR}/Aligned.sortedByCoord.out.bam" >/dev/null 2>&1; then
      echo -e "sample\tstatus\tattempt\treason" > "${STATUS_TSV}"
      echo -e "${SAMPLE}\tSUCCESS\t0\talready_present" >> "${STATUS_TSV}"
      exit 0
    fi
  fi
fi

[[ -s "${STAR_INDEX}/Genome" ]] || { echo "ERROR: STAR index missing at ${STAR_INDEX}" >&2; exit 1; }

# -------------------------
# Helpers
# -------------------------
append_failed_list() {
  local reason="$1"
  mkdir -p "${ALN_QC}"
  (
    flock 9
    if [[ ! -s "${FAILED_LIST}" ]]; then
      echo -e "sample\ttask_id\treason" > "${FAILED_LIST}"
    fi
    echo -e "${SAMPLE}\t${SLURM_ARRAY_TASK_ID}\t${reason}" >> "${FAILED_LIST}"
  ) 9>"${FAILED_LOCK}"
}

validate_outputs() {
  local bam="${SDIR}/Aligned.sortedByCoord.out.bam"
  local bai="${SDIR}/Aligned.sortedByCoord.out.bam.bai"
  local logf="${SDIR}/Log.final.out"

  [[ -s "${bam}" ]] || return 10
  [[ -s "${bai}" ]] || return 11
  [[ -s "${logf}" ]] || return 12
  grep -q "Uniquely mapped reads %" "${logf}" || return 13
  samtools quickcheck -v "${bam}" >/dev/null 2>&1 && return 0 || return 14
}

run_star_once() {
  local attempt="$1"

  # Prefer node-local scratch if available; else fall back to SDIR/_STARtmp
  local tmpbase=""
  if [[ -n "${TMPDIR:-}" && -d "${TMPDIR}" ]]; then
    tmpbase="${TMPDIR}/STARtmp_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}_a${attempt}"
  else
    tmpbase="${SDIR}/_STARtmp_a${attempt}"
  fi
  mkdir -p "${tmpbase}"

  if [[ "${FQ2}" == "." ]]; then
    STAR --runThreadN "${SLURM_CPUS_PER_TASK}" \
         --genomeDir "${STAR_INDEX}" \
         --readFilesIn "${FQ1}" \
         --readFilesCommand zcat \
         --outFileNamePrefix "${SDIR}/" \
         --outSAMtype BAM SortedByCoordinate \
         --outTmpDir "${tmpbase}"
  else
    STAR --runThreadN "${SLURM_CPUS_PER_TASK}" \
         --genomeDir "${STAR_INDEX}" \
         --readFilesIn "${FQ1}" "${FQ2}" \
         --readFilesCommand zcat \
         --outFileNamePrefix "${SDIR}/" \
         --outSAMtype BAM SortedByCoordinate \
         --outTmpDir "${tmpbase}"
  fi

  samtools index "${SDIR}/Aligned.sortedByCoord.out.bam"

  # Build alignment_row.tsv from Log.final.out (created only on completion)
  local logf="${SDIR}/Log.final.out"
  if [[ -s "${logf}" ]]; then
    NREADS=$(awk -F'\\|' '/Number of input reads/ {gsub(/[[:space:]]/,"",$2); print $2}' "${logf}" || true)
    UMAP=$(awk -F'\\|' '/Uniquely mapped reads %/ {gsub(/[%[:space:]]/,"",$2); print $2}' "${logf}" || true)
    MMAP=$(awk -F'\\|' '/% of reads mapped to multiple loci/ {gsub(/[%[:space:]]/,"",$2); print $2}' "${logf}" || true)
    TSHO=$(awk -F'\\|' '/% of reads unmapped: too short/ {gsub(/[%[:space:]]/,"",$2); print $2}' "${logf}" || true)
    echo -e "${SAMPLE}\t${FQ1}\t${UMAP}\t${MMAP}\t${TSHO}\t${NREADS}" > "${SDIR}/alignment_row.tsv"
  fi
}

clean_sample_dir() {
  # Remove partial STAR outputs before retry to avoid mixing attempt 1/2 artifacts
  rm -rf "${SDIR}"
  mkdir -p "${SDIR}"
}

# -------------------------
# Run + retry logic
# -------------------------
echo -e "sample\tstatus\tattempt\treason" > "${STATUS_TSV}"

ATTEMPT=1
REASON=""

set +e
run_star_once "${ATTEMPT}"
rc=$?
set -e

if [[ "${rc}" -ne 0 ]]; then
  REASON="star_exitcode_${rc}"
else
  if ! validate_outputs; then
    vrc=$?
    REASON="validate_fail_${vrc}"
  fi
fi

if [[ -z "${REASON}" ]]; then
  echo -e "${SAMPLE}\tSUCCESS\t1\tok" >> "${STATUS_TSV}"
  exit 0
fi

# Retry once
echo -e "${SAMPLE}\tRETRY\t1\t${REASON}" >> "${STATUS_TSV}"
clean_sample_dir

ATTEMPT=2
REASON2=""

set +e
run_star_once "${ATTEMPT}"
rc2=$?
set -e

if [[ "${rc2}" -ne 0 ]]; then
  REASON2="star_exitcode_${rc2}"
else
  if ! validate_outputs; then
    vrc2=$?
    REASON2="validate_fail_${vrc2}"
  fi
fi

if [[ -z "${REASON2}" ]]; then
  echo -e "${SAMPLE}\tSUCCESS\t2\tok_after_retry" >> "${STATUS_TSV}"
  exit 0
fi

# Still failed: record failure
echo -e "${SAMPLE}\tFAIL\t2\t${REASON2}" >> "${STATUS_TSV}"
append_failed_list "${REASON2}"
exit 1
