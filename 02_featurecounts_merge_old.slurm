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
#   sbatch 02_featurecounts_merge.slurm <FASTQ_DIR>
#   FASTQ_DIR is used ONLY to infer the dataset name (basename) and locate outputs:
#     - manifest: 03_counts/<DATASET>/manifest/fastq_manifest.tsv
#     - STAR outputs: 03_counts/<DATASET>/star/<SAMPLE>/
#
# OUTPUTS (per dataset)
#   03_counts/<DATASET>/featurecounts/gene_counts.tsv
#   03_counts/<DATASET>/featurecounts/gene_counts.tsv.summary
#   03_counts/<DATASET>/qc_alignment/alignment_stats.tsv            (from STAR alignment_row.tsv)
#   03_counts/<DATASET>/qc_alignment/log_final_metrics.tsv          (parsed from Log.final.out)
#   03_counts/<DATASET>/qc_alignment/alignment_outliers.tsv
#   03_counts/<DATASET>/qc_alignment/alignment_catastrophic.tsv
#   03_counts/<DATASET>/qc_alignment/alignment_summary.txt
#
# SLURM PARAMETERS (WHAT EACH ONE DOES)
#   --job-name        : name shown in squeue/sacct
#   --account         : which group account gets charged (important for fairshare/quota)
#   --output/--error  : log file paths; %j is the Slurm job ID
#   --time            : max walltime; kill job if exceeded
#   --cpus-per-task   : CPU threads available to featureCounts and parsing steps
#   --mem             : total RAM for this job (featureCounts + file parsing)
#   --mail-user       : email address for notifications (REMOVE if you do not want email)
#   --mail-type       : events that trigger notifications
#
# FEATURECOUNTS PARAMETERS (WHAT EACH ONE DOES)
#   -T <threads>      : threads used by featureCounts (use SLURM_CPUS_PER_TASK)
#   -a <gtf>          : gene annotation in GTF format
#   -o <outfile>      : output count matrix path
#   -g gene_id        : attribute used as gene identifier (gene_id is standard for GENCODE)
#   -t exon           : count reads overlapping exons (then summarize to gene level)
#
# QC THRESHOLDS (USED FOR FLAGS; DO NOT AUTO-DROP SAMPLES HERE)
#   UNIQ_MIN=50       : unique mapping % below this is flagged as an outlier
#   TOO_SHORT_MAX=20  : unmapped-too-short % above this is flagged
#   MULTI_MAX=20      : multi-mapped % above this is flagged
#   READS_MIN=1e6     : input reads below this is flagged
#   Catastrophic flags use stricter rules: uniq<30 OR too_short>40 OR reads<200k
#
# DEPENDENCIES
#   module load samtools
#   module load subread    (featureCounts)
#
# NOTE ON EMAIL
#   If you want to remove personal info, delete these lines:
#     #SBATCH --mail-user=...
#     #SBATCH --mail-type=...
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

if [[ $# -ne 1 ]]; then
  echo "USAGE: sbatch 02_featurecounts_merge.slurm <FASTQ_DIR>" >&2
  exit 1
fi

DATA_DIR=$(realpath "$1")
DATASET=$(basename "${DATA_DIR}")

PROJECT_ROOT="/orange/qsong1/Yansheng"
OUT_ROOT="${PROJECT_ROOT}/03_counts/${DATASET}"
MANIFEST="${OUT_ROOT}/manifest/fastq_manifest.tsv"
STAR_ROOT="${OUT_ROOT}/star"
ALN_QC="${OUT_ROOT}/qc_alignment"
FC_ROOT="${OUT_ROOT}/featurecounts"
LOG_DIR="${PROJECT_ROOT}/logs"

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
GTF_TMP="${REF_ROOT}/genes.gtf"

echo "Detected species: ${SPECIES}"
echo "Using reference: ${REF_ROOT}"
mkdir -p "${LOG_DIR}" "${ALN_QC}" "${FC_ROOT}"

module load samtools
module load subread

[[ -s "${MANIFEST}" ]] || { echo "ERROR: missing manifest: ${MANIFEST}" >&2; exit 1; }

NSAMPLES=$(( $(wc -l < "${MANIFEST}") - 1 ))
echo "Dataset: ${DATASET}"
echo "Expected samples: ${NSAMPLES}"
echo

# Collect BAMs and alignment rows
BAMS=()
echo -e "sample\tinput_fastq\tuniq_map_pct\tmulti_map_pct\ttoo_short_pct\tn_reads" > "${ALN_QC}/alignment_stats.tsv"

# Collect Log.final.out metrics independently (second source of truth)
LOGFINAL_TSV="${ALN_QC}/log_final_metrics.tsv"
echo -e "sample\tuniq_map_pct\tmulti_map_pct\ttoo_short_pct\tn_reads" > "${LOGFINAL_TSV}"

while IFS=$'\t' read -r sample fq1 fq2; do
  [[ "${sample}" == "sample_id" ]] && continue

  BAM="${STAR_ROOT}/${sample}/Aligned.sortedByCoord.out.bam"
  ROW="${STAR_ROOT}/${sample}/alignment_row.tsv"
  LOGFINAL="${STAR_ROOT}/${sample}/Log.final.out"

  # --- Sanity checks: required outputs from STAR ---
  [[ -s "${BAM}" ]] || { echo "ERROR: missing BAM for ${sample}: ${BAM}" >&2; exit 1; }
  samtools quickcheck -v "${BAM}" || { echo "ERROR: BAM fails quickcheck: ${BAM}" >&2; exit 1; }

  [[ -s "${ROW}" ]] || { echo "ERROR: missing alignment row for ${sample}: ${ROW}" >&2; exit 1; }

  [[ -s "${LOGFINAL}" ]] || { echo "ERROR: missing Log.final.out for ${sample}: ${LOGFINAL}" >&2; exit 1; }
  grep -q "Uniquely mapped reads %" "${LOGFINAL}" \
    || { echo "ERROR: Log.final.out incomplete/corrupt for ${sample}: ${LOGFINAL}" >&2; exit 1; }

  # Append the STAR-derived row (created during alignment)
  cat "${ROW}" >> "${ALN_QC}/alignment_stats.tsv"

  # Parse key metrics from Log.final.out (robust completion confirmation)
  UM=$(awk -F'\\|' '/Uniquely mapped reads %/ {gsub(/[%[:space:]]/,"",$2); print $2}' "${LOGFINAL}")
  MM=$(awk -F'\\|' '/% of reads mapped to multiple loci/ {gsub(/[%[:space:]]/,"",$2); print $2}' "${LOGFINAL}")
  TS=$(awk -F'\\|' '/% of reads unmapped: too short/ {gsub(/[%[:space:]]/,"",$2); print $2}' "${LOGFINAL}")
  RD=$(awk -F'\\|' '/Number of input reads/ {gsub(/[[:space:]]/,"",$2); print $2}' "${LOGFINAL}")

  [[ -n "${UM}" && -n "${MM}" && -n "${TS}" && -n "${RD}" ]] \
    || { echo "ERROR: failed to parse Log.final.out metrics for ${sample}" >&2; exit 1; }

  echo -e "${sample}\t${UM}\t${MM}\t${TS}\t${RD}" >> "${LOGFINAL_TSV}"

  BAMS+=("${BAM}")
done < "${MANIFEST}"

echo "BAMs found: ${#BAMS[@]}"
echo

# -------------------------
# Post-run sanity checks (alignment outliers)
# -------------------------
UNIQ_MIN=50
TOO_SHORT_MAX=20
MULTI_MAX=20
READS_MIN=1000000

OUTLIERS="${ALN_QC}/alignment_outliers.tsv"
SUMMARY="${ALN_QC}/alignment_summary.txt"
CATA="${ALN_QC}/alignment_catastrophic.tsv"

echo -e "sample\tuniq_map_pct\tmulti_map_pct\ttoo_short_pct\tn_reads\treasons" > "${OUTLIERS}"

awk -v umin="${UNIQ_MIN}" -v tmax="${TOO_SHORT_MAX}" -v mmax="${MULTI_MAX}" -v rmin="${READS_MIN}" -F'\t' '
NR==1 {next}
{
  s=$1; um=$3+0; mm=$4+0; ts=$5+0; rd=$6+0;
  r="";
  if (um < umin) r=r "low_unique;";
  if (ts > tmax) r=r "high_too_short;";
  if (mm > mmax) r=r "high_multimap;";
  if (rd < rmin) r=r "low_reads;";
  if (r != "") print s "\t" um "\t" mm "\t" ts "\t" rd "\t" r;
}' "${ALN_QC}/alignment_stats.tsv" >> "${OUTLIERS}"

awk -F'\t' '
function sort(a,n, i,j,t){for(i=1;i<=n;i++)for(j=i+1;j<=n;j++)if(a[i]>a[j]){t=a[i];a[i]=a[j];a[j]=t}}
NR==1{next}
{um[++n]=$3+0; mm[n]=$4+0; ts[n]=$5+0; rd[n]=$6+0}
END{
  sort(um,n); sort(mm,n); sort(ts,n); sort(rd,n);
  mid=int((n+1)/2);
  printf("Samples: %d\n", n);
  printf("UniqMap%%  min/med/max: %.2f / %.2f / %.2f\n", um[1], um[mid], um[n]);
  printf("MultiMap%% min/med/max: %.2f / %.2f / %.2f\n", mm[1], mm[mid], mm[n]);
  printf("TooShort%% min/med/max: %.2f / %.2f / %.2f\n", ts[1], ts[mid], ts[n]);
  printf("Reads     min/med/max: %.0f / %.0f / %.0f\n", rd[1], rd[mid], rd[n]);
}' "${ALN_QC}/alignment_stats.tsv" | tee "${SUMMARY}"

echo -e "sample\tuniq_map_pct\tmulti_map_pct\ttoo_short_pct\tn_reads\treasons" > "${CATA}"
awk -F'\t' '
NR==1{next}
{
  s=$1; um=$3+0; mm=$4+0; ts=$5+0; rd=$6+0; r="";
  if (um < 30) r=r "uniq<30;";
  if (ts > 40) r=r "too_short>40;";
  if (rd < 200000) r=r "reads<200k;";
  if (r!="") print s "\t" um "\t" mm "\t" ts "\t" rd "\t" r;
}' "${ALN_QC}/alignment_stats.tsv" >> "${CATA}"

CATA_COUNT=$(( $(wc -l < "${CATA}") - 1 ))
echo "Catastrophic flagged: ${CATA_COUNT} / ${NSAMPLES}"
echo

# -------------------------
# featureCounts
# -------------------------
echo "Running featureCounts..."

# Ensure temporary uncompressed GTF is removed even if job fails
trap 'rm -f "${GTF_TMP}"' EXIT

zcat "${GTF_GZ}" > "${GTF_TMP}"

featureCounts -T "${SLURM_CPUS_PER_TASK}" \
  -a "${GTF_TMP}" \
  -o "${FC_ROOT}/gene_counts.tsv" \
  -g gene_id \
  -t exon \
  "${BAMS[@]}"

echo
echo "DONE"
echo "Counts     : ${FC_ROOT}/gene_counts.tsv"
echo "STAR row QC : ${ALN_QC}/alignment_stats.tsv"
echo "Log.final QC: ${LOGFINAL_TSV}"
echo "Outliers   : ${OUTLIERS}"
echo "Summary    : ${SUMMARY}"
date
