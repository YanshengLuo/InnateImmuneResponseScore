#!/bin/bash
#SBATCH --job-name=featurecounts
#SBATCH --account=qsong1
#SBATCH --output=/orange/qsong1/Yansheng/logs/featurecounts_%j.out
#SBATCH --error=/orange/qsong1/Yansheng/logs/featurecounts_%j.err
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
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

REF_ROOT="${PROJECT_ROOT}/00_metadata/reference/mm39_gencodeM38"
GTF_GZ="${REF_ROOT}/genes.gtf.gz"

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

while IFS=$'\t' read -r sample fq1 fq2; do
  [[ "${sample}" == "sample_id" ]] && continue
  BAM="${STAR_ROOT}/${sample}/Aligned.sortedByCoord.out.bam"
  ROW="${STAR_ROOT}/${sample}/alignment_row.tsv"

  [[ -s "${BAM}" ]] || { echo "ERROR: missing BAM for ${sample}: ${BAM}" >&2; exit 1; }
  samtools quickcheck -v "${BAM}" || { echo "ERROR: BAM fails quickcheck: ${BAM}" >&2; exit 1; }

  [[ -s "${ROW}" ]] || { echo "ERROR: missing alignment row for ${sample}: ${ROW}" >&2; exit 1; }
  cat "${ROW}" >> "${ALN_QC}/alignment_stats.tsv"

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
zcat "${GTF_GZ}" > "${REF_ROOT}/genes.gtf"

featureCounts -T ${SLURM_CPUS_PER_TASK} \
  -a "${REF_ROOT}/genes.gtf" \
  -o "${FC_ROOT}/gene_counts.tsv" \
  -g gene_id \
  -t exon \
  "${BAMS[@]}"

rm -f "${REF_ROOT}/genes.gtf"

echo
echo "DONE"
echo "Counts: ${FC_ROOT}/gene_counts.tsv"
echo "Stats : ${ALN_QC}/alignment_stats.tsv"
echo "Outliers: ${OUTLIERS}"
echo "Summary : ${SUMMARY}"
date
