#!/bin/bash
#SBATCH --job-name=prep_star_index
#SBATCH --account=qsong1
#SBATCH --output=/orange/qsong1/Yansheng/logs/prep_star_index_%j.out
#SBATCH --error=/orange/qsong1/Yansheng/logs/prep_star_index_%j.err
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=80G
#SBATCH --mail-user=yanshengluo@ufl.edu
#SBATCH --mail-type=BEGIN,END,FAIL

set -euo pipefail

PROJECT_ROOT="/orange/qsong1/Yansheng"
LOG_DIR="${PROJECT_ROOT}/logs"

REF_ROOT="${PROJECT_ROOT}/00_metadata/reference/mm39_gencodeM38"
FA_GZ="${REF_ROOT}/genome.fa.gz"
GTF_GZ="${REF_ROOT}/genes.gtf.gz"
STAR_INDEX="${REF_ROOT}/STAR_index"

mkdir -p "${LOG_DIR}" "${REF_ROOT}" "${STAR_INDEX}"

module load star

echo "=== Prep STAR index ==="
date
echo "REF_ROOT: ${REF_ROOT}"

# Pinned reference: GENCODE Mouse M38 (GRCm39 / mm39)
if [[ ! -s "${FA_GZ}" || ! -s "${GTF_GZ}" ]]; then
  echo "Downloading reference FASTA+GTF (GENCODE M38 / GRCm39)"
  wget -O "${FA_GZ}"  "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M38/GRCm39.primary_assembly.genome.fa.gz"
  wget -O "${GTF_GZ}" "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M38/gencode.vM38.primary_assembly.annotation.gtf.gz"
fi

if [[ -s "${STAR_INDEX}/Genome" ]]; then
  echo "STAR index already exists: ${STAR_INDEX}"
  exit 0
fi

# Build index. Use sjdbOverhang=100 as a robust default across datasets.
# (Overhang mainly affects splice junction DB; this is fine for mixed read lengths.)
echo "Building STAR index..."
zcat "${FA_GZ}" > "${REF_ROOT}/genome.fa"
zcat "${GTF_GZ}" > "${REF_ROOT}/genes.gtf"

STAR --runMode genomeGenerate \
     --runThreadN ${SLURM_CPUS_PER_TASK} \
     --genomeDir "${STAR_INDEX}" \
     --genomeFastaFiles "${REF_ROOT}/genome.fa" \
     --sjdbGTFfile "${REF_ROOT}/genes.gtf" \
     --sjdbOverhang 100

rm -f "${REF_ROOT}/genome.fa" "${REF_ROOT}/genes.gtf"

echo "DONE: STAR index built at ${STAR_INDEX}"
date
