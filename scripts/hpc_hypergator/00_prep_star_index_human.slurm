#!/bin/bash
#SBATCH --job-name=prep_star_index_human
#SBATCH --account=qsong1
#SBATCH --output=/orange/qsong1/Yansheng/logs/prep_star_index_human_%j.out
#SBATCH --error=/orange/qsong1/Yansheng/logs/prep_star_index_human_%j.err
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --mail-user=yanshengluo@ufl.edu
#SBATCH --mail-type=BEGIN,END,FAIL

set -euo pipefail

PROJECT_ROOT="/orange/qsong1/Yansheng"
LOG_DIR="${PROJECT_ROOT}/logs"

# Human reference: GENCODE Human release 49 (GRCh38 / GRCh38.p14)
REF_ROOT="${PROJECT_ROOT}/00_metadata/reference/GRCh38_gencode49"
FA_GZ="${REF_ROOT}/genome.fa.gz"
GTF_GZ="${REF_ROOT}/genes.gtf.gz"
STAR_INDEX="${REF_ROOT}/STAR_index"

mkdir -p "${LOG_DIR}" "${REF_ROOT}" "${STAR_INDEX}"

module load star

echo "=== Prep STAR index (Human: GRCh38 + GENCODE v49) ==="
date
echo "REF_ROOT: ${REF_ROOT}"

# Pinned reference files (GENCODE Human release 49)
FA_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/GRCh38.primary_assembly.genome.fa.gz"
GTF_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.primary_assembly.annotation.gtf.gz"

if [[ ! -s "${FA_GZ}" || ! -s "${GTF_GZ}" ]]; then
  echo "Downloading reference FASTA+GTF (GENCODE Human v49 / GRCh38 primary assembly)"
  wget -O "${FA_GZ}"  "${FA_URL}"
  wget -O "${GTF_GZ}" "${GTF_URL}"
fi

if [[ -s "${STAR_INDEX}/Genome" ]]; then
  echo "STAR index already exists: ${STAR_INDEX}"
  exit 0
fi

# Build index. Use sjdbOverhang=100 as a robust default across datasets.
echo "Building STAR index..."
zcat "${FA_GZ}"  > "${REF_ROOT}/genome.fa"
zcat "${GTF_GZ}" > "${REF_ROOT}/genes.gtf"

STAR --runMode genomeGenerate \
     --runThreadN "${SLURM_CPUS_PER_TASK}" \
     --genomeDir "${STAR_INDEX}" \
     --genomeFastaFiles "${REF_ROOT}/genome.fa" \
     --sjdbGTFfile "${REF_ROOT}/genes.gtf" \
     --sjdbOverhang 100

rm -f "${REF_ROOT}/genome.fa" "${REF_ROOT}/genes.gtf"

echo "DONE: STAR index built at ${STAR_INDEX}"
date
