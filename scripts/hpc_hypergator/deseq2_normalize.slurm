#!/bin/bash
#SBATCH --job-name=deseq2_norm
#SBATCH --account=qsong1
#SBATCH --partition=hpg-default
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=/orange/qsong1/Yansheng/logs/deseq2_norm_%j.out
#SBATCH --error=/orange/qsong1/Yansheng/logs/deseq2_norm_%j.err

set -euo pipefail

DATASET_ID="${1:?Usage: sbatch 04_deseq2_normalize.slurm <GSE_ID>}"

PROJECT_ROOT="/orange/qsong1/Yansheng"
RSCRIPT="/orange/qsong1/Yansheng/Hypergator_scripts/InnateImmuneResponseScore/R/04_deseq2_normalize.R"

module purge
module load R/4.3

mkdir -p /orange/qsong1/Yansheng/logs

Rscript "${RSCRIPT}" "${DATASET_ID}" "${PROJECT_ROOT}"
