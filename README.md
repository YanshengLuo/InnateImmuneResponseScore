# InnateImmuneResponseScore
This is a project contains python script for processing bulk-RNA data and R codes for data analysis and modeling.

HOW TO RUN: STAR alignment (array) + featureCounts merge
========================================================

This pipeline is designed for HiPerGator with a hard resource cap of:
- TOTAL CPU available to you at once: 20
- Memory available: up to 168 GB
- Account to charge: qsong1

It uses:
1) One-time STAR index build (mm39 / GRCm39 + GENCODE M38)
2) STAR alignment as a Slurm job array (parallelized, capped at 20 CPUs total)
3) featureCounts once (after all BAMs are produced)

----------------------------------------
Directory assumptions (your project tree)
----------------------------------------

PROJECT_ROOT:
  /orange/qsong1/Yansheng

Input FASTQs (example):
  /orange/qsong1/Yansheng/01_raw/fastq/GSE264344

Scripts live here:
  /orange/qsong1/Yansheng/Hypergator_scripts

Logs:
  /orange/qsong1/Yansheng/logs

Outputs per dataset (example for GSE264344):
  /orange/qsong1/Yansheng/03_counts/GSE264344/
    manifest/fastq_manifest.tsv
    star/<sample>/Aligned.sortedByCoord.out.bam
    star/<sample>/Aligned.sortedByCoord.out.bam.bai
    star/<sample>/Log.final.out
    star/<sample>/alignment_row.tsv
    qc_alignment/alignment_stats.tsv
    qc_alignment/alignment_outliers.tsv
    qc_alignment/alignment_catastrophic.tsv
    qc_alignment/alignment_summary.txt
    featurecounts/gene_counts.tsv
    featurecounts/gene_counts.tsv.summary

Shared reference + STAR index:
  /orange/qsong1/Yansheng/00_metadata/reference/mm39_gencodeM38/
    genome.fa.gz
    genes.gtf.gz
    STAR_index/

----------------------------------------
Prerequisites
----------------------------------------

1) Make sure these scripts exist:

  /orange/qsong1/Yansheng/Hypergator_scripts/00_prep_star_index.slurm
  /orange/qsong1/Yansheng/Hypergator_scripts/01_star_align_array.slurm
  /orange/qsong1/Yansheng/Hypergator_scripts/02_featurecounts_merge.slurm

2) FASTQs must be gzip-compressed and in ONE directory:
   <FASTQ_DIR>/*.fastq.gz

   Example:
   /orange/qsong1/Yansheng/01_raw/fastq/GSE264344/*.fastq.gz

----------------------------------------
How to run (copy/paste)
----------------------------------------

cd /orange/qsong1/Yansheng/Hypergator_scripts

# Set your dataset FASTQ directory here
FASTQ_DIR="/orange/qsong1/Yansheng/01_raw/fastq/GSE264344"

# (1) Build STAR index once (mm39 + GENCODE M38)
jid1=$(sbatch 00_prep_star_index.slurm | awk '{print $4}')
echo "prep_star_index jobid: ${jid1}"

# (2) STAR align as an array
# We cap total CPU at 20 by using:
#   cpus-per-task = 4
#   max concurrent tasks = 5
# => 4 * 5 = 20 CPUs total
#
# We set array as 0-999 to avoid needing to know N in advance.
# Extra tasks will self-exit if they exceed the manifest sample count.
jid2=$(sbatch --dependency=afterok:${jid1} --array=0-999%5 01_star_align_array.slurm "${FASTQ_DIR}" | awk '{print $4}')
echo "star_align array jobid: ${jid2}"

# (3) Run featureCounts after ALL array tasks finish successfully
jid3=$(sbatch --dependency=afterok:${jid2} 02_featurecounts_merge.slurm "${FASTQ_DIR}" | awk '{print $4}')
echo "featureCounts jobid: ${jid3}"

----------------------------------------
Monitoring
----------------------------------------

# Check queue status
squeue -u $USER

# Check output logs
ls -lh /orange/qsong1/Yansheng/logs | tail -n 20

# After alignment starts, check one sample output folder
DATASET=$(basename "${FASTQ_DIR}")
ls -lh "/orange/qsong1/Yansheng/03_counts/${DATASET}/star" | head

----------------------------------------
Expected final deliverable files
----------------------------------------

After the pipeline completes, confirm these exist:

DATASET=$(basename "${FASTQ_DIR}")

# Alignment stats table (all samples)
ls -lh "/orange/qsong1/Yansheng/03_counts/${DATASET}/qc_alignment/alignment_stats.tsv"

# featureCounts gene-level matrix
ls -lh "/orange/qsong1/Yansheng/03_counts/${DATASET}/featurecounts/gene_counts.tsv"
ls -lh "/orange/qsong1/Yansheng/03_counts/${DATASET}/featurecounts/gene_counts.tsv.summary"

----------------------------------------
Common failure modes
----------------------------------------

1) STAR index missing:
   - You skipped step (1) or it failed.
   - Fix: re-run 00_prep_star_index.slurm and check logs.

2) Manifest missing:
   - array task 0 failed before creating manifest
   - Fix: check star_align logs for task 0; re-run the array.

3) featureCounts missing BAM:
   - some STAR tasks failed; featureCounts merge job will stop
   - Fix: inspect the failed array task logs, re-run the STAR array.

----------------------------------------
Notes on read layout (SE vs PE)
----------------------------------------

The alignment script uses filename heuristics:
- Paired-end if it finds matching *_1/*_2 or *R1*/*R2* files
- Otherwise treats as single-end

This is fine for typical SRA fastq-dump output for bulk RNA-seq.
