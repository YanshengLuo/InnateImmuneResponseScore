# InnateImmuneResponseScore
This is a project contains python script for processing bulk-RNA data and R codes for data analysis and modeling.

HOW TO RUN: STAR (array) + featureCounts (merge)
================================================

Goal:
- Build STAR index once (mm39/GRCm39 + GENCODE M38)
- Align samples in parallel using a Slurm array (capped to 20 CPUs total)
- Run featureCounts once after all alignments finish

You must set these 3 variables for your environment:
- ACCOUNT      : Slurm account to charge (e.g., qsong1)
- PROJECT_ROOT : where your project lives under /orange/<ACCOUNT>/...
- FASTQ_DIR    : directory that contains *.fastq.gz for one dataset (e.g., .../GSE264344)

Expected scripts (already created earlier):
- 00_prep_star_index.slurm
- 01_star_align_array.slurm
- 02_featurecounts_merge.slurm

------------------------------------------------
COPY/PASTE COMMANDS
------------------------------------------------

# Before start, you may need to run : dos2unix *.slurm
# 0) Set variables (EDIT THESE THREE)
ACCOUNT="qsong1"
PROJECT_ROOT="/orange/qsong1/PROJECT_NAME"      # <-- change PROJECT_NAME to your folder
SCRIPTS_DIR="$PROJECT_ROOT/Hypergator_scripts"
FASTQ_DIR="$PROJECT_ROOT/01_raw/fastq/GSE264344" # <-- change dataset folder as needed

cd "$SCRIPTS_DIR"

# 1) Build STAR index (run once)
```bash
sbatch --account="$ACCOUNT" 00_prep_star_index.slurm 
```

# 2) STAR alignment array
 CPU cap = 20 total CPUs:
   cpus-per-task = 4  (defined inside 01_star_align_array.slurm)
   concurrency   = 5  (the %5 below)
 => 4 * 5 = 20 CPUs total
jid2=$(sbatch --account="$ACCOUNT" --dependency=afterok:$jid1 --array=0-999%5 \
  01_star_align_array.slurm "$FASTQ_DIR" | awk '{print $4}')
echo "star_align array jobid = $jid2"

------------------------------------------------
or just copy:
```bash
sbatch --array=0-9%5 01_star_align_array.slurm /orange/qsong1/Yansheng/01_raw/fastq/GSE264344
```
------------------------------------------------

# 3) featureCounts merge (runs after ALL array tasks succeed)
```bash
jid3=$(sbatch --account="$ACCOUNT" --dependency=afterok:$jid2 \
  02_featurecounts_merge.slurm "$FASTQ_DIR" | awk '{print $4}')
echo "featureCounts jobid = $jid3"
```
------------------------------------------------
MONITORING
------------------------------------------------
```bash
squeue -u "$USER"
ls -lh "$PROJECT_ROOT/logs" | tail -n 20

DATASET=$(basename "$FASTQ_DIR")
echo "Outputs:"
echo "$PROJECT_ROOT/03_counts/$DATASET/"
```
------------------------------------------------
EXPECTED OUTPUTS (per dataset)
------------------------------------------------

Counts:
- $PROJECT_ROOT/03_counts/<DATASET>/featurecounts/gene_counts.tsv
- $PROJECT_ROOT/03_counts/<DATASET>/featurecounts/gene_counts.tsv.summary

Alignment QC tables:
- $PROJECT_ROOT/03_counts/<DATASET>/qc_alignment/alignment_stats.tsv
- $PROJECT_ROOT/03_counts/<DATASET>/qc_alignment/alignment_outliers.tsv
- $PROJECT_ROOT/03_counts/<DATASET>/qc_alignment/alignment_catastrophic.tsv
- $PROJECT_ROOT/03_counts/<DATASET>/qc_alignment/alignment_summary.txt

Per-sample STAR outputs:
- $PROJECT_ROOT/03_counts/<DATASET>/star/<SAMPLE>/Log.final.out
- $PROJECT_ROOT/03_counts/<DATASET>/star/<SAMPLE>/Aligned.sortedByCoord.out.bam
- $PROJECT_ROOT/03_counts/<DATASET>/star/<SAMPLE>/Aligned.sortedByCoord.out.bam.bai

------------------------------------------------
NOTES
------------------------------------------------

- FASTQ_DIR must contain .fastq.gz files directly under it (no nested subfolders).
- The alignment script detects PE vs SE by filename patterns (*_1/*_2 or *R1/*R2).
- If your cluster requires the account to be set ONLY in the script headers,
  then remove '--account=...' from the sbatch commands above and ensure each
  slurm script contains:
    #SBATCH --account=
