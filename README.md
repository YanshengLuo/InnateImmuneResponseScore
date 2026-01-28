# Innate Immune Response Score
This is a project contains python and shell script for processing bulk-RNA data and R codes for data analysis and modeling.

# Step 1: Download cDNA files using srr_download script
```bash
sbatch srr_download_dabaseinput.slurm <database>[start end]
```
# Step 2: Run sequence alignment and featurecounts to generate gene expression table

### HOW TO RUN: STAR (array) + featureCounts (merge)

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

# Before start, you may need to run : 
```bash
dos2unix *.slurm
```
# 0) Set variables (EDIT THESE THREE)
```bash
ACCOUNT="qsong1"
PROJECT_ROOT="/orange/qsong1/PROJECT_NAME"      # <-- change PROJECT_NAME to your folder
SCRIPTS_DIR="$PROJECT_ROOT/Hypergator_scripts"
FASTQ_DIR="$PROJECT_ROOT/01_raw/fastq/GSE264344" # <-- change dataset folder as needed
cd "$SCRIPTS_DIR"
```
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


# QC STEPS (MANDATORY) — COPY / PASTE SECTION
 This pipeline uses three explicit QC stages.
 Each QC step targets a different failure mode and must be run
 in order before downstream modeling.

 QC-0 : FASTQ-level QC (raw reads sanity check)
 QC-1 : Alignment-level QC (STAR output validation)
 QC-2 : Count-matrix QC (library size & matrix integrity)


# QC-0) FASTQ-level QC (qc_all.slurm)

# WHY:
   - Detect corrupted or incomplete FASTQ downloads
   - Check read length consistency and gross quality issues
  - Ensure raw reads are valid BEFORE alignment
  - Diagnostic only (no trimming or filtering)

# INPUT:
```bash
FASTQ_DIR containing *.fastq.gz for ONE dataset
```
# OUTPUT:
```bash
   /orange/qsong1/Yansheng/02_qc/<DATASET>/fastqc/
   /orange/qsong1/Yansheng/02_qc/<DATASET>/multiqc/multiqc_report.html
   /orange/qsong1/Yansheng/02_qc/<DATASET>/summary/

PROJECT_ROOT="/orange/qsong1/Yansheng"
FASTQ_DIR="$PROJECT_ROOT/01_raw/fastq/GSE264344"
```
```bash
cd "$PROJECT_ROOT/Hypergator_scripts/InnateImmuneResponseScore"

sbatch qc_all.slurm "$FASTQ_DIR"
```


# QC-1) Alignment-level QC (02_featurecounts_merge.slurm)

 WHY:
   - Verify STAR alignment completed successfully
   - Detect samples with pathological mapping rates
   - Catch silent alignment failures BEFORE quantification
   - Prevent broken samples from contaminating gene counts
------------------------------------------------------------

# INPUT:
```bash
   FASTQ_DIR (same as above)
```
```bash
# OUTPUT:
   $PROJECT_ROOT/03_counts/<DATASET>/qc_alignment/alignment_summary.txt
   $PROJECT_ROOT/03_counts/<DATASET>/qc_alignment/alignment_outliers.tsv
   $PROJECT_ROOT/03_counts/<DATASET>/qc_alignment/alignment_catastrophic.tsv
   $PROJECT_ROOT/03_counts/<DATASET>/featurecounts/gene_counts.tsv
```
```bash   
sbatch 02_featurecounts_merge.slurm "$FASTQ_DIR"
```


# QC-2) Count-matrix QC (10_validate_counts.slurm)

# WHY:
   - Validate integrity of gene-level count matrix
   - Detect extremely low-depth samples
   - Ensure matrix is safe for normalization and DE analysis
   - Enforce traceable sample exclusion rules

 INPUT:
 ```bash
   gene_counts.tsv generated by featureCounts
```
 OUTPUT:
 ```bash
   $OUTDIR/library_sizes.tsv
   $OUTDIR/qc_report.txt
   $OUTDIR/gene_counts_clean.tsv

DATASET=$(basename "$FASTQ_DIR")
COUNTS="$PROJECT_ROOT/03_counts/${DATASET}/featurecounts/gene_counts.tsv"
OUTDIR="$PROJECT_ROOT/03_counts/${DATASET}/featurecounts/validation"
```
```bash
sbatch 10_validate_counts.slurm "$COUNTS" "$OUTDIR"
```


# QC DECISION RULES (SUMMARY)

 - FASTQ QC (QC-0):
     Re-download only if reads are clearly broken or truncated

 - Alignment QC (QC-1):
    alignment_catastrophic.tsv  -> DROP sample
    alignment_outliers.tsv      -> FLAG for review

 - Count QC (QC-2):
     Extreme low-depth samples   -> DROP or reprocess

## All QC actions occur BEFORE normalization or modeling.
 #QC philosophy (important)

- QC is pre-analysis only

- No batch correction or biological tuning is performed

- Each QC step targets a different failure mode

- All exclusions are traceable and reviewer-defensible
============================================================

