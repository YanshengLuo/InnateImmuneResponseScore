#!/bin/bash
#SBATCH --job-name=fc_merge
#SBATCH --account=qsong1
#SBATCH --output=/orange/qsong1/Yansheng/logs/fc_merge_%j.out
#SBATCH --error=/orange/qsong1/Yansheng/logs/fc_merge_%j.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --mail-user=yanshengluo@ufl.edu
#SBATCH --mail-type=BEGIN,END,FAIL

set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "USAGE: sbatch 02c_merge_fc_batches.slurm <DATASET>" >&2
  exit 1
fi

DATASET="$1"
PROJECT_ROOT="/orange/qsong1/Yansheng"
FC_ROOT="${PROJECT_ROOT}/03_counts/${DATASET}/featurecounts"

module load python

export FC_ROOT

python - <<'PY'
import os
import glob
import pandas as pd

fc_root = os.environ["FC_ROOT"]
files = sorted(glob.glob(os.path.join(fc_root, "gene_counts.batch_*.tsv")))
if not files:
    raise SystemExit("ERROR: no batch count files found")

merged = None

for f in files:
    df = pd.read_csv(f, sep="\t", comment="#")

    meta = df.iloc[:, :6].copy()
    counts = df.iloc[:, 6:].copy()

    if merged is None:
        merged = pd.concat([meta, counts], axis=1)
    else:
        merged = merged.merge(
            pd.concat([meta[["Geneid"]], counts], axis=1),
            on="Geneid",
            how="inner"
        )

out = os.path.join(fc_root, "gene_counts.tsv")
merged.to_csv(out, sep="\t", index=False)
print(f"Wrote merged counts: {out}")
PY