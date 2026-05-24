#!/usr/bin/env python3
import argparse
import sys
import re
from pathlib import Path

import pandas as pd
import numpy as np


NON_GENE_PREFIXES = ("__",)  # featureCounts special rows typically start with "__"


def is_integer_series(s: pd.Series) -> bool:
    """
    Returns True if all non-NA values are integers (or integer-like floats).
    """
    # If dtype is integer already, good.
    if pd.api.types.is_integer_dtype(s):
        return True
    # If float, check integer-like
    if pd.api.types.is_float_dtype(s):
        vals = s.dropna().to_numpy()
        return np.all(np.isclose(vals, np.round(vals)))
    # If object/string, attempt numeric conversion
    try:
        vals = pd.to_numeric(s.dropna(), errors="raise").to_numpy()
        return np.all(np.isclose(vals, np.round(vals)))
    except Exception:
        return False


def main():
    ap = argparse.ArgumentParser(
        description="Validate and clean a gene count matrix (genes x samples)."
    )
    ap.add_argument("--counts", required=True, help="Input counts TSV/CSV (genes in rows, samples in columns).")
    ap.add_argument("--outdir", required=True, help="Output directory.")
    ap.add_argument("--sep", default="\t", help="Delimiter (default: tab). Use ',' for CSV.")
    ap.add_argument("--gene-col", default=None, help="Gene ID column name. Default: first column.")
    ap.add_argument("--drop-non-gene-rows", action="store_true",
                    help="Drop featureCounts meta rows like __no_feature.")
    ap.add_argument("--collapse-duplicates", action="store_true",
                    help="If duplicate gene IDs exist, collapse by summing counts.")
    ap.add_argument("--min-libsize", type=float, default=None,
                    help="Optional hard minimum library size. Samples below are flagged.")
    ap.add_argument("--outlier-factor", type=float, default=5.0,
                    help="Flag samples with libsize > factor*median or < median/factor (default: 5).")
    args = ap.parse_args()

    counts_path = Path(args.counts)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    report_lines = []
    report_lines.append(f"INPUT: {counts_path}")
    report_lines.append(f"SEP: {repr(args.sep)}")

    # Read
    try:
        df = pd.read_csv(counts_path, sep=args.sep, dtype=str)  # read as str first to validate
    except Exception as e:
        print(f"ERROR: failed to read counts file: {e}", file=sys.stderr)
        sys.exit(2)

    if df.shape[1] < 2:
        print("ERROR: counts table must have at least 2 columns (gene + >=1 sample).", file=sys.stderr)
        sys.exit(3)

    # Determine gene column
    gene_col = args.gene_col or df.columns[0]
    if gene_col not in df.columns:
        print(f"ERROR: gene column '{gene_col}' not found. Available: {list(df.columns)}", file=sys.stderr)
        sys.exit(4)

    # Basic shape
    report_lines.append(f"RAW_SHAPE: rows={df.shape[0]} cols={df.shape[1]}")
    report_lines.append(f"GENE_COL: {gene_col}")

    # Remove empty gene IDs
    n_empty_gene = (df[gene_col].isna() | (df[gene_col].astype(str).str.strip() == "")).sum()
    if n_empty_gene > 0:
        report_lines.append(f"WARNING: empty gene IDs found: {n_empty_gene}. Dropping them.")
        df = df[~(df[gene_col].isna() | (df[gene_col].astype(str).str.strip() == ""))].copy()

    # Drop featureCounts meta rows
    if args.drop_non_gene_rows:
        is_meta = df[gene_col].astype(str).str.startswith(NON_GENE_PREFIXES)
        meta_count = int(is_meta.sum())
        report_lines.append(f"DROP_META_ROWS: {meta_count}")
        df = df[~is_meta].copy()

    # Convert sample columns to numeric
    sample_cols = [c for c in df.columns if c != gene_col]
    report_lines.append(f"N_SAMPLES: {len(sample_cols)}")

    # Check sample column uniqueness
    if len(sample_cols) != len(set(sample_cols)):
        report_lines.append("ERROR: duplicate sample column names detected.")
        (outdir / "qc_report.txt").write_text("\n".join(report_lines) + "\n")
        sys.exit(5)

    # Try numeric conversion
    numeric = df[sample_cols].apply(pd.to_numeric, errors="coerce")

    # Identify non-numeric cells
    non_numeric_cells = int(numeric.isna().sum().sum()) - int(df[sample_cols].isna().sum().sum())
    # The subtraction prevents double-counting true NAs.
    report_lines.append(f"NON_NUMERIC_CELLS_COERCED_TO_NA: {non_numeric_cells}")

    # Check integer-like counts
    # We'll evaluate on the full matrix
    all_integer_like = True
    for c in sample_cols:
        if not is_integer_series(numeric[c]):
            all_integer_like = False
            report_lines.append(f"WARNING: column '{c}' has non-integer-like values (possible TPM/FPKM/log).")
    report_lines.append(f"ALL_INTEGER_LIKE: {all_integer_like}")

    # Replace df sample cols with numeric
    df_numeric = pd.concat([df[[gene_col]].reset_index(drop=True), numeric.reset_index(drop=True)], axis=1)

    # Duplicate genes
    dup_mask = df_numeric[gene_col].duplicated(keep=False)
    n_dup_rows = int(dup_mask.sum())
    n_dup_genes = int(df_numeric.loc[dup_mask, gene_col].nunique()) if n_dup_rows > 0 else 0
    report_lines.append(f"DUPLICATE_GENE_ROWS: {n_dup_rows}")
    report_lines.append(f"DUPLICATE_GENE_IDS: {n_dup_genes}")

    if n_dup_rows > 0 and args.collapse_duplicates:
        report_lines.append("ACTION: collapsing duplicate gene IDs by summing counts.")
        df_numeric = (
            df_numeric
            .groupby(gene_col, as_index=False)[sample_cols]
            .sum()
        )

    # Library sizes
    libsize = df_numeric[sample_cols].sum(axis=0, numeric_only=True)
    median_lib = float(np.median(libsize.values))
    report_lines.append(f"MEDIAN_LIBSIZE: {median_lib:.3f}")

    factor = float(args.outlier_factor)
    low_thr = median_lib / factor if median_lib > 0 else 0
    high_thr = median_lib * factor
    flags = []
    for s, v in libsize.items():
        flag_reasons = []
        if v <= 0:
            flag_reasons.append("ZERO_OR_NEGATIVE_LIBSIZE")
        if median_lib > 0 and (v < low_thr):
            flag_reasons.append(f"LOW_OUTLIER(<median/{factor})")
        if median_lib > 0 and (v > high_thr):
            flag_reasons.append(f"HIGH_OUTLIER(>{factor}*median)")
        if args.min_libsize is not None and v < args.min_libsize:
            flag_reasons.append(f"BELOW_MIN_LIBSIZE({args.min_libsize})")
        flags.append(";".join(flag_reasons) if flag_reasons else "OK")

    lib_df = pd.DataFrame({
        "sample": libsize.index,
        "library_size": libsize.values,
        "flag": flags
    }).sort_values("library_size")

    n_flagged = int((lib_df["flag"] != "OK").sum())
    report_lines.append(f"FLAGGED_SAMPLES: {n_flagged}")

    # Write outputs
    cleaned_path = outdir / "gene_counts_clean.tsv"
    lib_path = outdir / "library_sizes.tsv"
    report_path = outdir / "qc_report.txt"

    df_numeric.to_csv(cleaned_path, sep="\t", index=False)
    lib_df.to_csv(lib_path, sep="\t", index=False)
    report_path.write_text("\n".join(report_lines) + "\n")

    print(f"[OK] Wrote: {cleaned_path}")
    print(f"[OK] Wrote: {lib_path}")
    print(f"[OK] Wrote: {report_path}")

    # Exit code: if severe issues
    # If not integer-like, return non-zero so pipelines can stop if desired.
    if not all_integer_like:
        print("[WARN] Counts are not all integer-like; verify this is raw counts.", file=sys.stderr)
        sys.exit(10)

    sys.exit(0)


if __name__ == "__main__":
    main()
