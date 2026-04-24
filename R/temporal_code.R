## ============================================================
## TEMP_FIX_GSE119119_COUNT_COLUMN_NAMES.R
## Remove trailing _1 from GSE119119 count matrix sample columns
## ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(stringr)
})

counts_path <- "D:/IMRS_Project/03_counts/GSE119119/featurecounts/validation/gene_counts_clean.tsv"

counts <- read_tsv(counts_path, show_col_types = FALSE)

old_names <- colnames(counts)
new_names <- old_names

# Keep first column as Geneid/gene_id.
# Remove only trailing _1 from sample columns.
new_names[-1] <- str_replace(new_names[-1], "_1$", "")

cat("\nBefore:\n")
print(head(old_names, 20))

cat("\nAfter:\n")
print(head(new_names, 20))

# Safety checks
if (any(is.na(new_names)) || any(new_names == "")) {
  stop("Bad column names after cleaning.")
}

if (any(duplicated(new_names))) {
  print(new_names[duplicated(new_names)])
  stop("Duplicate column names after removing _1.")
}

colnames(counts) <- new_names

write_tsv(counts, counts_path)

cat("\nFixed count column names in:\n")
cat(counts_path, "\n")