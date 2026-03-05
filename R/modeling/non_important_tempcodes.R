

# ============================================================
# IMRS Core Gene Set Diagnostics
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

# ---- CHANGE PATH IF NEEDED ----
core_path <- "D:/IMRS_Project/05_score/anchors/core_gene_set.tsv"

if (!file.exists(core_path)) {
  stop("File not found: ", core_path)
}

core <- read_tsv(core_path, show_col_types = FALSE)

cat("\n==============================\n")
cat("CORE GENE SET DIAGNOSTICS\n")
cat("==============================\n\n")

# 1️⃣ Total size
cat("Total genes in core set:\n")
cat(nrow(core), "\n\n")

# 2️⃣ Datasets support distribution
if ("support_datasets" %in% names(core)) {
  cat("Support_datasets distribution:\n")
  print(table(core$support_datasets))
  cat("\n")
} else {
  cat("Column 'support_datasets' not found.\n\n")
}

# 3️⃣ Support fraction summary
if ("support_fraction" %in% names(core)) {
  cat("Support_fraction summary:\n")
  print(summary(core$support_fraction))
  cat("\n")
} else {
  cat("Column 'support_fraction' not found.\n\n")
}

# 4️⃣ Strong support genes (≥4 datasets)
if ("support_datasets" %in% names(core)) {
  strong_support <- sum(core$support_datasets >= 4, na.rm = TRUE)
  cat("Genes supported by ≥4 datasets:\n")
  cat(strong_support, "\n\n")
}

# 5️⃣ Preview top genes
cat("Top 20 genes:\n")
print(head(core, 20))

cat("\n==============================\n")
cat("DONE\n")
cat("==============================\n")

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) args[1] else "D:/IMRS_Project"

power_path <- file.path(project_root, "05_score", "anchors", "gene_power.tsv")

if (!file.exists(power_path)) {
  stop("File not found: ", power_path)
}

power <- read_tsv(power_path, show_col_types = FALSE)

cat("\n==============================\n")
cat("POWER DIAGNOSTICS (ANCHOR)\n")
cat("==============================\n\n")

total_genes <- nrow(power)
total_core  <- sum(power$is_core, na.rm = TRUE)
low_power_core <- sum(power$is_core & power$low_power_flag, na.rm = TRUE)

cat("Total genes:", total_genes, "\n")
cat("Total core genes:", total_core, "\n")
cat("Low-power core genes (<0.8):", low_power_core, "\n\n")

if (total_core > 0) {
  core_power <- power$power_delta[power$is_core]

  cat("Power summary (core genes):\n")
  print(summary(core_power))
  
  cat("\nMin power (core):", min(core_power, na.rm = TRUE), "\n")
  cat("Median power (core):", median(core_power, na.rm = TRUE), "\n")
  cat("Mean power (core):", mean(core_power, na.rm = TRUE), "\n")
}

cat("\n==============================\n")
cat("DONE\n")
cat("==============================\n\n")

