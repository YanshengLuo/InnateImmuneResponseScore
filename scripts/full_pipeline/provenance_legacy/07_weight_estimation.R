#!/usr/bin/env Rscript
# ============================================================
# IMRS Step 07 — Gene Weight Estimation (ANCHOR ONLY, FROZEN)
#
# Inputs:
#   05_score/anchors/core_gene_set.tsv
#   05_score/anchors/gene_heterogeneity.tsv   (beta_meta, se_meta, I2)
#   05_score/anchors/gene_power.tsv           (power column)
#
# Cap rule:
#   c = 95th percentile of |beta_meta| (core genes only)
#   beta_capped = clamp(beta_meta, [-c, +c])
#
# Penalties:
#   penalty_het   = 1 - I2
#   penalty_power = 0.2 if power < 0.80 else 1.0
#   penalty_final = max(penalty_het * penalty_power, 0.2)
#
# Final weight:
#   weight = beta_capped * penalty_final
#
# Output:
#   05_score/anchors/gene_weights.tsv
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) args[1] else "D:/IMRS_Project"

phase_root <- file.path(project_root, "05_score", "anchors")

if (!dir.exists(phase_root)) {
  stop("Missing folder: ", phase_root)
}

# -------------------------
# Required files
# -------------------------
core_path <- file.path(phase_root, "core_gene_set.tsv")
het_path  <- file.path(phase_root, "gene_heterogeneity.tsv")
pow_path  <- file.path(phase_root, "gene_power.tsv")

if (!file.exists(core_path)) stop("Missing core_gene_set.tsv: ", core_path)
if (!file.exists(het_path))  stop("Missing gene_heterogeneity.tsv: ", het_path)
if (!file.exists(pow_path))  stop("Missing gene_power.tsv: ", pow_path)

# -------------------------
# Load inputs
# -------------------------
core <- read_tsv(core_path, show_col_types = FALSE)
het  <- read_tsv(het_path,  show_col_types = FALSE)
pow  <- read_tsv(pow_path,  show_col_types = FALSE)

if (!("gene_id" %in% names(core))) stop("core_gene_set.tsv must contain gene_id")
if (!all(c("gene_id","beta_meta","se_meta") %in% names(het)))
  stop("gene_heterogeneity.tsv must contain gene_id, beta_meta, se_meta")

# detect I2 column
i2_col <- if ("I2" %in% names(het)) "I2" else if ("i2" %in% names(het)) "i2" else NA_character_
if (is.na(i2_col)) stop("gene_heterogeneity.tsv missing I2 column")

# detect power column (tolerant)
power_col <- dplyr::case_when(
  "power" %in% names(pow)       ~ "power",
  "power_delta" %in% names(pow) ~ "power_delta",
  "power_delta1" %in% names(pow)~ "power_delta1",
  TRUE ~ NA_character_
)
if (is.na(power_col))
  stop("gene_power.tsv missing recognizable power column")

core_genes <- unique(as.character(core$gene_id))

df <- core %>%
  transmute(gene_id = as.character(gene_id)) %>%
  left_join(
    het %>%
      transmute(
        gene_id   = as.character(gene_id),
        beta_meta = as.numeric(beta_meta),
        se_meta   = as.numeric(se_meta),
        I2        = as.numeric(.data[[i2_col]])
      ),
    by = "gene_id"
  ) %>%
  left_join(
    pow %>%
      transmute(
        gene_id = as.character(gene_id),
        power   = as.numeric(.data[[power_col]])
      ),
    by = "gene_id"
  ) %>%
  filter(!is.na(beta_meta), !is.na(se_meta))

if (nrow(df) == 0) {
  stop("No valid rows after merging core + 6B + 6C.")
}

# -------------------------
# Cap rule
# -------------------------
cap_quantile <- 0.95
cap_c <- as.numeric(
  quantile(abs(df$beta_meta), probs = cap_quantile,
           na.rm = TRUE, names = FALSE)
)

if (!is.finite(cap_c) || cap_c <= 0) {
  stop("Invalid cap value computed.")
}

clamp <- function(x, lo, hi) pmax(pmin(x, hi), lo)

# -------------------------
# Penalties
# -------------------------
power_cutoff <- 0.80
penalty_power_low <- 0.2
penalty_floor <- 0.2

df <- df %>%
  mutate(
    beta_capped = clamp(beta_meta, -cap_c, cap_c),
    capped_flag = abs(beta_meta) > cap_c,

    penalty_het   = ifelse(is.na(I2), 1.0, pmax(0, 1 - I2)),
    penalty_power = ifelse(!is.na(power) & power < power_cutoff,
                           penalty_power_low, 1.0),

    penalty_combined = penalty_het * penalty_power,
    penalty_final    = pmax(penalty_combined, penalty_floor),

    weight = beta_capped * penalty_final
  ) %>%
  arrange(desc(abs(weight)))

# -------------------------
# Output
# -------------------------
out_path <- file.path(phase_root, "gene_weights.tsv")

out <- df %>%
  transmute(
    gene = gene_id,
    weight = weight,
    beta_meta = beta_meta,
    se_meta = se_meta,
    I2 = I2,
    power = power,
    penalty_het = penalty_het,
    penalty_power = penalty_power,
    penalty_final = penalty_final,
    beta_capped = beta_capped,
    cap_c = cap_c,
    capped_flag = capped_flag
  )

write_tsv(out, out_path)

cat("\n==============================\n")
cat("IMRS Step 07 — Anchor Weights Built\n")
cat("==============================\n")
cat("Core genes: ", nrow(out), "\n", sep = "")
cat("Cap (95th percentile): ", signif(cap_c, 4), "\n", sep = "")
cat("Output: ", out_path, "\n", sep = "")
cat("\nTop 10 by |weight|:\n")
print(head(out, 10))
cat("\nDONE\n\n")