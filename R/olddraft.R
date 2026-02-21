suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(purrr)
  library(tidyr)
})

project_root <- "C:/Users/john/Desktop/IMRS_Project"

# INPUT: where DE result TSVs are
de_out_root  <- file.path(project_root, "04_de", "comparison")

# OUTPUT: write ALL outputs here (per your request)
score_root   <- file.path(project_root, "05_score")
dir.create(score_root, showWarnings = FALSE, recursive = TRUE)

# -------------------------
# Parameters (tune if needed)
# -------------------------
padj_cutoff <- 0.05
lfc_cutoff  <- 1.0          # abs(log2FC) >= 1
support_cutoff <- 0.5       # core gene must pass in >=50% anchors
use_stat <- "median"        # "median" recommended; "mean" also ok

# anchor/calibration time rules
is_anchor_time <- function(t) !is.na(t) & t < 24
is_calib_time  <- function(t) !is.na(t) & t >= 24 & t <= 72

# -------------------------
# Helper: read one DE file (Option A)
#   - Keep ONLY the columns needed for modeling
#   - Drop any meta columns that can collide with meta_use columns
# -------------------------
read_de <- function(path) {
  df <- read_tsv(path, show_col_types = FALSE)

  # Harmonize gene column name
  if (!("gene_id" %in% names(df))) {
    if ("gene" %in% names(df)) {
      df <- df %>% rename(gene_id = gene)
    } else {
      names(df)[1] <- "gene_id"
    }
  }

  # Keep only columns used downstream (prevents duplicate name collisions)
  keep <- intersect(
    c("gene_id", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"),
    names(df)
  )
  df <- df %>% select(all_of(keep))

  df %>%
    mutate(
      log2FoldChange = as.numeric(log2FoldChange),
      padj = as.numeric(padj)
    )
}

# -------------------------
# Helper: parse metadata from filename
#   expects DE filename includes "...__H=24__..." etc.
# -------------------------
parse_meta_from_filename <- function(path) {
  bn <- basename(path)

  dataset_id <- str_extract(bn, "^GSE[^_]+(?:_[A-Za-z0-9]+)?")  # catches GSE190850_HUMAN too
  if (is.na(dataset_id)) dataset_id <- str_extract(bn, "^GSE\\d+")

  # time in hours: __H=24__
  time_h <- suppressWarnings(as.numeric(str_match(bn, "__H=([0-9.]+)__")[,2]))

  tissue <- str_match(bn, "__T=([^_]+)__")[,2]
  if (is.na(tissue)) tissue <- NA_character_

  # group label: __G=delivery_lnp__
  group <- str_match(bn, "__G=([^_]+)__")[,2]
  if (is.na(group)) group <- NA_character_

  # vs label: __VS=baseline_0
  vs <- str_match(bn, "__VS=([^\\.]+)")[,2]
  if (is.na(vs)) vs <- NA_character_

  contrast_label <- ifelse(!is.na(group) & !is.na(vs),
                           paste0(group, "__vs__", vs),
                           bn)

  tibble(
    dataset_id = dataset_id,
    tissue = tissue,
    time_h = time_h,
    contrast_label = contrast_label,
    file = path
  )
}

# -------------------------
# Collect all DE files
# -------------------------
de_files <- list.files(de_out_root, pattern = "__DE\\.tsv$", recursive = TRUE, full.names = TRUE)

if (length(de_files) == 0) {
  stop("No DE files found under: ", de_out_root, "\nExpected files ending with __DE.tsv")
}

meta <- map_dfr(de_files, parse_meta_from_filename) %>%
  mutate(
    time_bin = case_when(
      is_anchor_time(time_h) ~ "anchor_lt24h",
      is_calib_time(time_h)  ~ "calib_24to72h",
      TRUE ~ "other_or_unknown"
    ),
    # IMPORTANT: flag HUMAN datasets so they never enter anchor/calibration modeling
    is_human = str_detect(dataset_id, "_HUMAN$")
  )

# -------------------------
# Summary: contrast counts by dataset (includes HUMAN flag)
# NOTE: this is contrast-level counts (DE files), not raw sample counts.
# -------------------------
meta_summary_by_dataset <- meta %>%
  mutate(
    use_flag = case_when(
      time_bin == "anchor_lt24h" ~ "anchor",
      time_bin == "calib_24to72h" ~ "calibration",
      TRUE ~ "excluded"
    )
  ) %>%
  group_by(dataset_id, is_human) %>%
  summarise(
    total_contrasts = n(),
    anchor_contrasts = sum(use_flag == "anchor"),
    calibration_contrasts = sum(use_flag == "calibration"),
    excluded_contrasts = sum(use_flag == "excluded"),
    excluded_missing_time = sum(is.na(time_h)),
    excluded_gt72h = sum(!is.na(time_h) & time_h > 72),
    .groups = "drop"
  ) %>%
  arrange(desc(total_contrasts))

meta_summary_overall <- meta %>%
  mutate(
    use_flag = case_when(
      time_bin == "anchor_lt24h" ~ "anchor",
      time_bin == "calib_24to72h" ~ "calibration",
      TRUE ~ "excluded"
    )
  ) %>%
  summarise(
    dataset_n = n_distinct(dataset_id),
    dataset_n_human = n_distinct(dataset_id[is_human]),
    dataset_n_mouse = n_distinct(dataset_id[!is_human]),
    total_contrasts = n(),
    anchor_contrasts = sum(use_flag == "anchor"),
    calibration_contrasts = sum(use_flag == "calibration"),
    excluded_contrasts = sum(use_flag == "excluded"),
    excluded_missing_time = sum(is.na(time_h)),
    excluded_gt72h = sum(!is.na(time_h) & time_h > 72)
  )

message("\n=== DE file / contrast counts by dataset (HUMAN flagged) ===")
print(meta_summary_by_dataset, n = Inf)
message("\n=== Overall totals ===")
print(meta_summary_overall)

write_tsv(meta_summary_by_dataset, file.path(score_root, "contrast_counts_by_dataset.tsv"))
write_tsv(meta_summary_overall, file.path(score_root, "contrast_counts_overall.tsv"))

# -------------------------
# Build model using MOUSE ONLY (anchors + calibration time window)
# -------------------------
meta_use_mouse <- meta %>%
  filter(
    time_bin %in% c("anchor_lt24h", "calib_24to72h"),
    !is_human
  )

# Load mouse DE tables
all_de_mouse <- meta_use_mouse %>%
  mutate(de = map(file, read_de)) %>%
  select(-file) %>%
  tidyr::unnest(de)

# -------------------------
# Step 6B: build core gene set from MOUSE anchors only
# -------------------------
anchors_mouse <- all_de_mouse %>% filter(time_bin == "anchor_lt24h")

if (nrow(anchors_mouse) == 0) {
  stop("No MOUSE anchor contrasts found (time_h < 24). Check filename parsing for __H=..__")
}

anchors_flagged <- anchors_mouse %>%
  mutate(pass = !is.na(padj) & padj <= padj_cutoff &
                !is.na(log2FoldChange) & abs(log2FoldChange) >= lfc_cutoff)

n_anchor_contrasts <- n_distinct(anchors_flagged %>% select(dataset_id, tissue, time_h, contrast_label))
message("MOUSE anchor contrasts used: ", n_anchor_contrasts)

gene_support <- anchors_flagged %>%
  group_by(gene_id) %>%
  summarise(
    support = mean(pass, na.rm = TRUE),
    n_anchor_contrasts = n_distinct(paste(dataset_id, tissue, time_h, contrast_label, sep="|")),
    pos_frac = ifelse(sum(pass, na.rm=TRUE) > 0,
                      mean(log2FoldChange[pass] > 0, na.rm = TRUE),
                      NA_real_),
    .groups = "drop"
  ) %>%
  mutate(
    direction = case_when(
      !is.na(pos_frac) & pos_frac >= 0.8 ~ "UP",
      !is.na(pos_frac) & pos_frac <= 0.2 ~ "DOWN",
      TRUE ~ "MIXED"
    )
  )

core <- gene_support %>%
  filter(support >= support_cutoff) %>%
  arrange(desc(support))

write_tsv(core, file.path(score_root, "core_gene_set.tsv"))
message("Wrote core gene set: ", file.path(score_root, "core_gene_set.tsv"), " (n=", nrow(core), ")")

# -------------------------
# Step 6C: gene weights from MOUSE anchors (robust median by default)
# -------------------------
core_genes <- core$gene_id
anchor_core <- anchors_mouse %>% filter(gene_id %in% core_genes)

weights <- anchor_core %>%
  group_by(gene_id) %>%
  summarise(
    weight = if (use_stat == "median") median(log2FoldChange, na.rm = TRUE) else mean(log2FoldChange, na.rm = TRUE),
    mean_lfc = mean(log2FoldChange, na.rm = TRUE),
    median_lfc = median(log2FoldChange, na.rm = TRUE),
    sd_lfc = sd(log2FoldChange, na.rm = TRUE),
    n = sum(!is.na(log2FoldChange)),
    .groups = "drop"
  )

write_tsv(weights, file.path(score_root, "gene_weights.tsv"))
message("Wrote gene weights: ", file.path(score_root, "gene_weights.tsv"))

# -------------------------
# Step 6D: Score MOUSE anchors + calibrations
# -------------------------
w <- weights %>% select(gene_id, weight)

mouse_scored <- all_de_mouse %>%
  filter(time_bin %in% c("anchor_lt24h", "calib_24to72h")) %>%
  inner_join(w, by = "gene_id") %>%
  mutate(contrib = weight * log2FoldChange)

mouse_scores <- mouse_scored %>%
  group_by(dataset_id, tissue, time_h, time_bin, contrast_label) %>%
  summarise(
    imrs_score = sum(contrib, na.rm = TRUE),
    imrs_score_norm = sum(contrib, na.rm = TRUE) / sum(abs(weight), na.rm = TRUE),
    n_core_genes_used = sum(!is.na(log2FoldChange)),
    .groups = "drop"
  ) %>%
  arrange(dataset_id, time_h, tissue, contrast_label)

write_tsv(mouse_scores, file.path(score_root, "imrs_scores_mouse_anchor_and_calibration.tsv"))
message("Wrote MOUSE scores: ", file.path(score_root, "imrs_scores_mouse_anchor_and_calibration.tsv"))

mouse_calib_summary <- mouse_scores %>%
  filter(time_bin == "calib_24to72h") %>%
  group_by(dataset_id, contrast_label) %>%
  summarise(
    n = n(),
    score_mean = mean(imrs_score_norm, na.rm = TRUE),
    score_sd   = sd(imrs_score_norm, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(score_mean))

write_tsv(mouse_calib_summary, file.path(score_root, "calibration_summary_mouse.tsv"))
message("Wrote MOUSE calibration summary: ", file.path(score_root, "calibration_summary_mouse.tsv"))

# -------------------------
# OPTIONAL: Score HUMAN datasets separately (no influence on model)
# -------------------------
meta_use_human <- meta %>%
  filter(
    time_bin %in% c("anchor_lt24h", "calib_24to72h"),
    is_human
  )

if (nrow(meta_use_human) > 0) {

  all_de_human <- meta_use_human %>%
    mutate(de = map(file, read_de)) %>%
    select(-file) %>%
    tidyr::unnest(de)

  human_scored <- all_de_human %>%
    inner_join(w, by = "gene_id") %>%
    mutate(contrib = weight * log2FoldChange)

  human_scores <- human_scored %>%
    group_by(dataset_id, tissue, time_h, time_bin, contrast_label) %>%
    summarise(
      imrs_score = sum(contrib, na.rm = TRUE),
      imrs_score_norm = sum(contrib, na.rm = TRUE) / sum(abs(weight), na.rm = TRUE),
      n_core_genes_used = sum(!is.na(log2FoldChange)),
      .groups = "drop"
    ) %>%
    arrange(dataset_id, time_h, tissue, contrast_label)

  write_tsv(human_scores, file.path(score_root, "imrs_scores_human.tsv"))
  message("Wrote HUMAN scores: ", file.path(score_root, "imrs_scores_human.tsv"))

  human_calib_summary <- human_scores %>%
    filter(time_bin == "calib_24to72h") %>%
    group_by(dataset_id, contrast_label) %>%
    summarise(
      n = n(),
      score_mean = mean(imrs_score_norm, na.rm = TRUE),
      score_sd   = sd(imrs_score_norm, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(score_mean))

  write_tsv(human_calib_summary, file.path(score_root, "calibration_summary_human.tsv"))
  message("Wrote HUMAN calibration summary: ", file.path(score_root, "calibration_summary_human.tsv"))

} else {
  message("No HUMAN datasets detected for scoring (dataset_id ending with _HUMAN).")
}

message("\nDONE: mouse-only model built; human scored separately (no leakage).")