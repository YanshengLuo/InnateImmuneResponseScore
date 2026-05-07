#!/usr/bin/env Rscript
# ============================================================
# TEMPORARY FIX — GSE119119-only split-design generator
#
# Purpose:
#   Regenerate ONLY GSE119119 split-design files with genotype-matched controls.
#
# Fix:
#   WT delivery groups        vs WT_uninfected_whole_liver
#   TRIM21_KO delivery groups vs TRIM21_KO_uninfected_whole_liver
#
# This script does NOT touch any other dataset folder.
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
})

# ============================
# USER SETTINGS
# ============================
in_file <- "D:/IMRS_Project/00_metadata/verified_metadata/GSE119119_design.tsv"

out_root <- "D:/IMRS_Project/00_metadata/verified_metadata/splited"
dataset_folder <- file.path(out_root, "GSE119119_design")

# Set TRUE to remove old mixed 4/8 GSE119119 splits before writing corrected 4/4 splits
CLEAN_OLD_GSE119119_SPLITS <- TRUE

# GSE119119 has matched controls at 4 h, so strict is safest.
# Kept here only for compatibility.
CONTROL_FALLBACK_MODE <- "strict"  # "strict" or "time0"

# ----------------------------
# Helpers
# ----------------------------
safe_token <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- "NA"
  x <- str_trim(x)
  x <- str_replace_all(x, "[\\\\/:*?\"<>|]+", "_")
  x <- str_replace_all(x, "\\s+", "_")
  x <- str_replace_all(x, "_+", "_")
  x <- str_replace_all(x, "^_+|_+$", "")
  ifelse(nchar(x) == 0, "NA", x)
}

short_hash <- function(x) {
  ints <- utf8ToInt(as.character(x))
  if (length(ints) == 0) return("00000000")
  h <- sum(ints * seq_along(ints)) %% 100000000
  sprintf("%08d", h)
}

short_group_token <- function(g, max_len = 35) {
  g2 <- safe_token(g)
  if (nchar(g2) <= max_len) return(g2)
  paste0("GID", short_hash(g2))
}

to_hours <- function(x) {
  x <- as.character(x)
  x <- str_replace(x, "h$", "")
  suppressWarnings(as.numeric(x))
}

infer_genotype <- function(group_raw) {
  g <- as.character(group_raw)

  case_when(
    str_detect(g, "^TRIM21_KO_") ~ "TRIM21_KO",
    str_detect(g, "^WT_") ~ "WT",
    TRUE ~ "NA"
  )
}

is_gse119119_control <- function(group_raw) {
  group_raw %in% c(
    "WT_uninfected_whole_liver",
    "TRIM21_KO_uninfected_whole_liver"
  )
}

# ============================
# SAFETY CHECKS
# ============================
if (!file.exists(in_file)) {
  stop("Missing input design file: ", in_file)
}

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

if (CLEAN_OLD_GSE119119_SPLITS && dir.exists(dataset_folder)) {
  message("Removing old GSE119119 split folder only: ", dataset_folder)
  unlink(dataset_folder, recursive = TRUE, force = TRUE)
}

dir.create(dataset_folder, recursive = TRUE, showWarnings = FALSE)

# ============================
# LOAD DESIGN
# ============================
df <- read_tsv(in_file, show_col_types = FALSE)

req <- c("sample_id", "group", "tissue", "timepoint_hr", "batch")
miss <- setdiff(req, names(df))

if (length(miss) > 0) {
  stop("GSE119119 design missing columns: ", paste(miss, collapse = ", "))
}

df2 <- df %>%
  mutate(
    sample_id = as.character(sample_id),
    group_raw = as.character(group),
    tissue_raw = as.character(tissue),

    time_raw = as.character(timepoint_hr),
    time_h = to_hours(time_raw),
    time_label = ifelse(!is.na(time_h), paste0(time_h, "h"), safe_token(time_raw)),

    batch_raw = as.character(batch),
    batch_raw = ifelse(is.na(batch_raw) | batch_raw == "", "NA", batch_raw),

    genotype_raw = infer_genotype(group_raw),
    is_control = is_gse119119_control(group_raw)
  )

if (any(df2$genotype_raw == "NA")) {
  warning("Some GSE119119 groups did not parse into WT or TRIM21_KO:")
  print(unique(df2$group_raw[df2$genotype_raw == "NA"]))
}

message("Detected controls:")
print(df2 %>% filter(is_control) %>% count(genotype_raw, group_raw))

written <- character(0)
group_map_rows <- list()

# ============================
# SPLIT BY tissue + time + batch + genotype
# ============================
strata <- df2 %>%
  distinct(tissue_raw, time_label, batch_raw, genotype_raw) %>%
  arrange(tissue_raw, time_label, batch_raw, genotype_raw)

for (s in seq_len(nrow(strata))) {
  tissue_s <- strata$tissue_raw[s]
  time_s <- strata$time_label[s]
  batch_s <- strata$batch_raw[s]
  genotype_s <- strata$genotype_raw[s]

  block <- df2 %>%
    filter(
      tissue_raw == tissue_s,
      time_label == time_s,
      batch_raw == batch_s,
      genotype_raw == genotype_s
    )

  ctrl_block <- block %>%
    filter(is_control)

  control_time_policy <- "matched_time"

  if (nrow(ctrl_block) == 0 && CONTROL_FALLBACK_MODE == "time0") {
    ctrl_block <- df2 %>%
      filter(
        tissue_raw == tissue_s,
        batch_raw == batch_s,
        genotype_raw == genotype_s,
        is_control,
        !is.na(time_h),
        time_h == 0
      )

    if (nrow(ctrl_block) > 0) {
      control_time_policy <- "fallback_time0"
    }
  }

  if (nrow(ctrl_block) == 0) {
    message("No controls for stratum: ",
            tissue_s, " | ", time_s, " | ", batch_s, " | ", genotype_s)
    next
  }

  delivery_groups <- block %>%
    filter(!is_control) %>%
    pull(group_raw) %>%
    unique() %>%
    sort()

  if (length(delivery_groups) == 0) next

  for (g in delivery_groups) {
    deliv_block <- block %>%
      filter(!is_control, group_raw == g)

    if (nrow(deliv_block) == 0) next

    actual_control_label <- paste(unique(ctrl_block$group_raw), collapse = ";")

    sub <- bind_rows(ctrl_block, deliv_block) %>%
      mutate(
        condition_simple = ifelse(is_control, "CONTROL", "DELIVERY"),
        control_label = actual_control_label,
        contrast_label = paste0(g, "__vs__", actual_control_label),
        control_time_policy = control_time_policy,
        genotype_split = genotype_s
      )

    tissue_tok <- safe_token(tissue_s)

    time_tok_num <- to_hours(time_s)
    time_tok <- ifelse(!is.na(time_tok_num), as.character(time_tok_num), safe_token(time_s))

    batch_tok <- safe_token(batch_s)
    g_tok <- short_group_token(g)
    ctrl_tok <- safe_token(actual_control_label)

    group_map_rows[[length(group_map_rows) + 1]] <- tibble(
      tissue_raw = tissue_s,
      time_label = time_s,
      batch_raw = batch_s,
      genotype_raw = genotype_s,
      group_raw = g,
      group_token = g_tok,
      control_raw = actual_control_label,
      n_control = nrow(ctrl_block),
      n_delivery = nrow(deliv_block)
    )

    out_name <- sprintf(
      "GSE119119_design__T=%s__H=%s__B=%s__G=%s__VS=%s.tsv",
      tissue_tok,
      time_tok,
      batch_tok,
      g_tok,
      ctrl_tok
    )

    out_path <- file.path(dataset_folder, out_name)

    write_tsv(sub, out_path)
    written <- c(written, out_path)
  }
}

# ============================
# WRITE GSE119119-ONLY INDEX FILES
# ============================
if (length(written) == 0) {
  stop("No GSE119119 split files were written.")
}

group_map <- bind_rows(group_map_rows) %>%
  distinct() %>%
  arrange(genotype_raw, group_raw)

index <- tibble(design_file = written) %>%
  mutate(
    file_base = basename(design_file),
    tissue = str_match(file_base, "__T=([^_]+)__H=")[, 2],
    time_h = str_match(file_base, "__H=([^_]+)__B=")[, 2],
    batch = str_match(file_base, "__B=([^_]+)__G=")[, 2],
    group_token = str_match(file_base, "__G=([^_]+)__VS=")[, 2],
    control_token = str_match(file_base, "__VS=([^\\.]+)\\.tsv$")[, 2]
  ) %>%
  left_join(
    group_map %>%
      select(group_token, group_raw, genotype_raw, control_raw, n_control, n_delivery) %>%
      distinct(),
    by = "group_token"
  )

idx_path <- file.path(dataset_folder, "GSE119119_design__contrast_index.tsv")
map_path <- file.path(dataset_folder, "GSE119119_design__group_map.tsv")

write_tsv(index, idx_path)
write_tsv(group_map, map_path)

message("\nDone. Only GSE119119 was processed.")
message("Wrote ", length(written), " corrected GSE119119 split files.")
message("Index: ", idx_path)
message("Group map: ", map_path)

message("\nCheck n_control / n_delivery:")
print(group_map %>% select(genotype_raw, group_raw, control_raw, n_delivery, n_control))