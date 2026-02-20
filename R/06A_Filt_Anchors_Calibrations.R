suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(purrr)
  library(tibble)
})

# ----------------------------
# Paths
# ----------------------------
project_root <- "C:/Users/john/Desktop/IMRS_Project"
comparison_root <- file.path(project_root, "04_de", "comparison")
de_root <- file.path(project_root, "04_de")

# ----------------------------
# 1) Find all DE result files
# ----------------------------
de_files <- list.files(
  comparison_root,
  pattern = "__DE\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(de_files) == 0) stop("No __DE.tsv files found under: ", comparison_root)

# ----------------------------
# Helpers
# ----------------------------
time_bin_from_h <- function(time_h) {
  dplyr::case_when(
    !is.na(time_h) & time_h < 24 ~ "anchor_lt24h",
    !is.na(time_h) & time_h >= 24 & time_h <= 72 ~ "calib_24to72h",
    !is.na(time_h) & time_h > 72 ~ "late_gt72h",
    TRUE ~ "time_unknown"
  )
}

parse_meta_from_filename <- function(path) {
  bn <- basename(path)
  
  # dataset_id is folder under comparison_root: .../comparison/<DATASET>/deseq2_contrasts/...
  dataset_id <- basename(dirname(dirname(path)))
  
  tissue <- str_match(bn, "__T=([^_]+)__")[,2]
  time_h <- suppressWarnings(as.numeric(str_match(bn, "__H=([0-9.]+)__")[,2]))
  group  <- str_match(bn, "__G=([^_]+)__")[,2]
  vs     <- str_match(bn, "__VS=([^_\\.]+)")[,2]
  
  contrast_label <- ifelse(!is.na(group) & !is.na(vs),
                           paste0(group, "__vs__", vs),
                           NA_character_)
  
  tibble(
    dataset_id = dataset_id,
    tissue = tissue,
    time_h = time_h,
    time_bin = time_bin_from_h(time_h),
    group = group,
    vs = vs,
    contrast_label = contrast_label,
    file = path,
    file_name = bn
  )
}

# Read minimal info from a DE file to get n_control / n_delivery
# (they are repeated per row, so reading a few rows is enough)
read_de_sample_counts <- function(path) {
  df <- read_tsv(path, show_col_types = FALSE, n_max = 5)
  n_control  <- if ("n_control"  %in% names(df)) suppressWarnings(as.integer(df$n_control[1])) else NA_integer_
  n_delivery <- if ("n_delivery" %in% names(df)) suppressWarnings(as.integer(df$n_delivery[1])) else NA_integer_
  tibble(n_control = n_control, n_delivery = n_delivery)
}

# Total samples in dataset (from normalized counts columns)
get_total_samples_in_dataset <- function(dataset_id) {
  norm_dir <- file.path(de_root, dataset_id, "normalized")
  counts_candidates <- c(
    file.path(norm_dir, "gene_counts_normalized.tsv"),
    file.path(norm_dir, "gene_counts_normailzed.tsv")
  )
  counts_path <- counts_candidates[file.exists(counts_candidates)][1]
  if (is.na(counts_path)) return(NA_integer_)
  
  # Read header only
  hdr <- names(read_tsv(counts_path, show_col_types = FALSE, n_max = 0))
  if (length(hdr) == 0) return(NA_integer_)
  
  # assume first column is gene_id-ish
  as.integer(max(0, length(hdr) - 1))
}

# ----------------------------
# 2) Build manifest (one row per contrast file)
# ----------------------------
manifest <- map_dfr(de_files, parse_meta_from_filename) %>%
  mutate(sample_counts = map(file, read_de_sample_counts)) %>%
  tidyr::unnest(sample_counts) %>%
  mutate(n_samples_in_contrast = n_control + n_delivery)

# attach total samples available per dataset
total_samples_tbl <- manifest %>%
  distinct(dataset_id) %>%
  mutate(total_samples_in_dataset = map_int(dataset_id, get_total_samples_in_dataset))

manifest <- manifest %>%
  left_join(total_samples_tbl, by = "dataset_id") %>%
  mutate(
    selection_flag = case_when(
      time_bin %in% c("anchor_lt24h", "calib_24to72h") ~ "USED_FOR_MODEL",
      time_bin == "late_gt72h" ~ "NOT_FOR_MODEL_BUT_SCORE_LATER",
      TRUE ~ "UNKNOWN_TIME_CHECK_FILENAME"
    )
  )

# ----------------------------
# 3) Summaries you asked for
# ----------------------------

# (A) How many contrast files are selected to what, and sample counts
selection_summary <- manifest %>%
  group_by(dataset_id, time_bin, selection_flag) %>%
  summarise(
    n_contrast_files = n(),
    total_samples_in_dataset = first(total_samples_in_dataset),
    sum_n_control = sum(n_control, na.rm = TRUE),
    sum_n_delivery = sum(n_delivery, na.rm = TRUE),
    median_n_control = median(n_control, na.rm = TRUE),
    median_n_delivery = median(n_delivery, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(dataset_id, match(time_bin, c("anchor_lt24h","calib_24to72h","late_gt72h","time_unknown")))

# (B) Overall totals across the whole project
overall_summary <- manifest %>%
  group_by(time_bin, selection_flag) %>%
  summarise(
    n_contrast_files = n(),
    sum_n_control = sum(n_control, na.rm = TRUE),
    sum_n_delivery = sum(n_delivery, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(match(time_bin, c("anchor_lt24h","calib_24to72h","late_gt72h","time_unknown")))

print(overall_summary)

# ----------------------------
# 4) Write outputs
# ----------------------------
manifest_out <- file.path(comparison_root, "step6A_file_manifest.tsv")
summary_out  <- file.path(comparison_root, "step6A_selection_summary.tsv")
overall_out  <- file.path(comparison_root, "step6A_overall_summary.tsv")

write_tsv(manifest, manifest_out)
write_tsv(selection_summary, summary_out)
write_tsv(overall_summary, overall_out)

cat("\nWrote:\n", manifest_out, "\n", summary_out, "\n", overall_out, "\n")

# ----------------------------
# 5) Extract the subsets (in memory) if you want
# ----------------------------
anchor_files <- manifest %>% filter(time_bin == "anchor_lt24h") %>% pull(file)
calib_files  <- manifest %>% filter(time_bin == "calib_24to72h") %>% pull(file)
late_files   <- manifest %>% filter(time_bin == "late_gt72h") %>% pull(file)

cat("\nCounts of files by bin:\n")
cat("anchor (<24h): ", length(anchor_files), "\n")
cat("calib (24-72h): ", length(calib_files), "\n")
cat("late (>72h): ", length(late_files), "\n")
