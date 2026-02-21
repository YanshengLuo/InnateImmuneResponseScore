suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
})

in_path  <- "C:/Users/john/Desktop/IMRS_Project/05_score/imrs_scores_samples_mouse.tsv"
out_dir  <- "C:/Users/john/Desktop/IMRS_Project/05_score"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

df <- read_tsv(in_path, show_col_types = FALSE) %>%
  filter(is_human == FALSE) %>%
  mutate(
    is_control = delivery_group %in% c("baseline_0", "control", "ctrl", "vehicle", "mock")
  )

# 1) Quick QC: core gene coverage per dataset+tissue stratum
qc_core <- df %>%
  group_by(dataset_id, tissue, stratum) %>%
  summarise(
    n_samples = n(),
    core_presence_fraction = first(core_presence_fraction),
    n_core_genes_present = first(n_core_genes_present),
    .groups = "drop"
  ) %>%
  arrange(dataset_id, tissue)

write_tsv(qc_core, file.path(out_dir, "mouse_qc_core_by_stratum.tsv"))

# 2) Baseline vs delivery summary within each stratum
# (only meaningful if stratum has at least 2 controls and 2 non-controls)
effect_tbl <- df %>%
  group_by(dataset_id, tissue, stratum) %>%
  mutate(
    n_ctrl = sum(is_control),
    n_trt  = sum(!is_control)
  ) %>%
  ungroup() %>%
  filter(n_ctrl >= 2, n_trt >= 2) %>%
  group_by(dataset_id, tissue, stratum) %>%
  summarise(
    n_ctrl = sum(is_control),
    n_trt  = sum(!is_control),
    ctrl_mean = mean(imrs_norm[is_control], na.rm = TRUE),
    trt_mean  = mean(imrs_norm[!is_control], na.rm = TRUE),
    delta_mean = trt_mean - ctrl_mean,
    ctrl_sd = sd(imrs_norm[is_control], na.rm = TRUE),
    trt_sd  = sd(imrs_norm[!is_control], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(delta_mean))

write_tsv(effect_tbl, file.path(out_dir, "mouse_effect_baseline_vs_delivery_by_stratum.tsv"))

# 3) Delivery-group level summary (most actionable)
group_tbl <- df %>%
  group_by(dataset_id, tissue, stratum, delivery_group) %>%
  summarise(
    n = n(),
    mean_imrs = mean(imrs_norm, na.rm = TRUE),
    sd_imrs   = sd(imrs_norm, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(dataset_id, tissue, desc(mean_imrs))

write_tsv(group_tbl, file.path(out_dir, "mouse_group_summary_by_stratum.tsv"))

message("DONE (mouse-only QC + summaries).")
message("Wrote:")
message(" - ", file.path(out_dir, "mouse_qc_core_by_stratum.tsv"))
message(" - ", file.path(out_dir, "mouse_effect_baseline_vs_delivery_by_stratum.tsv"))
message(" - ", file.path(out_dir, "mouse_group_summary_by_stratum.tsv"))