#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(tibble)
  library(forcats)
  library(patchwork)
  library(purrr)
})

# ============================================================
# Figure 2 — IMRS captures strong, weak, and heterogeneous
# innate immune activation across independent datasets
#
# Inputs:
#   D:/IMRS_Project/05_score/failure_diagnosis/all_dataset_imrs_summary.tsv
#   D:/IMRS_Project/05_score/failure_diagnosis/diagnosis/diagnosis_table.tsv
#   D:/IMRS_Project/05_score/noise_rescore/rescore_summary.tsv
#   D:/IMRS_Project/05_score/transfer/scores/*__imrs_scores.tsv
#
# Outputs:
#   D:/IMRS_Project/06_figures/Figure2_IMRS_response_regimes.png
#   D:/IMRS_Project/06_figures/Figure2_IMRS_response_regimes.pdf
#   D:/IMRS_Project/06_figures/Figure2_legend.txt
# ============================================================

# -------------------------
# CONFIG
# -------------------------
project_root <- "D:/IMRS_Project"

all_summary_file <- file.path(
  project_root, "05_score", "failure_diagnosis", "all_dataset_imrs_summary.tsv"
)

diagnosis_file <- file.path(
  project_root, "05_score", "failure_diagnosis", "diagnosis", "diagnosis_table.tsv"
)

rescore_summary_file <- file.path(
  project_root, "05_score", "noise_rescore", "rescore_summary.tsv"
)

scores_dir <- file.path(
  project_root, "05_score", "transfer", "scores"
)

figure_dir <- file.path(project_root, "06_figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

out_png <- file.path(figure_dir, "Figure2_IMRS_response_regimes.png")
out_pdf <- file.path(figure_dir, "Figure2_IMRS_response_regimes.pdf")
out_legend <- file.path(figure_dir, "Figure2_legend.txt")

anchor_ids <- c("GSE39129", "GSE167521", "GSE264344", "GSE279744")

# thresholds for display classification
strong_threshold <- 10
weak_threshold <- 0.5

# -------------------------
# HELPERS
# -------------------------
stop_if_missing <- function(path) {
  if (!file.exists(path)) stop("Missing required file: ", path)
}

safe_read_tsv <- function(path) {
  read_tsv(path, show_col_types = FALSE, progress = FALSE)
}

classify_display_category <- function(mean_val, sd_val, diagnosis = NA_character_) {
  # If diagnosis explicitly says heterogeneous, honor that
  if (!is.na(diagnosis) && diagnosis %in% c("noisy_dataset", "heterogeneous_response")) {
    return("heterogeneous")
  }

  if (is.na(mean_val)) return("unknown")

  # strong = clearly high signal
  if (mean_val > strong_threshold) return("strong")

  # weak = low signal
  if (mean_val <= strong_threshold) {
    # use SD to split weak vs heterogeneous
    if (!is.na(sd_val) && sd_val >= 3) {
      return("heterogeneous")
    } else {
      return("weak")
    }
  }

  "unknown"
}

read_score_file_from_id <- function(id, scores_dir) {
  path <- file.path(scores_dir, paste0(id, "__imrs_scores.tsv"))
  if (!file.exists(path)) return(NULL)
  df <- safe_read_tsv(path)

  req <- c("sample_id", "condition_simple", "imrs_z")
  if (!all(req %in% names(df))) return(NULL)

  df %>%
    mutate(
      condition_simple = toupper(as.character(condition_simple)),
      imrs_z = as.numeric(imrs_z)
    )
}

pick_rep_dataset <- function(df, class_name) {
  x <- df %>% filter(display_category == class_name)

  if (nrow(x) == 0) return(NA_character_)

  if (class_name == "strong") {
    return(as.character(x %>% arrange(desc(delivery_mean_imrs_z)) %>% slice(1) %>% pull(dataset_id)))
  }
  if (class_name == "weak") {
    return(as.character(x %>% arrange(delivery_mean_imrs_z, delivery_sd_imrs_z) %>% slice(1) %>% pull(dataset_id)))
  }
  if (class_name == "heterogeneous") {
    return(as.character(x %>% arrange(desc(delivery_sd_imrs_z), desc(delivery_mean_imrs_z)) %>% slice(1) %>% pull(dataset_id)))
  }

  as.character(x %>% slice(1) %>% pull(dataset_id))
}

# -------------------------
# LOAD DATA
# -------------------------
stop_if_missing(all_summary_file)
stop_if_missing(diagnosis_file)
stop_if_missing(rescore_summary_file)

all_df <- safe_read_tsv(all_summary_file)
diag_df <- safe_read_tsv(diagnosis_file)
rescore_df <- safe_read_tsv(rescore_summary_file)

required_all <- c("dataset_id", "id", "delivery_mean_imrs_z", "delivery_sd_imrs_z")
missing_all <- setdiff(required_all, names(all_df))
if (length(missing_all) > 0) {
  stop("all_dataset_imrs_summary.tsv missing columns: ", paste(missing_all, collapse = ", "))
}

if (!("diagnosis" %in% names(diag_df))) {
  diag_df <- diag_df %>% mutate(diagnosis = NA_character_)
}

# join diagnosis onto full summary
panel_a_df <- all_df %>%
  left_join(
    diag_df %>% select(dataset_id, diagnosis) %>% distinct(),
    by = "dataset_id"
  ) %>%
  mutate(
    dataset_type = ifelse(dataset_id %in% anchor_ids, "anchor", "external"),
    display_category = purrr::pmap_chr(
      list(delivery_mean_imrs_z, delivery_sd_imrs_z, diagnosis),
      classify_display_category
    )
  ) %>%
  filter(display_category %in% c("strong", "weak", "heterogeneous")) %>%
  mutate(
    dataset_id = as.character(dataset_id),
    dataset_id_f = fct_reorder(dataset_id, delivery_mean_imrs_z)
  )

# -------------------------
# PANEL A
# -------------------------
panel_a <- ggplot(
  panel_a_df,
  aes(x = dataset_id_f, y = delivery_mean_imrs_z, fill = display_category)
) +
  geom_col(width = 0.75) +
  geom_errorbar(
    aes(
      ymin = pmax(delivery_mean_imrs_z - delivery_sd_imrs_z, 0),
      ymax = delivery_mean_imrs_z + delivery_sd_imrs_z
    ),
    width = 0.2,
    linewidth = 0.45
  ) +
  coord_flip() +
  scale_fill_manual(
    values = c(
      "strong" = "#C43C39",
      "weak" = "#3C78D8",
      "heterogeneous" = "#E69138"
    )
  ) +
  labs(
    title = "A  Dataset-level IMRS response summary",
    x = "Dataset",
    y = "Mean IMRS in delivery samples",
    fill = "Response class"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

# -------------------------
# PANEL B
# use original score files, not only noise-rescore samples
# -------------------------
rep_strong <- pick_rep_dataset(panel_a_df, "strong")
rep_weak <- pick_rep_dataset(panel_a_df, "weak")
rep_hetero <- pick_rep_dataset(panel_a_df, "heterogeneous")

message("rep_strong: ", ifelse(is.na(rep_strong), "NA", rep_strong))
message("rep_weak: ", ifelse(is.na(rep_weak), "NA", rep_weak))
message("rep_hetero: ", ifelse(is.na(rep_hetero), "NA", rep_hetero))

rep_map <- panel_a_df %>%
  filter(dataset_id %in% c(rep_strong, rep_weak, rep_hetero)) %>%
  select(dataset_id, id, display_category) %>%
  distinct() %>%
  mutate(dataset_id = as.character(dataset_id))

panel_b_list <- lapply(seq_len(nrow(rep_map)), function(i) {
  ds_id <- rep_map$dataset_id[i]
  obj_id <- rep_map$id[i]
  catg <- rep_map$display_category[i]

  sc <- read_score_file_from_id(obj_id, scores_dir)
  if (is.null(sc)) return(NULL)

  sc %>%
    mutate(
      dataset_id = ds_id,
      display_category = catg,
      facet_label = case_when(
        ds_id == rep_strong ~ paste0(ds_id, "\nStrong"),
        ds_id == rep_weak ~ paste0(ds_id, "\nWeak"),
        ds_id == rep_hetero ~ paste0(ds_id, "\nHeterogeneous"),
        TRUE ~ ds_id
      )
    )
})

panel_b_df <- bind_rows(panel_b_list)

if (nrow(panel_b_df) == 0) {
  panel_b <- ggplot() +
    annotate("text", x = 1, y = 1, label = "No sample-level data available for Panel B", size = 5) +
    theme_void() +
    labs(title = "B  Representative sample-level IMRS distributions")
} else {
  panel_b <- ggplot(
    panel_b_df,
    aes(x = condition_simple, y = imrs_z, color = condition_simple)
  ) +
    geom_violin(trim = FALSE, alpha = 0.15, linewidth = 0.4) +
    geom_jitter(width = 0.12, height = 0, size = 2, alpha = 0.85) +
    facet_wrap(~facet_label, scales = "free_y", nrow = 1) +
    scale_color_manual(
      values = c("CONTROL" = "#4D4D4D", "DELIVERY" = "#B22222")
    ) +
    labs(
      title = "B  Representative sample-level IMRS distributions",
      x = "Condition",
      y = "IMRS z-score"
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold"),
      legend.position = "none"
    )
}

# -------------------------
# PANEL C
# -------------------------
required_rescore <- c("dataset_id", "before_delta_imrs", "after_delta_imrs")
missing_rescore <- setdiff(required_rescore, names(rescore_df))
if (length(missing_rescore) > 0) {
  stop("rescore_summary.tsv missing columns: ", paste(missing_rescore, collapse = ", "))
}

panel_c_df <- rescore_df %>%
  left_join(
    panel_a_df %>% select(dataset_id, display_category) %>% distinct(),
    by = "dataset_id"
  ) %>%
  mutate(
    display_category = ifelse(is.na(display_category), "heterogeneous", display_category)
  )

xy_max <- max(c(panel_c_df$before_delta_imrs, panel_c_df$after_delta_imrs), na.rm = TRUE)
xy_min <- min(c(panel_c_df$before_delta_imrs, panel_c_df$after_delta_imrs), na.rm = TRUE)

panel_c <- ggplot(
  panel_c_df,
  aes(x = before_delta_imrs, y = after_delta_imrs, color = display_category)
) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.5) +
  geom_point(size = 3) +
  geom_text(
    aes(label = dataset_id),
    nudge_y = 0.03 * (xy_max - xy_min),
    size = 3,
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = c(
      "strong" = "#C43C39",
      "weak" = "#3C78D8",
      "heterogeneous" = "#E69138"
    )
  ) +
  labs(
    title = "C  Robustness to outlier removal",
    x = expression(Delta*"IMRS before outlier removal"),
    y = expression(Delta*"IMRS after outlier removal"),
    color = "Response class"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "none"
  )

# -------------------------
# COMBINE
# -------------------------
figure2 <- (panel_a / panel_b / panel_c) +
  plot_layout(heights = c(1.25, 1, 1)) &
  theme(plot.margin = margin(5.5, 8, 5.5, 5.5))

ggsave(
  filename = out_png,
  plot = figure2,
  width = 10.5,
  height = 14,
  dpi = 400,
  bg = "white"
)

ggsave(
  filename = out_pdf,
  plot = figure2,
  width = 10.5,
  height = 14,
  bg = "white"
)

# -------------------------
# LEGEND
# -------------------------
figure_legend <- paste(
  "Figure 2. IMRS captures strong, weak, and heterogeneous innate immune activation across independent datasets.",
  "",
  "(A) Dataset-level summary of IMRS across evaluated datasets. For each dataset, the mean IMRS of delivery-treated samples is shown with error bars representing the standard deviation across delivery samples. Datasets are grouped into strong, weak, and heterogeneous response classes based on score magnitude and dispersion, with heterogeneous datasets additionally informed by diagnostic classification. This panel shows that IMRS distinguishes distinct response regimes across studies while preserving the project’s per-dataset interpretation framework rather than relying on joint cross-dataset normalization.",
  "",
  "(B) Representative sample-level IMRS distributions for one strong, one weak, and one heterogeneous dataset. Each point represents a sample, with violin overlays showing the within-dataset score distribution for CONTROL and DELIVERY groups. Strong-response datasets show clear separation between control and delivery samples. Weak-response datasets show limited upward shift. Heterogeneous-response datasets show broad dispersion among delivery samples without discrete extreme points, indicating distributed variability rather than isolated outliers.",
  "",
  "(C) Noise robustness analysis based on outlier-removal rescoring. For datasets selected for noise-sensitivity testing, Delta IMRS (mean DELIVERY minus mean CONTROL) is plotted before and after removal of samples flagged using a robust median absolute deviation criterion. Points clustering near the identity line indicate that outlier removal does not materially alter the score, supporting the conclusion that elevated variance in these datasets reflects intrinsic biological heterogeneity rather than technical artifacts.",
  "",
  "All scores were generated using frozen anchor-derived gene weights and the same dataset-internal scoring workflow. Controls were used to define the within-dataset baseline, consistent with the IMRS framework in which absolute values are interpreted within, rather than across, datasets.",
  sep = "\n"
)

writeLines(figure_legend, out_legend)

message("Representative datasets used in Panel B:")
message("  strong: ", ifelse(is.na(rep_strong), "NA", rep_strong))
message("  weak: ", ifelse(is.na(rep_weak), "NA", rep_weak))
message("  heterogeneous: ", ifelse(is.na(rep_hetero), "NA", rep_hetero))

message("Saved:")
message("  ", out_png)
message("  ", out_pdf)
message("  ", out_legend)