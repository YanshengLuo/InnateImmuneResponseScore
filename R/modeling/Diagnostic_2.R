library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

# -------------------------
# INPUT
# -------------------------
contrib_path <- "D:/IMRS_Project/05_score/human_transfer/qc/GSE190850_HUMAN__featurecounts_validation__top_contributors.tsv"
scores_path  <- "D:/IMRS_Project/05_score/human_transfer/scores/GSE190850_HUMAN__featurecounts_validation__imrs_scores.tsv"

contrib <- read_tsv(contrib_path, show_col_types = FALSE)
scores  <- read_tsv(scores_path, show_col_types = FALSE)

# -------------------------
# REBUILD ORIGINAL SCORE FROM CONTRIBUTIONS
# -------------------------
# top_contributors file contains ranked per-sample contributions.
# Since your file keeps top 20 only, this test evaluates stability within the
# reported contribution space, not all 137 genes.
orig_from_contrib <- contrib %>%
  group_by(sample_id) %>%
  summarise(
    score_from_top_contrib = sum(w_times_z, na.rm = TRUE),
    .groups = "drop"
  )

score_df <- scores %>%
  select(sample_id, condition_simple, imrs_raw, imrs_z) %>%
  inner_join(orig_from_contrib, by = "sample_id")

# -------------------------
# LEAVE-ONE-GENE-OUT
# -------------------------
all_genes <- sort(unique(contrib$gene_id))

loog_results <- lapply(all_genes, function(g) {
  tmp <- contrib %>%
    filter(gene_id != g) %>%
    group_by(sample_id) %>%
    summarise(
      loo_score = sum(w_times_z, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    right_join(score_df, by = "sample_id")

  # Compare against original contribution-based score
  pearson_r <- suppressWarnings(cor(tmp$score_from_top_contrib, tmp$loo_score, method = "pearson", use = "complete.obs"))
  spearman_r <- suppressWarnings(cor(tmp$score_from_top_contrib, tmp$loo_score, method = "spearman", use = "complete.obs"))

  diff_vec <- tmp$loo_score - tmp$score_from_top_contrib

  tibble(
    removed_gene = g,
    pearson_r = pearson_r,
    spearman_r = spearman_r,
    mean_abs_shift = mean(abs(diff_vec), na.rm = TRUE),
    max_abs_shift = max(abs(diff_vec), na.rm = TRUE)
  )
})

loog_tbl <- bind_rows(loog_results) %>%
  arrange(desc(max_abs_shift))

print(loog_tbl)

# -------------------------
# SAVE TABLE
# -------------------------
out_dir <- "D:/IMRS_Project/05_score/human_transfer/eval"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write_tsv(loog_tbl, file.path(out_dir, "step09_human_leave_one_gene_out.tsv"))

# -------------------------
# PLOTS
# -------------------------
p1 <- ggplot(loog_tbl, aes(x = pearson_r)) +
  geom_histogram(bins = 30) +
  theme_bw() +
  labs(
    title = "Leave-one-gene-out stability",
    x = "Pearson correlation with original score",
    y = "Number of removed genes"
  )

ggsave(file.path(out_dir, "step09_human_leave_one_gene_out_corr.png"), p1, width = 8, height = 6, dpi = 300)
ggsave(file.path(out_dir, "step09_human_leave_one_gene_out_corr.pdf"), p1, width = 8, height = 6)

p2 <- ggplot(loog_tbl, aes(x = reorder(removed_gene, max_abs_shift), y = max_abs_shift)) +
  geom_col() +
  coord_flip() +
  theme_bw() +
  labs(
    title = "Maximum score shift after removing one gene",
    x = "Removed gene",
    y = "Max absolute score shift"
  )

ggsave(file.path(out_dir, "step09_human_leave_one_gene_out_shift.png"), p2, width = 9, height = 7, dpi = 300)
ggsave(file.path(out_dir, "step09_human_leave_one_gene_out_shift.pdf"), p2, width = 9, height = 7)

# -------------------------
# QUICK SUMMARY
# -------------------------
summary_tbl <- loog_tbl %>%
  summarise(
    n_genes_tested = n(),
    min_pearson_r = min(pearson_r, na.rm = TRUE),
    median_pearson_r = median(pearson_r, na.rm = TRUE),
    mean_pearson_r = mean(pearson_r, na.rm = TRUE),
    max_score_shift = max(max_abs_shift, na.rm = TRUE),
    median_score_shift = median(max_abs_shift, na.rm = TRUE)
  )

print(summary_tbl)
write_tsv(summary_tbl, file.path(out_dir, "step09_human_leave_one_gene_out_summary.tsv"))