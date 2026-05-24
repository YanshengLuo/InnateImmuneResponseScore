library(readr)
library(dplyr)
library(ggplot2)

# path to your scoring file
scores <- read_tsv(
"D:/IMRS_Project/05_score/transfer/scores/GSE262515__featurecounts_validation__imrs_scores.tsv",
show_col_types = FALSE
)

# ensure ordering
scores$condition_simple <- factor(
scores$condition_simple,
levels = c("CONTROL","DELIVERY")
)

# plot
p <- ggplot(scores, aes(x = condition_simple, y = imrs_z, fill = condition_simple)) +
  geom_violin(trim = FALSE, alpha = 0.4, color = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA) +
  geom_jitter(width = 0.05, size = 2, alpha = 0.8) +
  theme_bw(base_size = 14) +
  labs(
    title = "IMRS Score Distribution — GSE262515",
    x = "",
    y = "IMRS z-score"
  ) +
  scale_fill_manual(values = c("CONTROL" = "#4C72B0", "DELIVERY" = "#DD8452")) +
  theme(legend.position = "none")

print(p)

wilcox.test(imrs_z ~ condition_simple, data = scores)