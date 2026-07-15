library(readr)
library(dplyr)
library(ggplot2)

contrib <- read_tsv(
  "D:/IMRS_Project/05_score/human_transfer/qc/GSE190850_HUMAN__featurecounts_validation__top_contributors.tsv"
)

# compute contribution share
contrib2 <- contrib %>%
  group_by(sample_id) %>%
  mutate(
    total_abs = sum(abs_w_times_z),
    frac = abs_w_times_z / total_abs
  )

# cumulative contribution
cum <- contrib2 %>%
  arrange(sample_id, rank) %>%
  group_by(sample_id) %>%
  mutate(cum_frac = cumsum(frac))

# average across samples
avg <- cum %>%
  group_by(rank) %>%
  summarise(mean_cum_frac = mean(cum_frac), .groups="drop")

ggplot(avg, aes(rank, mean_cum_frac)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  labs(
    title="IMRS contribution concentration",
    x="Top contributing genes (rank)",
    y="Cumulative fraction of score"
  )

  