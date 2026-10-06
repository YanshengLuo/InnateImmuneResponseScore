# Final figure polish QC

1. Permutation result: 49/68 contrasts outside the strict 95% null interval.
2. BH result: 43/68 two-sided empirical permutation p-values remain significant after BH adjustment over 68 tests.
3. Discrepancy resolved: the stale 70-contrast run consumed random draws for two later-excluded HepG2 contrasts; post hoc filtering yielded 48/42, while correct pre-permutation exclusion yields 49/43.
4. Manuscript correction: not needed because the validated result matches 49/43.
5. Figure 3A: Primary acute validation and Extended validation are labeled directly on the x axis; the redundant fill legend is removed; annotations and values are unchanged.
6. Internal terminology: `Interpretation support level` and audit-row wording are absent from publication figures.
7. Zero baselines: all six publication bar/histogram value axes use lower limit 0 and lower expansion 0; see `bargraph_zero_baseline_qc.tsv`.
8. Comparator zero versus missing: Secondary-support ISG directionality is measured 0/3 (0%) with zero missing values, not NA, and is labeled 0% at the baseline.
9. Figure 5D: the dashed reference line retains the 3.3% overall mean, and the point-linked annotation reads `Maximum observed = 8.7%`; the underlying fraction is unchanged.
10. Figure 1B and Supplementary Figure S1A: dataset, tissue, and verified numeric time are shown in a consistent order; Figure 1B no longer mixes `acute` or omitted times with numeric time labels; the GSE264344 aggregate is labeled `1-24 h`.
11. Supplementary Figure S2: Gene count precedes the mathematical-minus log10(FDR) legend in A/B/C; term content, ranking, sizes, and blue gradient are unchanged.
12. Time course: redundant 72 h point labels were removed because 72 h is already an x-axis tick; values and trajectories are unchanged.
13. Native-size and approximately 6.5-inch visual QC: Figures 1-6 and Supplementary Figures S1-S2 PASS for clipping, collisions, legends, labels, points, zero/NA clarity, percentages, and publication terminology.
