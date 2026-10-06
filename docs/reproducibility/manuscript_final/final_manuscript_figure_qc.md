# Final manuscript figure QC

The checks below refer to the actual regenerated final composites. Visual status is finalized after native-size and 6.5-inch rendered inspection.

## Figure 1 - PASS
- Manuscript alignment: A workflow; B compact global dataset/context landscape.
- Numerical state: 68 contrasts (40/14/11/3).
- Exclusion: no scored human GSE262515/HepG2 rows.
- Readability/clipping/redundancy/panel order: PASS; Figure 1B uses dataset | tissue | verified numeric time consistently, including `1-24 h` for the aggregated GSE264344 acute window.

## Figure 2 - PASS
- Manuscript alignment and panel order: A support counts; B missing anchor; C coefficient distribution; D top absolute coefficients.
- Numerical state: 287 genes; 178 all-five; 109 four-of-five; missing-anchor counts 18/13/33/10/35.
- Exclusion/redundancy: not applicable to scored human contrasts; no duplicated panel role.
- Readability/clipping: PASS.

## Figure 3 - PASS
- Manuscript alignment: A split-level primary-versus-extended distribution; B nested splits within three independent primary datasets.
- Numerical state: primary 14/14 positive; dataset means 10.739, 11.920, and 8.148.
- Exclusion/readability/clipping/panel order: PASS; A labels both analysis groups directly on the x axis and has no redundant legend.
- Redundancy: PASS; B is not a subsetted Figure 1B forest.

## Figure 4 - PASS
- Manuscript alignment: A biological-context boundary view; B GSE264344 time course.
- Numerical state: 13 context rows with category totals 6/2/2/1/1/1.
- Exclusion: no human HepG2 evidence.
- Readability/clipping/redundancy/panel order: PASS; all context points use data-driven scale expansion.

## Figure 5 - PASS
- Manuscript alignment and panel order: permutation intervals; group summary; leave-one-gene-out; maximum single-gene contribution.
- Numerical state: 68 permutation contrasts with 1,000 permutations each; 65 positive; 25 genes x 68 contrasts = 1,700 leave-one-out rows; 1,699 preserve direction.
- Permutation cross-check: the corrected final-68 run produces 49/68 outside the strict 95% null interval and 43/68 with BH-adjusted two-sided p < 0.05, matching manuscript text values 49/68 and 43/68.
- The old 48/42 figure state filtered two invalid human contrasts only after a 70-contrast permutation run; the corrected run excludes them before random draws and applies BH to exactly 68 tests.
- Panel D metric: mean_max_contribution_fraction (Overall mean = 3.3%; Maximum observed = 8.7%).
- Exclusion/readability/clipping/redundancy/panel order: PASS; complete compact group legend retained.

## Figure 6 - PASS
- Manuscript alignment: comparator benchmarking is a two-panel main figure.
- Numerical state/exclusion: comparator summaries were recomputed from the synchronized 68-contrast table; secondary support n=3.
- Readability/clipping/redundancy/panel order: PASS; display labels are ISG signature, Chemokine/inflammatory signature, Generic innate signature, and IMRS.
- Secondary-support ISG directionality is a measured 0/3 (0%) with zero missing values, not NA; the zero is explicitly labeled at the baseline.

## Supplementary Figure S1 - PASS
- Manuscript alignment: exactly two panels: A detailed faceted provenance; B context-category counts.
- Numerical state/exclusion: 68 valid scored contrasts and 13 context rows; no human HepG2 scoring.
- Readability/clipping/redundancy/panel order: PASS; A is more detailed and differently encoded than Figure 1B; B counts unique context-shifted contrasts.

## Supplementary Figure S2 - PASS
- Manuscript alignment: GO Biological Process, Reactome, and MSigDB Hallmark remain present.
- Numerical state/exclusion: enrichment results and mapped-gene inputs are scientifically unchanged.
- Readability/clipping/redundancy/panel order: PASS; all panels show Gene count before the mathematical-minus FDR legend.
