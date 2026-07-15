#!/usr/bin/env Rscript

# Full IMRS pipeline rebuild template from public data.
# This is not the reviewer manuscript-output runner. It outlines the full workflow and may require HiPerGator/HPC.

# 1. Retrieve public data from accessions.
# TODO: use data/accessions/public_dataset_accessions.tsv and locally configured retrieval/alignment scripts.
# Account-specific HPC submission scripts are intentionally excluded from the public release candidate.

# 2. Prepare metadata.
# TODO: run metadata processing scripts after config/path cleanup.

# 3. Generate delivery-versus-control split designs.
# TODO: use curated metadata and split design scripts to create data/split_designs/.

# 4. Run anchor differential expression.
# TODO: run anchor delivery-versus-control DE scripts, likely on HPC for large count matrices.

# 5. Select retained genes.
# TODO: run locked-anchor retained-gene selection and support/heterogeneity/power scripts.

# 6. Estimate frozen coefficients.
# TODO: estimate frozen IMRS gene coefficients and write frozen_gene_weights.

# 7. Score samples with frozen IMRS.
# TODO: apply frozen scoring to dataset-internal normalized/control-referenced matrices.

# 8. Evaluate split-level delivery-minus-control ΔIMRSz.
# TODO: evaluate split-level delivery-minus-control ΔIMRSz, directionality, and secondary AUC.

# 9. Run robustness analyses.
# TODO: run label permutation, leave-one-gene-out, gene dominance, threshold sensitivity, and leave-one-anchor-out analyses where source scripts are available.

# 10. Generate figures/tables.
# TODO: call run_all_manuscript_outputs_v6.R after derived inputs are created and active scripts are path-portable.

message('This is a full-pipeline template. It intentionally contains TODO markers for raw-data and HPC-dependent steps.')
