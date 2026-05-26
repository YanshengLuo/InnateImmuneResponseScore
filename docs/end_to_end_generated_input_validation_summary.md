# End-To-End Generated-Input Validation Summary

## Validation Record

The generated-input route was validated in a local validation environment on
May 26, 2026, completing at 01:32:38 EDT.

## Command Sequence

```sh
Rscript scripts/portable_full_pipeline/run_count_to_v6_outputs.R --config config/local_d_drive_config_TEMPLATE.yml --mode all_scored --force
Rscript scripts/portable_full_pipeline/run_step09_to_layer2_inputs.R --config config/local_d_drive_config_TEMPLATE.yml --mode all_scored --force
Rscript run_all_manuscript_outputs_v6.R --config <output_root>/layer2_generated_inputs/layer2_generated_inputs_config.yml
```

For the final Layer 2 validation execution, active-script execution was
enabled in the generated local configuration. The repository-facing
generated config remains suitable for a dry-run/checklist invocation unless
the user explicitly enables output generation.

## Stages Completed

- Layer 1 preflight and design/split preparation completed.
- Anchor-only DESeq2 and Steps 06A-07 reconstructed frozen weights from the five locked anchors.
- The regenerated frozen weight table matched the released canonical table within tolerance.
- Step08 and Step09 applied the frozen weights to available manuscript datasets in `all_scored` mode.
- The ported distributed Step09-to-Layer2 bridge generated manuscript-facing input tables and a manifest.
- The generated-input Layer 2 run completed manuscript figures, Priority3 gene-program enrichment, and supplementary tables with exit status 0.
- `internal_readiness` was disabled and not run.

## Step09 Validation Datasets

In addition to the locked anchors, generated Step09 evaluation included:

```text
GSE119119
GSE139529
GSE166655
GSE178313
GSE279743
GSE314070
GSE262515_cell_line
GSE262515_tissue
```

Validation datasets were scored with frozen locked-anchor coefficients; they
did not participate in Steps 06A-07 and did not refit the model.

## Generated Layer 2 Outputs

The active generated-input run wrote its manuscript-output products beneath:

```text
<output_root>/layer2_generated_inputs/manuscript_outputs/
```

Generated bridge inputs are staged separately from committed release inputs.
Released `data/derived/` files were not overwritten. The bridge manifest may
label individual Layer 2-required inputs as released references when no
direct original generator was available in the ported bridge components.

## Excluded Operations

This validation did not run raw-data retrieval, alignment or feature
counting, HPC/SLURM jobs, archive scripts, or human exploratory workflows.
