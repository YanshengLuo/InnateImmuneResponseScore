# Configuration

Copy `config_template.yml` to `config.yml` before running repository scripts. The template uses repository-relative paths; `project_root: "."` is appropriate when running from the repository root.

The Layer 2 manuscript-output runner expects released inputs under `data/derived/`, including the explicit `data/derived/figure_inputs/` bundle used by the repo-contained figure implementation. It writes generated artifacts beneath `results_release_templates/`.

Keep `execute_active_scripts: false` to run the preflight/checklist only. Set it to `true` only for controlled regeneration of manuscript figures, Priority3 enrichment outputs, and supplementary tables.

`run_internal_readiness: false` excludes optional NAR G&B readiness packaging from the default reviewer-facing run. Full computational provenance and HiPerGator/HPC operations remain separate from Layer 2.
