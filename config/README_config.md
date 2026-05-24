# Configuration

Copy `config_template.yml` to `config.yml` before running repository scripts. Edit `project_root` to the local repository checkout and confirm the data, results, scripts, and logs directories.

The manuscript-output runner expects derived inputs, curated metadata, split definitions, and frozen gene weights to be present at the configured paths. HiPerGator/HPC steps are disabled by default with `allow_hpc_steps: false`.
