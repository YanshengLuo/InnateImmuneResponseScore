# Historical Layer 1 Scripts

This directory preserves the original Layer 1 count-to-score workflow scripts for computational provenance.

These scripts document the path from gene counts and verified metadata through contrast construction, locked-anchor retained-gene construction, frozen coefficient estimation, frozen IMRS scoring, and Step09 split-level evaluation. They may retain historical local path assumptions and workflow conventions and are not intended to be executed as-is from a fresh clone.

The validated reviewer-facing executable route is the repository-root `run_all_manuscript_outputs_v6.R` Layer 2 runner, which regenerates manuscript outputs from released derived inputs and does not invoke this directory.
