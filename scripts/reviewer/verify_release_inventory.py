#!/usr/bin/env python3
"""Verify the release's canonical Git blobs independently of checkout line endings."""

import argparse
import csv
import hashlib
import io
from pathlib import Path
import subprocess
import zipfile

ROOT = Path(__file__).resolve().parents[2]
INVENTORY = "docs/reproducibility/reproduce_package_inventory.tsv"


def git(*args):
    return subprocess.check_output(["git", *args], cwd=ROOT)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--staged", action="store_true", help="Verify the Git index instead of HEAD")
    parser.add_argument("--write-inventory", action="store_true", help="Write checksums from the staged release")
    args = parser.parse_args()
    staged = args.staged or args.write_inventory
    if staged:
        names = git("ls-files", "-z").decode().split("\0")
    else:
        names = git("ls-tree", "-r", "--name-only", "-z", "HEAD").decode().split("\0")
    names = sorted(name for name in names if name and name != INVENTORY)
    prefix = ":" if staged else "HEAD:"
    records = []
    for name in names:
        content = git("show", prefix + name)
        records.append({"path": name, "bytes": str(len(content)), "sha256": hashlib.sha256(content).hexdigest()})
    if args.write_inventory:
        with (ROOT / INVENTORY).open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=["path", "bytes", "sha256"], delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(records)
    else:
        manifest = git("show", prefix + INVENTORY).decode("utf-8-sig")
        expected = list(csv.DictReader(io.StringIO(manifest), delimiter="\t"))
        if records != expected:
            raise SystemExit("Release inventory differs from canonical Git content")

    stems = [f"Figure{i}_main_v5" for i in range(1, 7)] + ["FigureS1_main_v5", "FigureS2_gene_program_enrichment_combined"]
    for stem in stems:
        for extension in ("png", "pdf", "svg"):
            name = f"publication_outputs/figures/{stem}.{extension}"
            if name not in names or not git("show", prefix + name):
                raise SystemExit(f"Missing release figure: {name}")

    name = "data/derived/figure_inputs/label_permutation_null_summary.tsv"
    rows = list(csv.DictReader(io.StringIO(git("show", prefix + name).decode()), delimiter="\t"))
    outside = sum(row["observed_outside_95pct_null"].upper() == "TRUE" for row in rows)
    significant = sum(float(row["empirical_p_two_sided_fdr"]) < 0.05 for row in rows)
    if (len(rows), outside, significant) != (68, 49, 43):
        raise SystemExit(f"Permutation release invariant failed: {len(rows)}/{outside}/{significant}")
    if any(int(row["n_permutations"]) != 1000 for row in rows):
        raise SystemExit("Expected 1,000 permutations per contrast")
    mirror = "data/derived/label_permutation_null_summary.tsv"
    if git("show", prefix + name) != git("show", prefix + mirror):
        raise SystemExit("Permutation source copies differ")

    ppt_name = "publication_outputs/IMRS_Final_Figures_R_Assembled_PublicationReady.pptx"
    with zipfile.ZipFile(io.BytesIO(git("show", prefix + ppt_name))) as ppt:
        if ppt.testzip() is not None:
            raise SystemExit("PowerPoint ZIP integrity failed")
        slides = [name for name in ppt.namelist() if name.startswith("ppt/slides/slide") and name.endswith(".xml") and "/_rels/" not in name]
        if len(slides) != 8:
            raise SystemExit(f"Expected 8 PowerPoint slides; found {len(slides)}")

    forbidden = ("tmp/", ".codex-build/", "node_modules/", "results_release_templates/", "renv/library/")
    if any(name.startswith(forbidden) for name in names):
        raise SystemExit("Local build or generated working files are tracked")
    oversized = [record["path"] for record in records if int(record["bytes"]) >= 100 * 1024 * 1024]
    if oversized:
        raise SystemExit("Files exceed GitHub's normal 100 MiB limit: " + ", ".join(oversized))
    print(f"PASS: {len(records)} release files; 24 figure files; 8 PPT slides; permutation counts 68/49/43")


if __name__ == "__main__":
    main()
