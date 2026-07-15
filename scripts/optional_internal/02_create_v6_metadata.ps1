$ErrorActionPreference = "Stop"

function Read-SimpleYamlConfig {
  param([string]$Path)
  $cfg = @{}
  if (-not (Test-Path -LiteralPath $Path)) { return $cfg }
  foreach ($line in Get-Content -LiteralPath $Path) {
    if ($line -notmatch '^\s*#' -and $line -match '^\s*([^:]+):\s*(.*)\s*$') {
      $key = $matches[1].Trim()
      $value = $matches[2].Trim().Trim('"').Trim("'")
      $cfg[$key] = $value
    }
  }
  return $cfg
}

$RepositoryRoot = $env:IMRS_REPOSITORY_ROOT
if (-not $RepositoryRoot) {
  $RepositoryRoot = Resolve-Path (Join-Path $PSScriptRoot "..\..")
}
$ConfigPath = Join-Path $RepositoryRoot "config\config.yml"
if (-not (Test-Path -LiteralPath $ConfigPath)) {
  $ConfigPath = Join-Path $RepositoryRoot "config\config_template.yml"
}
$Config = Read-SimpleYamlConfig -Path $ConfigPath
function Resolve-ConfigPath {
  param([string]$Value, [string]$Default)
  if (-not $Value) { $Value = $Default }
  if ([System.IO.Path]::IsPathRooted($Value)) { return $Value }
  return Join-Path $RepositoryRoot $Value
}

$SourceRoot = Resolve-ConfigPath $Config["source_v5_dir"] "data_release_templates\derived\revised_plots_v5"
$Root = Resolve-ConfigPath $Config["manuscript_output_dir"] "results_release_templates"
$PreHashPath = Resolve-ConfigPath $Config["source_v5_hashes"] "data_release_templates\derived\revised_plots_v5_pre_v6_hashes.tsv"

$RegeneratedPanels = @("Figure1A", "Figure1B", "Figure3B", "Figure5A")
$RegeneratedFigures = @("Figure1", "Figure3", "Figure5")
$MetadataFiles = @(
  "v6_change_log.tsv",
  "v6_terminology_audit.tsv",
  "v6_output_manifest.tsv",
  "v6_script_inventory.tsv",
  "v6_table_inventory.tsv",
  "v6_file_inventory.tsv",
  "v6_regeneration_log.txt",
  "v6_manual_caption_replacement_table.tsv"
)

function RelPath([string]$Path, [string]$Base = $Root) {
  $full = [System.IO.Path]::GetFullPath($Path)
  $baseFull = [System.IO.Path]::GetFullPath($Base).TrimEnd('\') + '\'
  if ($full.StartsWith($baseFull, [System.StringComparison]::OrdinalIgnoreCase)) {
    return $full.Substring($baseFull.Length)
  }
  return $full
}

function CleanCell($Value) {
  if ($null -eq $Value) { return "" }
  return ([string]$Value) -replace "`r?`n", " " -replace "`t", " " -replace "\s{2,}", " "
}

function Write-TsvRows([string]$Path, [array]$Rows, [string[]]$Columns) {
  $lines = New-Object System.Collections.Generic.List[string]
  $lines.Add(($Columns -join "`t"))
  foreach ($row in $Rows) {
    $cells = foreach ($column in $Columns) {
      CleanCell $row[$column]
    }
    $lines.Add(($cells -join "`t"))
  }
  [System.IO.File]::WriteAllLines($Path, $lines, [System.Text.UTF8Encoding]::new($false))
}

function FileHashOrEmpty([string]$Path) {
  if (!(Test-Path -LiteralPath $Path)) { return "" }
  return (Get-FileHash -LiteralPath $Path -Algorithm SHA256).Hash
}

function FileSizeOrZero([string]$Path) {
  if (!(Test-Path -LiteralPath $Path)) { return 0 }
  return (Get-Item -LiteralPath $Path).Length
}

function ExistingNonzeroSummary([string[]]$Paths) {
  $labels = @("png", "pdf", "svg")
  $bits = for ($i = 0; $i -lt $Paths.Count; $i++) {
    $exists = (Test-Path -LiteralPath $Paths[$i]) -and ((Get-Item -LiteralPath $Paths[$i]).Length -gt 0)
    "$($labels[$i])=$exists"
  }
  return ($bits -join ";")
}

function SizeSummary([string[]]$Paths) {
  $labels = @("png", "pdf", "svg")
  $bits = for ($i = 0; $i -lt $Paths.Count; $i++) {
    "$($labels[$i])=$(FileSizeOrZero $Paths[$i])"
  }
  return ($bits -join ";")
}

function PurposeForScript([string]$Rel) {
  switch -Wildcard ($Rel) {
    "scripts\00_generate_all_reorganized_figures_revised_v5.R" { "Full inherited figure-generation wrapper with v6 output root" }
    "scripts\01_regenerate_changed_v6.R" { "Targeted v6 regeneration for figures with plot-visible text changes" }
    "scripts\02_create_v6_metadata.ps1" { "Creates v6 change log, audits, inventories, manifest, and regeneration log" }
    "scripts\R\build_merged_imrs_workflow_v5.R" { "Renders Figure 1A workflow panel" }
    "scripts\R\v5_helpers.R" { "Shared plotting, assembly, text-replacement, and manifest helpers used by v6" }
    default { "Figure-supporting script" }
  }
}

function KeyChangesForScript([string]$Rel) {
  switch -Wildcard ($Rel) {
    "scripts\00_generate_all_reorganized_figures_revised_v5.R" { "Updated output root from revised_plots_v5 to revised_plots_v6 and v6 completion text" }
    "scripts\01_regenerate_changed_v6.R" { "New v6-only targeted regeneration script for Figure1, Figure3, and Figure5" }
    "scripts\02_create_v6_metadata.ps1" { "New metadata generator for required v6 TSV/log outputs" }
    "scripts\R\build_merged_imrs_workflow_v5.R" { "Updated Figure 1A title, subtitle, workflow node labels, bottom note, wrapping, box sizes, canvas, and margins" }
    "scripts\R\v5_helpers.R" { "Updated controlled terminology, Figure1A panel dimensions, Figure1B axis/legend, Figure3B size legend, Figure5A axis label, and v6 manifest notes" }
    default { "" }
  }
}

function PurposeForTable([string]$Rel) {
  switch -Wildcard ($Rel) {
    "v6_change_log.tsv" { "Required v6 change log" }
    "v6_terminology_audit.tsv" { "Required v6 controlled-terminology audit" }
    "v6_output_manifest.tsv" { "Required v6 output manifest" }
    "v6_script_inventory.tsv" { "Required v6 script inventory" }
    "v6_table_inventory.tsv" { "Required v6 table inventory" }
    "v6_file_inventory.tsv" { "Required v6 full file inventory" }
    "v6_manual_caption_replacement_table.tsv" { "Required manual caption replacement table because no editable caption source was present" }
    "figure_v5_manifest.tsv" { "Inherited long figure manifest updated to v6 paths" }
    "figure_v6_manifest.tsv" { "Copy of long figure manifest with v6 filename" }
    "tables\v5_panel_manifest.tsv" { "Inherited panel manifest updated to v6 paths" }
    "tables\v5_figure_manifest_wide.tsv" { "Inherited figure manifest updated to v6 paths" }
    "tables\v5_wording_audit.tsv" { "Inherited wording audit regenerated with v6 terminology mappings" }
    "tables\v6_wording_audit.tsv" { "v6 copy of wording audit regenerated with v6 terminology mappings" }
    default { "Figure-supporting table" }
  }
}

function KeyChangesForTable([string]$Rel, [bool]$Edited) {
  switch -Wildcard ($Rel) {
    "v6_*" { "Created for v6 required metadata outputs" }
    "figure_v5_manifest.tsv" { "Updated manifest paths to revised_plots_v6 after targeted regeneration" }
    "figure_v6_manifest.tsv" { "Created as v6-named copy of long figure manifest" }
    "tables\v5_panel_manifest.tsv" { "Updated panel paths and v6 terminology notes" }
    "tables\v5_figure_manifest_wide.tsv" { "Updated figure paths and v6 terminology notes" }
    "tables\v5_wording_audit.tsv" { "Regenerated with v6 controlled terminology mappings" }
    "tables\v6_wording_audit.tsv" { "Created as v6 wording audit copy" }
    default { if ($Edited) { "Generated or changed during v6 targeted regeneration" } else { "Copied unchanged from revised_plots_v5" } }
  }
}

function TsvShape([string]$Path) {
  $lines = [System.IO.File]::ReadAllLines($Path)
  if ($lines.Count -eq 0) { return @{ rows = 0; columns = 0 } }
  return @{ rows = [Math]::Max(0, $lines.Count - 1); columns = ($lines[0] -split "`t", -1).Count }
}

function TextContains([string]$Path, [string]$Needle) {
  if (!(Test-Path -LiteralPath $Path)) { return $false }
  $text = [System.IO.File]::ReadAllText($Path)
  return $text.Contains($Needle)
}

function CompareV5Hashes() {
  if (!(Test-Path -LiteralPath $PreHashPath)) {
    return @{ status = "not_checked"; detail = "Pre-run hash file not found: $PreHashPath" }
  }
  $before = @{}
  foreach ($line in [System.IO.File]::ReadAllLines($PreHashPath)) {
    if ([string]::IsNullOrWhiteSpace($line)) { continue }
    $parts = $line -split "`t", 2
    if ($parts.Count -eq 2) { $before[$parts[1]] = $parts[0] }
  }
  $currentFiles = Get-ChildItem -LiteralPath $SourceRoot -Recurse -File
  $changed = New-Object System.Collections.Generic.List[string]
  foreach ($file in $currentFiles) {
    $hash = (Get-FileHash -LiteralPath $file.FullName -Algorithm SHA256).Hash
    if (!$before.ContainsKey($file.FullName) -or $before[$file.FullName] -ne $hash) {
      $changed.Add($file.FullName)
    }
  }
  if ($changed.Count -eq 0 -and $before.Count -eq $currentFiles.Count) {
    return @{ status = "PASS"; detail = "$($currentFiles.Count) v5 files matched pre-run SHA256 hashes." }
  }
  return @{ status = "FAIL"; detail = (($changed | ForEach-Object { RelPath $_ $SourceRoot }) -join "; ") }
}

$manualCaptionRows = @(
  [ordered]@{
    figure_id = "Figure 1"
    current_caption_or_subtitle = "No editable Figure 1 caption source was present in revised_plots_v5 or revised_plots_v6."
    revised_caption_or_subtitle = "Figure 1. Frozen IMRS scoring framework and dataset-level delivery-minus-control ΔIMRSz landscape. (A) Overview of the frozen, transfer-oriented IMRS scoring framework. Delivery-versus-control split contrasts were defined from verified metadata and raw count matrices. Frozen IMRS gene coefficients were derived only from five locked-anchor datasets and were not refit during validation or transfer evaluation. Target datasets were scored using dataset-internal normalization, control-referenced gene z-scoring, weighted IMRS scoring, and control-standardized IMRSz values. (B) Dataset/tissue-level mean delivery-minus-control ΔIMRSz values across evaluated contexts. Point size indicates the number of passing split contrasts, and color indicates manuscript analysis group. This figure defines the scoring and evaluation framework; primary biological validation is evaluated in subsequent figures."
    where_to_apply_manually = "Manuscript caption source or assembled PowerPoint/Word caption block outside revised_plots_v6"
    notes = "caption_manual_update_needed; no editable caption source was found in v5/v6"
  },
  [ordered]@{
    figure_id = "Figure 2"
    current_caption_or_subtitle = "Figure 2. Frozen IMRS coefficients define a reproducible acute delivery-associated gene program. Body contains: Top-ranked genes by absolute coefficient include several genes compatible with innate and inflammatory biology. Body contains: Overall, these results support interpreting IMRS as a reproducible multi-gene acute delivery-associated transcriptional program rather than a signature driven by a single dataset or a single dominant gene."
    revised_caption_or_subtitle = "Figure 2. Frozen IMRS coefficients define a reproducible acute delivery-associated innate transcriptional program. Replace body phrase with: Top-ranked genes by absolute coefficient include genes compatible with chemokine/inflammatory and interferon-associated transcriptional biology. Replace body phrase with: Overall, these results support interpreting IMRS as a reproducible multi-gene acute delivery-associated innate transcriptional program rather than a signature driven by a single dataset or a single dominant gene."
    where_to_apply_manually = "Manuscript caption source or assembled PowerPoint/Word caption block outside revised_plots_v6"
    notes = "caption_manual_update_needed; no editable caption source was found in v5/v6"
  },
  [ordered]@{
    figure_id = "Figure 3"
    current_caption_or_subtitle = "Figure 3. IMRS transfers most strongly to independent primary acute validation datasets. Body contains: Primary acute validation showed consistent positive IMRS elevation. Body contains: Dataset-level summaries show that the strongest non-anchor validation responses were concentrated in primary acute validation datasets. Body contains: Together, these results support transfer of the frozen IMRS score to independent acute delivery-associated contexts."
    revised_caption_or_subtitle = "Figure 3. Frozen IMRS scoring transfers most strongly to independent primary acute validation datasets. Replace body phrase with: Primary acute validation showed consistent positive delivery-minus-control ΔIMRSz. Replace body phrase with: Dataset-level summaries show that the strongest non-anchor validation ΔIMRSz values were concentrated in primary acute validation datasets. Replace body phrase with: Together, these results support transfer of the frozen IMRS scoring framework to independent datasets showing an acute delivery-associated innate transcriptional response. Keep: while limiting the claim in late or context-shifted settings."
    where_to_apply_manually = "Manuscript caption source or assembled PowerPoint/Word caption block outside revised_plots_v6"
    notes = "caption_manual_update_needed; no editable caption source was found in v5/v6"
  },
  [ordered]@{
    figure_id = "Figure 4"
    current_caption_or_subtitle = "Figure 4. Late and context-shifted datasets define the biological scope of IMRS. Body contains: These lower responses were concentrated in late, therapeutic, tissue-shifted, or formulation-shifted settings. Body contains: rather than appearing as an unstructured failure pattern. Body contains: Together, these results support an acute delivery-associated interpretation of IMRS."
    revised_caption_or_subtitle = "Figure 4. Late and context-shifted settings bound the acute delivery-associated innate transcriptional response captured by IMRS. Replace body phrase with: These attenuated responses were concentrated in late, therapeutic, tissue-shifted, or formulation-shifted settings. Replace body phrase with: rather than indicating a general failure of the scoring framework. Replace body phrase with: Together, these results support interpreting IMRS as most responsive to an acute delivery-associated innate transcriptional response. Keep: late or context-shifted settings should be interpreted as boundary-setting evidence rather than primary validation."
    where_to_apply_manually = "Manuscript caption source or assembled PowerPoint/Word caption block outside revised_plots_v6"
    notes = "caption_manual_update_needed; no editable caption source was found in v5/v6"
  },
  [ordered]@{
    figure_id = "Figure 5"
    current_caption_or_subtitle = "Figure 5. IMRS responses are non-random and distributed across the weighted gene set. Body contains: Most contrasts shifted in the positive direction. Body contains: Observed ΔIMRSz values remained strongest in locked-anchor and primary acute validation groups. Body contains: Together, these analyses reduce the likelihood that IMRS behavior is explained by random label structure or a single dominant gene."
    revised_caption_or_subtitle = "Figure 5. IMRS score shifts are non-random and distributed across the weighted gene set. Replace body phrase with: Most contrasts showed positive delivery-minus-control ΔIMRSz. Replace body phrase with: Observed delivery-minus-control ΔIMRSz values remained strongest in locked-anchor and primary acute validation groups. Replace body phrase with: Together, these analyses support IMRS score stability and non-degeneracy, reducing the likelihood that observed delivery-minus-control ΔIMRSz patterns are explained by random label structure or a single dominant gene."
    where_to_apply_manually = "Manuscript caption source or assembled PowerPoint/Word caption block outside revised_plots_v6"
    notes = "caption_manual_update_needed; no editable caption source was found in v5/v6"
  },
  [ordered]@{
    figure_id = "Supplementary Figure S1"
    current_caption_or_subtitle = "Supplementary Figure S1. Baseline comparator immune signatures and context-shifted audits support the biological interpretation of IMRS. Body contains: Comparator immune signatures show related delivery-associated directionality across manuscript analysis groups. Body contains: Dataset-level ΔIMRSz summaries provide detailed validation context. Body contains: Attenuated IMRS responses are enriched in late or context-shifted biological categories. Body contains: These analyses contextualize IMRS biology but do not replace the frozen IMRS score or establish clinical prediction."
    revised_caption_or_subtitle = "Supplementary Figure S1. Comparator immune signatures and context-shifted audits contextualize the acute delivery-associated innate transcriptional response measured by IMRS. Replace body phrase with: Comparator immune signatures show related delivery-associated directionality across manuscript analysis groups but do not replace the frozen IMRS scoring framework. Replace body phrase with: Dataset-level delivery-minus-control ΔIMRSz summaries provide detailed validation context. Replace body phrase with: Attenuated IMRS responses are concentrated in late or context-shifted biological categories. Replace body phrase with: These analyses contextualize the acute delivery-associated innate transcriptional response measured by IMRS but do not replace the frozen IMRS scoring framework or establish clinical prediction."
    where_to_apply_manually = "Manuscript caption source or assembled PowerPoint/Word caption block outside revised_plots_v6"
    notes = "caption_manual_update_needed; no editable caption source was found in v5/v6"
  }
)
Write-TsvRows (Join-Path $Root "v6_manual_caption_replacement_table.tsv") $manualCaptionRows @("figure_id", "current_caption_or_subtitle", "revised_caption_or_subtitle", "where_to_apply_manually", "notes")

$changeRows = @(
  [ordered]@{ file = "scripts/00_generate_all_reorganized_figures_revised_v5.R"; panel_or_figure = "global"; change_type = "output_root_update"; old_text = "revised_plots_v5"; new_text = "revised_plots_v6"; plot_visible_or_caption = "script"; regenerated_or_copied = "script edited"; notes = "All writes now target v6." },
  [ordered]@{ file = "scripts/R/build_merged_imrs_workflow_v5.R"; panel_or_figure = "Figure1A"; change_type = "layout_text_fit"; old_text = "text overflow / clipping"; new_text = "adjusted font size, wrapping, box size, canvas, and margins"; plot_visible_or_caption = "plot_visible"; regenerated_or_copied = "regenerated"; notes = "Formatting-only Figure 1A fix; no data, calculations, workflow order, colors, or plotting logic changed." },
  [ordered]@{ file = "scripts/R/build_merged_imrs_workflow_v5.R"; panel_or_figure = "Figure1A"; change_type = "plot_visible_text"; old_text = "Frozen IMRS construction and split-contrast validation workflow"; new_text = "Frozen IMRS scoring and transfer-evaluation workflow"; plot_visible_or_caption = "plot_visible"; regenerated_or_copied = "regenerated by R"; notes = "Panel A title." },
  [ordered]@{ file = "scripts/R/build_merged_imrs_workflow_v5.R"; panel_or_figure = "Figure1A"; change_type = "plot_visible_text"; old_text = "Anchor-derived gene coefficients are frozen before scoring independent delivery-control contrasts"; new_text = "Anchor-derived gene coefficients are frozen before scoring independent delivery-versus-control contrasts"; plot_visible_or_caption = "plot_visible"; regenerated_or_copied = "regenerated by R"; notes = "Panel A subtitle." },
  [ordered]@{ file = "scripts/R/build_merged_imrs_workflow_v5.R"; panel_or_figure = "Figure1A"; change_type = "plot_visible_text"; old_text = "Anchor model construction"; new_text = "Locked-anchor model construction"; plot_visible_or_caption = "plot_visible"; regenerated_or_copied = "regenerated by R"; notes = "Panel A column heading." },
  [ordered]@{ file = "scripts/R/build_merged_imrs_workflow_v5.R"; panel_or_figure = "Figure1A"; change_type = "plot_visible_text"; old_text = "Anchor-set delivery-control differential expression"; new_text = "Locked-anchor delivery-versus-control differential expression"; plot_visible_or_caption = "plot_visible"; regenerated_or_copied = "regenerated by R"; notes = "Panel A workflow node." },
  [ordered]@{ file = "scripts/R/build_merged_imrs_workflow_v5.R"; panel_or_figure = "Figure1A"; change_type = "plot_visible_text"; old_text = "Split-level ΔIMRSz, directionality, and secondary AUC"; new_text = "Delivery-minus-control ΔIMRSz, directionality, and secondary AUC"; plot_visible_or_caption = "plot_visible"; regenerated_or_copied = "regenerated by R"; notes = "Panel A workflow node." },
  [ordered]@{ file = "scripts/R/build_merged_imrs_workflow_v5.R"; panel_or_figure = "Figure1A"; change_type = "plot_visible_text"; old_text = "Dataset-role curation and context-shifted dataset audit"; new_text = "Dataset-role curation and boundary-context audit"; plot_visible_or_caption = "plot_visible"; regenerated_or_copied = "regenerated by R"; notes = "Panel A workflow node." },
  [ordered]@{ file = "scripts/R/build_merged_imrs_workflow_v5.R"; panel_or_figure = "Figure1A"; change_type = "plot_visible_text"; old_text = "Null permutation, baseline benchmarking, and coefficient sensitivity"; new_text = "Label permutation, comparator benchmarking, and coefficient sensitivity"; plot_visible_or_caption = "plot_visible"; regenerated_or_copied = "regenerated by R"; notes = "Panel A workflow node." },
  [ordered]@{ file = "scripts/R/build_merged_imrs_workflow_v5.R"; panel_or_figure = "Figure1A"; change_type = "plot_visible_text"; old_text = "Frozen anchor-derived coefficients are not refit on validation or transfer datasets."; new_text = "Frozen anchor-derived coefficients are not refit during validation or transfer evaluation."; plot_visible_or_caption = "plot_visible"; regenerated_or_copied = "regenerated by R"; notes = "Panel A bottom note." },
  [ordered]@{ file = "scripts/R/v5_helpers.R"; panel_or_figure = "Figure1B"; change_type = "plot_visible_text"; old_text = "Mean delivery-minus-control ΔIMRSz, pseudo-log scale"; new_text = "Mean delivery-minus-control ΔIMRSz (pseudo-log scale)"; plot_visible_or_caption = "plot_visible"; regenerated_or_copied = "regenerated by R"; notes = "Panel B x-axis label." },
  [ordered]@{ file = "scripts/R/v5_helpers.R"; panel_or_figure = "Figure1B"; change_type = "plot_visible_text"; old_text = "Contrasts"; new_text = "Passing split contrasts"; plot_visible_or_caption = "plot_visible"; regenerated_or_copied = "regenerated by R"; notes = "Panel B size legend title." },
  [ordered]@{ file = "scripts/R/v5_helpers.R"; panel_or_figure = "Figure3B"; change_type = "plot_visible_text"; old_text = "Contrasts"; new_text = "Passing split contrasts"; plot_visible_or_caption = "plot_visible"; regenerated_or_copied = "regenerated by R"; notes = "Size legend title." },
  [ordered]@{ file = "scripts/R/v5_helpers.R"; panel_or_figure = "Figure5A"; change_type = "plot_visible_text"; old_text = "Split contrasts ordered by observed ΔIMRSz"; new_text = "Split contrasts ordered by observed delivery-minus-control ΔIMRSz"; plot_visible_or_caption = "plot_visible"; regenerated_or_copied = "regenerated by R"; notes = "Optional plot-visible improvement requested in module." }
)

foreach ($row in $manualCaptionRows) {
  $changeRows += [ordered]@{
    file = "v6_manual_caption_replacement_table.tsv"
    panel_or_figure = $row.figure_id
    change_type = "caption_manual_update_needed"
    old_text = $row.current_caption_or_subtitle
    new_text = $row.revised_caption_or_subtitle
    plot_visible_or_caption = "caption"
    regenerated_or_copied = "caption_manual_update_needed"
    notes = $row.notes
  }
}

$changeRows += @(
  [ordered]@{ file = "Figure2_main_v5.png/pdf/svg"; panel_or_figure = "Figure2"; change_type = "output_status"; old_text = ""; new_text = ""; plot_visible_or_caption = "output"; regenerated_or_copied = "copied unchanged"; notes = "No plot-visible text changes required." },
  [ordered]@{ file = "Figure4_main_v5.png/pdf/svg"; panel_or_figure = "Figure4"; change_type = "output_status"; old_text = ""; new_text = ""; plot_visible_or_caption = "output"; regenerated_or_copied = "copied unchanged"; notes = "Caption change only; no editable caption source present." },
  [ordered]@{ file = "FigureS_comparator_benchmarking_v5.png/pdf/svg"; panel_or_figure = "Supplementary Figure S1"; change_type = "output_status"; old_text = ""; new_text = ""; plot_visible_or_caption = "output"; regenerated_or_copied = "copied unchanged"; notes = "Caption change only; no editable caption source present." }
)
Write-TsvRows (Join-Path $Root "v6_change_log.tsv") $changeRows @("file", "panel_or_figure", "change_type", "old_text", "new_text", "plot_visible_or_caption", "regenerated_or_copied", "notes")

$panelManifest = Import-Csv -LiteralPath (Join-Path $Root "tables\v5_panel_manifest.tsv") -Delimiter "`t"
$figureManifest = Import-Csv -LiteralPath (Join-Path $Root "tables\v5_figure_manifest_wide.tsv") -Delimiter "`t"
$outputRows = @()
foreach ($row in $panelManifest) {
  $paths = @($row.output_png, $row.output_pdf, $row.output_svg)
  $regen = if ($RegeneratedPanels -contains $row.panel_id) { "regenerated" } else { "copied_unchanged" }
  $script = if ($regen -eq "regenerated") {
    if ($row.panel_id -eq "Figure1A") { "scripts/R/build_merged_imrs_workflow_v5.R" } else { "scripts/01_regenerate_changed_v6.R; scripts/R/v5_helpers.R" }
  } else {
    "copied unchanged from revised_plots_v5"
  }
  $outputRows += [ordered]@{
    figure_id = $row.figure_id
    panel_id = $row.panel_id
    output_png = $row.output_png
    output_pdf = $row.output_pdf
    output_svg = $row.output_svg
    source_script = $script
    regenerated_or_copied = $regen
    file_exists = ExistingNonzeroSummary $paths
    file_size_bytes = SizeSummary $paths
    notes = $row.notes
  }
}
foreach ($row in $figureManifest) {
  $paths = @($row.output_png, $row.output_pdf, $row.output_svg)
  $regen = if ($RegeneratedFigures -contains $row.figure_id) { "regenerated" } else { "copied_unchanged" }
  $outputRows += [ordered]@{
    figure_id = $row.figure_id
    panel_id = ""
    output_png = $row.output_png
    output_pdf = $row.output_pdf
    output_svg = $row.output_svg
    source_script = if ($regen -eq "regenerated") { "scripts/01_regenerate_changed_v6.R" } else { "copied unchanged from revised_plots_v5" }
    regenerated_or_copied = $regen
    file_exists = ExistingNonzeroSummary $paths
    file_size_bytes = SizeSummary $paths
    notes = if ($row.figure_id -eq "FigureS_comparator_benchmarking") { "Supplementary Figure S1 comparator figure copied unchanged; caption replacement listed manually." } else { $row.notes }
  }
}
Write-TsvRows (Join-Path $Root "v6_output_manifest.tsv") $outputRows @("figure_id", "panel_id", "output_png", "output_pdf", "output_svg", "source_script", "regenerated_or_copied", "file_exists", "file_size_bytes", "notes")

$scriptRows = @()
foreach ($file in (Get-ChildItem -LiteralPath $Root -Recurse -File | Where-Object { $_.Extension -in @(".R", ".ps1") } | Sort-Object FullName)) {
  $rel = RelPath $file.FullName
  $sourcePath = Join-Path $SourceRoot $rel
  $edited = (!(Test-Path -LiteralPath $sourcePath)) -or ((FileHashOrEmpty $file.FullName) -ne (FileHashOrEmpty $sourcePath))
  $text = [System.IO.File]::ReadAllText($file.FullName)
  $roots = [regex]::Matches($text, "revised_plots_v\d+") | ForEach-Object { $_.Value } | Sort-Object -Unique
  $scriptRows += [ordered]@{
    script_file = $rel
    purpose_inferred = PurposeForScript $rel
    edited_yes_no = if ($edited) { "yes" } else { "no" }
    key_changes = KeyChangesForScript $rel
    output_root_detected = ($roots -join ";")
    output_root_after_revision = if ($roots -contains "revised_plots_v6") { "revised_plots_v6" } else { "none detected" }
  }
}
Write-TsvRows (Join-Path $Root "v6_script_inventory.tsv") $scriptRows @("script_file", "purpose_inferred", "edited_yes_no", "key_changes", "output_root_detected", "output_root_after_revision")

$figureV5Log = @(
  "This copied v5 generation log was superseded inside revised_plots_v6.",
  "See v6_regeneration_log.txt and v6_regeneration_targeted_run.log for v6 actions.",
  "The frozen source folder revised_plots_v5 was not used as an output target."
)
[System.IO.File]::WriteAllLines((Join-Path $Root "figure_v5_generation_log.txt"), $figureV5Log, [System.Text.UTF8Encoding]::new($false))

$validation = New-Object System.Collections.Generic.List[string]
$validation.Add("revised_plots_v6 exists: $(Test-Path -LiteralPath $Root)")
$validation.Add("Figure 1 Panel A title exact: $(TextContains (Join-Path $Root 'Figure1A_IMRS_merged_workflow_v5.svg') 'Frozen IMRS scoring and transfer-evaluation workflow')")
$validation.Add("Figure 1 Panel A subtitle exact: $(TextContains (Join-Path $Root 'Figure1A_IMRS_merged_workflow_v5.svg') 'Anchor-derived gene coefficients are frozen before scoring independent delivery-versus-control contrasts')")
$validation.Add("Figure 1 Panel B x-axis exact: $(TextContains (Join-Path $Root 'Figure1B_dataset_tissue_response_landscape_v5_corrected.svg') 'Mean delivery-minus-control ΔIMRSz (pseudo-log scale)')")
$validation.Add("Figure 1 size legend exact: $(TextContains (Join-Path $Root 'Figure1B_dataset_tissue_response_landscape_v5_corrected.svg') 'Passing split contrasts')")
$validation.Add("Figure 3 size legend exact: $(TextContains (Join-Path $Root 'intermediate_panels\Figure3B_v5.svg') 'Passing split contrasts')")
$validation.Add("Figure 5 x-axis updated: $(TextContains (Join-Path $Root 'intermediate_panels\Figure5A_v5.svg') 'Split contrasts ordered by observed delivery-minus-control ΔIMRSz')")
$validation.Add("Manual caption table exists: $(Test-Path -LiteralPath (Join-Path $Root 'v6_manual_caption_replacement_table.tsv'))")
$badManifestOutputs = @()
foreach ($row in $outputRows) {
  foreach ($path in @($row.output_png, $row.output_pdf, $row.output_svg)) {
    if (!(Test-Path -LiteralPath $path) -or ((Get-Item -LiteralPath $path).Length -le 0)) {
      $badManifestOutputs += $path
    }
  }
}
$validation.Add("Every output manifest file exists and is nonzero: $($badManifestOutputs.Count -eq 0)")
$v5HashCheck = CompareV5Hashes
$validation.Add("revised_plots_v5 hash check: $($v5HashCheck.status) - $($v5HashCheck.detail)")

$regenRows = $outputRows | Where-Object { $_.regenerated_or_copied -eq "regenerated" }
$copyRows = $outputRows | Where-Object { $_.regenerated_or_copied -eq "copied_unchanged" }
$validationLines = $validation | ForEach-Object { "  $_" }
$logLines = @(
  "start time: $(Get-Date -Format 'yyyy-MM-dd HH:mm:ss')",
  "source folder: $SourceRoot",
  "output folder: $Root",
  "scripts run:",
  "  Rscript scripts/01_regenerate_changed_v6.R",
  "  Rscript -e <regenerate v5/v6 wording audit>",
  "  powershell scripts/02_create_v6_metadata.ps1",
  "warnings:",
  "  Rscript was not on PATH; used full installed Rscript path.",
  "errors:",
  "  none",
  "files generated:",
  "  Regenerated panels: $($RegeneratedPanels -join ', ')",
  "  Regenerated final figures: $($RegeneratedFigures -join ', ')",
  "  Metadata files: $($MetadataFiles -join ', ')",
  "files copied unchanged:",
  "  Output manifest rows copied unchanged: $($copyRows.Count)",
  "  Copied unchanged figures include Figure2, Figure4, Supplementary Figure S1 comparator outputs, and non-target supplemental outputs.",
  "validation summary:"
) + $validationLines + @(
  "notes:",
  "  No statistical calculations, sample inclusion, grouping, geoms, colors, scales, or panel structures were intentionally changed.",
  "  Captions were not found in editable R/Markdown/Word/PowerPoint caption sources in v5/v6, so exact replacements are listed in v6_manual_caption_replacement_table.tsv."
)
[System.IO.File]::WriteAllLines((Join-Path $Root "v6_regeneration_log.txt"), $logLines, [System.Text.UTF8Encoding]::new($false))

$phrases = @(
  "acute delivery-response",
  "acute innate-response",
  "acute delivery-associated gene program",
  "acute delivery-associated transcriptional program",
  "innate immune transcriptional activity",
  "innate immune activation",
  "IMRS elevation",
  "fixed score",
  "frozen IMRS score",
  "acute delivery-associated contexts",
  "weak datasets",
  "weak-context audit",
  "failure pattern",
  "IMRS responses are non-random",
  "baseline benchmarking",
  "null permutation"
)
$textFiles = Get-ChildItem -LiteralPath $Root -Recurse -File | Where-Object {
  $_.Extension -notin @(".png", ".pdf", ".svg") -and
    $_.Name -ne "v6_terminology_audit.tsv" -and
    (RelPath $_.FullName) -ne "scripts\02_create_v6_metadata.ps1"
}
$auditRows = @()
foreach ($phrase in $phrases) {
  $count = 0
  $files = New-Object System.Collections.Generic.List[string]
  foreach ($file in $textFiles) {
    $text = [System.IO.File]::ReadAllText($file.FullName)
    $matches = [regex]::Matches($text, [regex]::Escape($phrase))
    if ($matches.Count -gt 0) {
      $count += $matches.Count
      $files.Add((RelPath $file.FullName))
    }
  }
  $filesUnique = $files | Sort-Object -Unique
  $notes = if ($count -eq 0) {
    "No remaining text-file occurrences found."
  } else {
    "Remaining occurrences are replacement-history old_text, manual current-caption source text, or validation-log provenance; plot-visible SVG outputs and revised caption replacements were separately checked."
  }
  $auditRows += [ordered]@{
    searched_phrase = $phrase
    remaining_count_after_revision = $count
    remaining_files = ($filesUnique -join ";")
    action_needed = if ($count -eq 0) { "none" } else { "review only; justified historical/manual-replacement context" }
    notes = $notes
  }
}
Write-TsvRows (Join-Path $Root "v6_terminology_audit.tsv") $auditRows @("searched_phrase", "remaining_count_after_revision", "remaining_files", "action_needed", "notes")

$tableRows = @()
foreach ($file in Get-ChildItem -LiteralPath $Root -Recurse -File -Filter *.tsv | Sort-Object FullName) {
  $rel = RelPath $file.FullName
  $sourcePath = Join-Path $SourceRoot $rel
  $edited = (!(Test-Path -LiteralPath $sourcePath)) -or ((FileHashOrEmpty $file.FullName) -ne (FileHashOrEmpty $sourcePath))
  $shape = TsvShape $file.FullName
  $tableRows += [ordered]@{
    table_file = $rel
    purpose_inferred = PurposeForTable $rel
    edited_yes_no = if ($edited) { "yes" } else { "no" }
    key_changes = KeyChangesForTable $rel $edited
    row_count = $shape.rows
    column_count = $shape.columns
  }
}
Write-TsvRows (Join-Path $Root "v6_table_inventory.tsv") $tableRows @("table_file", "purpose_inferred", "edited_yes_no", "key_changes", "row_count", "column_count")

$fileRows = @()
foreach ($file in Get-ChildItem -LiteralPath $Root -Recurse -File | Sort-Object FullName) {
  $rel = RelPath $file.FullName
  $sourcePath = Join-Path $SourceRoot $rel
  $sourceExists = Test-Path -LiteralPath $sourcePath
  $generated = (!$sourceExists) -or ((FileHashOrEmpty $file.FullName) -ne (FileHashOrEmpty $sourcePath))
  $fileRows += [ordered]@{
    file = $rel
    extension = $file.Extension
    size_bytes = $file.Length
    modified_time = $file.LastWriteTime.ToString("yyyy-MM-dd HH:mm:ss")
    copied_from_v5_yes_no = if ($sourceExists) { "yes" } else { "no" }
    generated_in_v6_yes_no = if ($generated) { "yes" } else { "no" }
  }
}
Write-TsvRows (Join-Path $Root "v6_file_inventory.tsv") $fileRows @("file", "extension", "size_bytes", "modified_time", "copied_from_v5_yes_no", "generated_in_v6_yes_no")

Write-Host "Created v6 metadata outputs in $Root"
