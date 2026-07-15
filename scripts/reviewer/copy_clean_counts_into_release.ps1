$ErrorActionPreference = "Stop"

$releaseRoot = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot "..\..")).Path
$sourceRoot = "D:\IMRS_Project\03_counts"
$datasets = @(
    "GSE119119"
    "GSE139529"
    "GSE166655"
    "GSE167521"
    "GSE178313"
    "GSE262515"
    "GSE264344"
    "GSE279372"
    "GSE279743"
    "GSE279744"
    "GSE314070"
    "GSE39129"
)

$mapping = $datasets | ForEach-Object {
    [PSCustomObject]@{
        Dataset = $_
        Source = Join-Path $sourceRoot "$_\featurecounts\validation\gene_counts_clean.tsv"
        Target = Join-Path $releaseRoot "data\counts\$_\featurecounts\validation\gene_counts_clean.tsv"
    }
}

$missingSources = $mapping | Where-Object { -not (Test-Path -LiteralPath $_.Source -PathType Leaf) }
if ($missingSources) {
    $missingList = ($missingSources.Source | ForEach-Object { " - $_" }) -join [Environment]::NewLine
    throw "Missing clean count source file(s):$([Environment]::NewLine)$missingList"
}

Write-Host "Copying cleaned featureCounts matrices into $releaseRoot"
foreach ($item in $mapping) {
    $targetDirectory = Split-Path -Parent $item.Target
    New-Item -ItemType Directory -Force -Path $targetDirectory | Out-Null
    Copy-Item -LiteralPath $item.Source -Destination $item.Target -Force
    $targetFile = Get-Item -LiteralPath $item.Target
    Write-Host ("{0}`t{1:N2} MB`t{2}" -f $item.Dataset, ($targetFile.Length / 1MB), $targetFile.FullName)
}

$missingTargets = $mapping | Where-Object { -not (Test-Path -LiteralPath $_.Target -PathType Leaf) }
if ($missingTargets) {
    $missingList = ($missingTargets.Target | ForEach-Object { " - $_" }) -join [Environment]::NewLine
    throw "Copy completed with missing expected target file(s):$([Environment]::NewLine)$missingList"
}

Write-Host ("Confirmed all {0} expected cleaned count matrices are present." -f $mapping.Count)
