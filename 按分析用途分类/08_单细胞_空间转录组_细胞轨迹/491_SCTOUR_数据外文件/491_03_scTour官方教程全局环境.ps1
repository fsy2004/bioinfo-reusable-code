# UTF-8 PowerShell script: run scTour with the global/user Python 3.12 environment.
param(
    [string]$InputFile = "",
    [string]$OutputDir = "",
    [int]$NTopGenes = 1000,
    [int]$NeEpoch = 0,
    [double]$Percent = -1,
    [switch]$ForceCpu
)

$ErrorActionPreference = "Stop"
[Console]::OutputEncoding = [System.Text.Encoding]::UTF8
$OutputEncoding = [System.Text.Encoding]::UTF8

$Root = Split-Path -Parent $MyInvocation.MyCommand.Path
if ([string]::IsNullOrWhiteSpace($InputFile)) {
    $InputFile = Join-Path $Root "EX_development_human_cortex_10X.h5ad"
}
if ([string]::IsNullOrWhiteSpace($OutputDir)) {
    $OutputDir = Join-Path $Root "sctour_results_official_tutorial_exact"
}

if (-not (Test-Path -LiteralPath $InputFile)) {
    throw "Input h5ad not found: $InputFile"
}

$PyArgs = @(
    "-3.12",
    (Join-Path $Root "06_run_official_tutorial_exact.py"),
    "--input", $InputFile,
    "--output-dir", $OutputDir,
    "--n-top-genes", "$NTopGenes"
)
if ($NeEpoch -gt 0) {
    $PyArgs += @("--nepoch", "$NeEpoch")
}
if ($Percent -ge 0) {
    $PyArgs += @("--percent", "$Percent")
}
if ($ForceCpu) {
    $PyArgs += "--force-cpu"
}

& py @PyArgs