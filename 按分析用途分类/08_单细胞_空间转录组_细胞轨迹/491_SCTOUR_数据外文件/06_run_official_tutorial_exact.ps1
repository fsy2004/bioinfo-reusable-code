# UTF-8 PowerShell script for reproducing the official scTour tutorial figure.
param(
    [int]$NTopGenes = 1000,
    [int]$NeEpoch = 0,
    [double]$Percent = -1,
    [switch]$ForceCpu
)

$ErrorActionPreference = "Stop"
[Console]::OutputEncoding = [System.Text.Encoding]::UTF8
$OutputEncoding = [System.Text.Encoding]::UTF8

$ScriptPath = $MyInvocation.MyCommand.Path
if ([string]::IsNullOrWhiteSpace($ScriptPath)) {
    $Root = (Get-Location).Path
}
else {
    $Root = Split-Path -Parent $ScriptPath
}
$VenvPython = Join-Path $Root ".venv_sctour\Scripts\python.exe"
$InputFile = Join-Path $Root "EX_development_human_cortex_10X.h5ad"

if (-not (Test-Path -LiteralPath $VenvPython)) {
    throw "Python environment not found. Run first: powershell -ExecutionPolicy Bypass -File .\00_install_sctour.ps1"
}
if (-not (Test-Path -LiteralPath $InputFile)) {
    throw "Input h5ad not found: $InputFile"
}

Push-Location $Root
try {
    $Args = @(
        ".\06_run_official_tutorial_exact.py",
        "--input", $InputFile,
        "--n-top-genes", "$NTopGenes"
    )
    if ($NeEpoch -gt 0) {
        $Args += @("--nepoch", "$NeEpoch")
    }
    if ($Percent -ge 0) {
        $Args += @("--percent", "$Percent")
    }
    if ($ForceCpu) {
        $Args += "--force-cpu"
    }

    & $VenvPython @Args
}
finally {
    Pop-Location
}
