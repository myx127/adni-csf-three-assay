param([switch]$VerifyOnly)
$ErrorActionPreference = 'Stop'
$stageRoot = Split-Path $PSScriptRoot -Parent
if ($VerifyOnly -or (Test-Path -LiteralPath (Join-Path $stageRoot 'SHA256_MANIFEST.csv'))) {
    python -B (Join-Path $PSScriptRoot 'verify_delivery.py')
    if ($LASTEXITCODE -ne 0) { throw 'Verification failed.' }
    exit
}
$steps = @('bootstrap','prepare_templates','validate_engine','run_calibration','analyze_real','summarize_and_describe','validate_science','make_figures','build_delivery')
foreach ($stepName in $steps) {
    python -B (Join-Path $PSScriptRoot ($stepName + '.py'))
    if ($LASTEXITCODE -ne 0) { throw ('Failed step: ' + $stepName) }
}
if (-not (Test-Path -LiteralPath (Join-Path $stageRoot 'results\FIGURE_REVIEW_RECORD.csv'))) {
    Write-Output 'Numerical analysis complete. Inspect the eight figures and record the visual review before delivery validation and sealing.'
    exit
}
python -B (Join-Path $PSScriptRoot 'validate_delivery.py')
if ($LASTEXITCODE -ne 0) { throw 'Delivery validation failed.' }
python -B (Join-Path $PSScriptRoot 'seal_delivery.py')
if ($LASTEXITCODE -ne 0) { throw 'Sealing failed.' }

