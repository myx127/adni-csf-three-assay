param()
$ErrorActionPreference = 'Stop'
$stage10Verifier = Join-Path $PSScriptRoot 'verify_delivery.py'
& python -X utf8 $stage10Verifier
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
