param()
$ErrorActionPreference = 'Stop'
$stage9VerifyScript = Join-Path $PSScriptRoot 'verify_delivery.py'
& python -X utf8 $stage9VerifyScript
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
