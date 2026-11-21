$ErrorActionPreference = 'Stop'
$stagePath = Split-Path -Parent $PSScriptRoot
& python -X utf8 (Join-Path $PSScriptRoot 'verify_delivery.py')
exit $LASTEXITCODE
