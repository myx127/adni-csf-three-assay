$ErrorActionPreference = 'Stop'
& python -X utf8 (Join-Path $PSScriptRoot 'verify_delivery.py')
exit $LASTEXITCODE
