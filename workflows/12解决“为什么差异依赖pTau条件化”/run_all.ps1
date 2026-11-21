$ErrorActionPreference = 'Stop'
$verifyScript = Join-Path $PSScriptRoot 'code\verify_delivery.py'
& python -B $verifyScript @args
exit $LASTEXITCODE
