param(
    [switch]$RecomputeResults,
    [switch]$RunSimulations
)
$ErrorActionPreference = 'Stop'
if ($RunSimulations -and -not $RecomputeResults) {
    throw 'RunSimulations requires RecomputeResults in a sibling reproduction copy.'
}
if ($RecomputeResults) {
    $ReproductionScript = Join-Path $PSScriptRoot 'reproduce_results.py'
    if ($RunSimulations) {
        & python -X utf8 $ReproductionScript --simulate
    } else {
        & python -X utf8 $ReproductionScript
    }
} else {
    $VerificationScript = Join-Path $PSScriptRoot 'verify_delivery.py'
    & python -X utf8 $VerificationScript
}
if ($LASTEXITCODE -ne 0) {
    throw "Stage6 verification/reproduction failed with exit code $LASTEXITCODE"
}
