$ErrorActionPreference = 'Stop'
$env:PYTHONIOENCODING = 'utf-8'

foreach ($auditScript in @('run_audit.py', 'validate_independently.py', 'build_deliverables.py', 'finalize_audit.py')) {
    $auditScriptPath = Join-Path -Path $PSScriptRoot -ChildPath $auditScript
    & python $auditScriptPath
    if ($LASTEXITCODE -ne 0) {
        throw "Audit step failed: $auditScript"
    }
}
