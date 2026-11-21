$ErrorActionPreference = 'Stop'
$env:PYTHONIOENCODING = 'utf-8'
function Invoke-CalibrationStep {
    param([string]$ScriptName, [string[]]$StepArguments = @())
    & python (Join-Path $PSScriptRoot $ScriptName) @StepArguments
    if ($LASTEXITCODE -ne 0) { throw "Calibration step failed: $ScriptName $StepArguments" }
}
Invoke-CalibrationStep 'validate_preflight.py'
Invoke-CalibrationStep 'validate_mean_preflight.py'
Invoke-CalibrationStep 'run_calibration.py' @('--phase','DEVELOPMENT')
Invoke-CalibrationStep 'analyze_calibration.py' @('--select-development')
Invoke-CalibrationStep 'run_calibration.py' @('--phase','CONFIRMATION')
Invoke-CalibrationStep 'run_calibration.py' @('--phase','N_GRID')
Invoke-CalibrationStep 'run_calibration.py' @('--phase','SENSITIVITY')
Invoke-CalibrationStep 'mean_followup.py' @('--phase','MEAN_DEV')
Invoke-CalibrationStep 'mean_followup.py' @('--select')
Invoke-CalibrationStep 'mean_followup.py' @('--phase','MEAN_CONFIRM')
Invoke-CalibrationStep 'mean_followup.py' @('--phase','MEAN_SENS_REPLAY')
Invoke-CalibrationStep 'bootstrap_resolution.py'
Invoke-CalibrationStep 'finish_analysis.py'
Invoke-CalibrationStep 'validate_completed.py'
Invoke-CalibrationStep 'build_delivery.py'
Invoke-CalibrationStep 'finalize_delivery.py'
