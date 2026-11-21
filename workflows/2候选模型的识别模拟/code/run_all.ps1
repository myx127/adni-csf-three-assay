$ErrorActionPreference = 'Stop'
$env:PYTHONIOENCODING = 'utf-8'
foreach ($simulationStep in @('run_simulations.py','summarize_simulations.py','validate_simulations.py','parameter_diagnostics.py','build_simulation_delivery.py','finalize_delivery.py')) {
    $simulationStepPath = Join-Path -Path $PSScriptRoot -ChildPath $simulationStep
    & python $simulationStepPath
    if ($LASTEXITCODE -ne 0) { throw "Simulation step failed: $simulationStep" }
}
