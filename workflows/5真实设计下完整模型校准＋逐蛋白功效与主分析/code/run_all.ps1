$ErrorActionPreference = 'Stop'
$env:PYTHONIOENCODING = 'utf-8'
$auditRoot = Split-Path -Parent $PSScriptRoot
$auditLogDirectory = Join-Path $auditRoot 'logs'
$jobs = @(
    @('verify_sources.py'),
    @('validate_preflight.py'),
    @('validate_diagnostics.py'),
    @('verify_reproducibility.py'),
    @('load_proteins.py'),
    @('analyze_main.py'),
    @('summarize_power.py'),
    @('post_analysis_diagnostics.py'),
    @('build_delivery.py','--tables-only'),
    @('validate_completed.py'),
    @('build_delivery.py'),
    @('finalize_delivery.py')
)
foreach ($job in $jobs) {
    $script = $job[0]
    Write-Host ('Running ' + $script)
    $startInfo = New-Object System.Diagnostics.ProcessStartInfo
    $startInfo.FileName = 'python'
    $startInfo.Arguments = '"' + (Join-Path $PSScriptRoot $script) + '"'
    if ($job.Count -gt 1) { $startInfo.Arguments += ' ' + $job[1] }
    $startInfo.UseShellExecute = $false
    $startInfo.CreateNoWindow = $true
    $startInfo.RedirectStandardOutput = $true
    $startInfo.RedirectStandardError = $true
    $process = New-Object System.Diagnostics.Process
    $process.StartInfo = $startInfo
    $process.Start() | Out-Null
    $stdoutTask = $process.StandardOutput.ReadToEndAsync()
    $stderrTask = $process.StandardError.ReadToEndAsync()
    $process.WaitForExit()
    $logText = $stdoutTask.Result + [Environment]::NewLine + $stderrTask.Result
    if ($script -eq 'finalize_delivery.py') { Write-Host $logText }
    else { [IO.File]::WriteAllText((Join-Path $auditLogDirectory ([IO.Path]::GetFileNameWithoutExtension($script) + '_REPRO.log')),$logText,[Text.UTF8Encoding]::new($false)) }
    if ($process.ExitCode -ne 0) { throw ('Failed ' + $script + '; see logs') }
}
Write-Host 'Frozen main analysis reproduced, verified and sealed.'
