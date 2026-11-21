$ErrorActionPreference = 'Stop'
$env:PYTHONIOENCODING = 'utf-8'
$auditRoot = Split-Path -Parent $PSScriptRoot
$auditLogs = Join-Path $auditRoot 'logs'
New-Item -ItemType Directory -Path $auditLogs -Force | Out-Null
$scripts = @('inspect_inputs.py','inspect_methods.py','inspect_metadata.py','inspect_formats.py','run_coverage.py','protein_coverage.py','supplement_coverage.py','prepare_tables.py','validate_independently.py','build_deliverables.py','finalize_audit.py')
foreach ($script in $scripts) {
    Write-Host ('Running ' + $script)
    $logName = [IO.Path]::GetFileNameWithoutExtension($script) + '_REPRO.log'
    $logPath = Join-Path $auditLogs $logName
    # Capture Python stderr without PowerShell treating nonfatal warnings as failures.
    $startInfo = New-Object System.Diagnostics.ProcessStartInfo
    $startInfo.FileName = 'python'
    $startInfo.Arguments = '"' + (Join-Path $PSScriptRoot $script) + '"'
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
    if ($script -eq 'finalize_audit.py') { Write-Host $logText }
    else { [IO.File]::WriteAllText($logPath,$logText,[Text.UTF8Encoding]::new($false)) }
    if ($process.ExitCode -ne 0) { throw ($script + ' failed; see ' + $logPath) }
}
Write-Host 'Stage4 audit completed and output hashes sealed.'
