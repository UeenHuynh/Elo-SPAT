# run_test.ps1
# Helper script to run MATLAB tests for RMS-ELO-SPAT
#
# Usage:
#   .\run_test.ps1 task5
#   .\run_test.ps1 task6
#   .\run_test.ps1 custom "your_matlab_code_here"

param(
    [Parameter(Mandatory=$true)]
    [string]$TestName,
    
    [Parameter(Mandatory=$false)]
    [string]$CustomCode = ""
)

# MATLAB configuration
$matlabExe = "C:\Program Files\MATLAB\R2026a\bin\matlab.exe"
$matRadPath = "D:\matRad"
$workspacePath = "d:\Elo-SPAT\userdata\elospat_time"

# Test definitions
$tests = @{
    "task5" = "cd('$workspacePath\planning'); addpath('$matRadPath'); matRad_rc; test_build_impt_plan"
    "task6" = "cd('$workspacePath\planning'); addpath('$matRadPath'); matRad_rc; test_build_discrete_arc_plan"
}

# Select test command
if ($TestName -eq "custom") {
    if ([string]::IsNullOrEmpty($CustomCode)) {
        Write-Host "Error: Custom code is required when using 'custom' test name" -ForegroundColor Red
        exit 1
    }
    $matlabCommand = "addpath('$matRadPath'); matRad_rc; $CustomCode"
} elseif ($tests.ContainsKey($TestName)) {
    $matlabCommand = $tests[$TestName]
} else {
    Write-Host "Error: Unknown test name '$TestName'" -ForegroundColor Red
    Write-Host "Available tests: $($tests.Keys -join ', '), custom" -ForegroundColor Yellow
    exit 1
}

# Display info
Write-Host "========================================" -ForegroundColor Cyan
Write-Host "RMS-ELO-SPAT Test Runner" -ForegroundColor Cyan
Write-Host "========================================" -ForegroundColor Cyan
Write-Host "Test: $TestName" -ForegroundColor Green
Write-Host "MATLAB: $matlabExe" -ForegroundColor Gray
Write-Host "matRad: $matRadPath" -ForegroundColor Gray
Write-Host "Command: $matlabCommand" -ForegroundColor Gray
Write-Host "========================================`n" -ForegroundColor Cyan

# Run MATLAB
& $matlabExe -batch $matlabCommand

# Check exit code
if ($LASTEXITCODE -eq 0) {
    Write-Host "`n========================================" -ForegroundColor Cyan
    Write-Host "Test completed successfully!" -ForegroundColor Green
    Write-Host "========================================" -ForegroundColor Cyan
} else {
    Write-Host "`n========================================" -ForegroundColor Cyan
    Write-Host "Test failed with exit code: $LASTEXITCODE" -ForegroundColor Red
    Write-Host "========================================" -ForegroundColor Cyan
    exit $LASTEXITCODE
}
