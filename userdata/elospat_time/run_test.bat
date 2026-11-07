@echo off
REM run_test.bat
REM Helper script to run MATLAB tests for RMS-ELO-SPAT
REM
REM Usage:
REM   run_test.bat task5
REM   run_test.bat task6
REM   run_test.bat task7

setlocal

set MATLAB_EXE="C:\Program Files\MATLAB\R2026a\bin\matlab.exe"
set MATRAD_PATH=D:\matRad
set WORKSPACE_PATH=d:\Elo-SPAT\userdata\elospat_time

if "%1"=="" (
    echo Error: Test name is required
    echo Usage: run_test.bat [task5^|task6^|task7]
    exit /b 1
)

if "%1"=="task5" (
    set "MATLAB_CMD=cd('d:\Elo-SPAT\userdata\elospat_time\planning'); addpath('D:\matRad'); matRad_rc; test_build_impt_plan"
) else if "%1"=="task6" (
    set "MATLAB_CMD=cd('d:\Elo-SPAT\userdata\elospat_time\planning'); addpath('D:\matRad'); matRad_rc; test_build_discrete_arc_plan"
) else if "%1"=="task7" (
    set "MATLAB_CMD=cd('d:\Elo-SPAT\userdata\elospat_time\planning'); addpath('D:\matRad'); matRad_rc; test_run_dense_arc_baseline"
) else (
    echo Error: Unknown test name '%1'
    echo Available tests: task5, task6, task7
    exit /b 1
)

echo ========================================
echo RMS-ELO-SPAT Test Runner
echo ========================================
echo Test: %1
echo MATLAB: %MATLAB_EXE%
echo matRad: %MATRAD_PATH%
echo ========================================
echo.

%MATLAB_EXE% -batch "%MATLAB_CMD%"

if %ERRORLEVEL% EQU 0 (
    echo.
    echo ========================================
    echo Test completed successfully!
    echo ========================================
) else (
    echo.
    echo ========================================
    echo Test failed with exit code: %ERRORLEVEL%
    echo ========================================
    exit /b %ERRORLEVEL%
)

endlocal
