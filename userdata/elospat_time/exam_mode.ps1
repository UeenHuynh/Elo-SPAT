param(
    [ValidateSet("enter", "restore", "status")]
    [string]$Mode = "status",

    [ValidateRange(1, 24)]
    [int]$ResumeAfterHours = 5
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

$scriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$statePath = Join-Path $scriptDir ".exam_mode_state.json"

$browserProcesses = @("chrome", "msedge")
$managedServices = @(
    [pscustomobject]@{
        Name = "sshd"
        ResumeStartupType = "Manual"
        TaskName = "EloSpat-Resume-sshd"
        ExtraProcesses = @()
    },
    [pscustomobject]@{
        Name = "Tailscale"
        ResumeStartupType = "Automatic"
        TaskName = "EloSpat-Resume-Tailscale"
        ExtraProcesses = @("tailscale", "tailscaled")
    }
)

function Test-IsAdministrator {
    $identity = [Security.Principal.WindowsIdentity]::GetCurrent()
    $principal = New-Object Security.Principal.WindowsPrincipal($identity)
    return $principal.IsInRole([Security.Principal.WindowsBuiltInRole]::Administrator)
}

function Assert-Administrator {
    if (-not (Test-IsAdministrator)) {
        throw "Run this script from an elevated PowerShell window (Run as Administrator)."
    }
}

function Get-ServiceRecord {
    param(
        [Parameter(Mandatory = $true)]
        [string]$Name
    )

    return Get-CimInstance -ClassName Win32_Service -Filter "Name='$Name'" -ErrorAction SilentlyContinue
}

function Get-StartupTypeForSetService {
    param(
        [AllowNull()]
        [string]$StartMode,

        [Parameter(Mandatory = $true)]
        [string]$Fallback
    )

    switch ($StartMode) {
        "Auto" { return "Automatic" }
        "Manual" { return "Manual" }
        "Disabled" { return "Disabled" }
        default { return $Fallback }
    }
}

function Load-State {
    if (-not (Test-Path -LiteralPath $statePath)) {
        return $null
    }

    return Get-Content -LiteralPath $statePath -Raw | ConvertFrom-Json
}

function Save-State {
    param(
        [Parameter(Mandatory = $true)]
        [pscustomobject]$State
    )

    $State | ConvertTo-Json -Depth 5 | Set-Content -LiteralPath $statePath -Encoding ASCII
}

function Remove-StateFile {
    if (Test-Path -LiteralPath $statePath) {
        Remove-Item -LiteralPath $statePath -Force
    }
}

function Get-ResumeTask {
    param(
        [Parameter(Mandatory = $true)]
        [string]$TaskName
    )

    return Get-ScheduledTask -TaskName $TaskName -ErrorAction SilentlyContinue
}

function Remove-ResumeTask {
    param(
        [Parameter(Mandatory = $true)]
        [string]$TaskName
    )

    if (Get-ResumeTask -TaskName $TaskName) {
        Unregister-ScheduledTask -TaskName $TaskName -Confirm:$false
    }
}

function Register-ResumeTask {
    param(
        [Parameter(Mandatory = $true)]
        [pscustomobject]$ServiceConfig,

        [Parameter(Mandatory = $true)]
        [datetime]$RunAt,

        [Parameter(Mandatory = $true)]
        [string]$StartupType
    )

    Remove-ResumeTask -TaskName $ServiceConfig.TaskName

    $command = "try { Set-Service -Name '$($ServiceConfig.Name)' -StartupType $StartupType -ErrorAction Stop; Start-Service -Name '$($ServiceConfig.Name)' -ErrorAction Stop } catch { exit 1 }"
    $action = New-ScheduledTaskAction -Execute "powershell.exe" -Argument "-NoProfile -WindowStyle Hidden -Command `"$command`""
    $trigger = New-ScheduledTaskTrigger -Once -At $RunAt
    $principal = New-ScheduledTaskPrincipal -UserId "SYSTEM" -LogonType ServiceAccount -RunLevel Highest

    Register-ScheduledTask `
        -TaskName $ServiceConfig.TaskName `
        -Description "Restore $($ServiceConfig.Name) after exam mode" `
        -Action $action `
        -Trigger $trigger `
        -Principal $principal `
        -Force | Out-Null
}

function Stop-BrowserProcesses {
    foreach ($processName in $browserProcesses) {
        $processes = @(Get-Process -Name $processName -ErrorAction SilentlyContinue)
        if ($processes.Count -gt 0) {
            $processes | Stop-Process -Force
        }
    }
}

function Stop-ExtraProcesses {
    param(
        [Parameter(Mandatory = $true)]
        [string[]]$ProcessNames
    )

    foreach ($processName in $ProcessNames) {
        $processes = @(Get-Process -Name $processName -ErrorAction SilentlyContinue)
        if ($processes.Count -gt 0) {
            $processes | Stop-Process -Force
        }
    }
}

function Enter-ExamMode {
    Assert-Administrator

    $resumeAt = (Get-Date).AddHours($ResumeAfterHours)
    $serviceState = @()

    foreach ($serviceConfig in $managedServices) {
        $record = Get-ServiceRecord -Name $serviceConfig.Name
        if (-not $record) {
            Write-Warning "Service '$($serviceConfig.Name)' was not found. Skipping."
            continue
        }

        $startupType = Get-StartupTypeForSetService -StartMode $record.StartMode -Fallback $serviceConfig.ResumeStartupType
        Register-ResumeTask -ServiceConfig $serviceConfig -RunAt $resumeAt -StartupType $startupType

        if ($record.State -ne "Stopped") {
            Stop-Service -Name $serviceConfig.Name -Force -ErrorAction SilentlyContinue
        }

        Set-Service -Name $serviceConfig.Name -StartupType Disabled

        if ($serviceConfig.ExtraProcesses.Count -gt 0) {
            Stop-ExtraProcesses -ProcessNames $serviceConfig.ExtraProcesses
        }

        $serviceState += [pscustomobject]@{
            Name = $serviceConfig.Name
            ResumeStartupType = $startupType
            TaskName = $serviceConfig.TaskName
        }
    }

    Stop-BrowserProcesses

    Save-State -State ([pscustomobject]@{
        Mode = "enter"
        EnteredAt = (Get-Date).ToString("o")
        ResumeAt = $resumeAt.ToString("o")
        BrowsersClosed = $browserProcesses
        Services = $serviceState
    })

    Write-Host ""
    Write-Host "Exam mode enabled." -ForegroundColor Green
    Write-Host "Resume time: $($resumeAt.ToString("yyyy-MM-dd HH:mm:ss"))"
    Write-Host ""
    Show-ExamModeStatus
}

function Restore-ExamMode {
    Assert-Administrator

    $savedState = Load-State
    $startupTypeByService = @{}

    if ($savedState -and $savedState.Services) {
        foreach ($serviceEntry in $savedState.Services) {
            $startupTypeByService[$serviceEntry.Name] = $serviceEntry.ResumeStartupType
        }
    }

    foreach ($serviceConfig in $managedServices) {
        $record = Get-ServiceRecord -Name $serviceConfig.Name
        if (-not $record) {
            continue
        }

        $startupType = $serviceConfig.ResumeStartupType
        if ($startupTypeByService.ContainsKey($serviceConfig.Name)) {
            $startupType = $startupTypeByService[$serviceConfig.Name]
        }

        Remove-ResumeTask -TaskName $serviceConfig.TaskName
        Set-Service -Name $serviceConfig.Name -StartupType $startupType
        Start-Service -Name $serviceConfig.Name -ErrorAction SilentlyContinue
    }

    Remove-StateFile

    Write-Host ""
    Write-Host "Exam mode restored." -ForegroundColor Green
    Write-Host ""
    Show-ExamModeStatus
}

function Show-ExamModeStatus {
    $browserStatus = foreach ($processName in $browserProcesses) {
        $running = @(Get-Process -Name $processName -ErrorAction SilentlyContinue).Count -gt 0
        [pscustomobject]@{
            Target = $processName
            Kind = "Process"
            Status = if ($running) { "Running" } else { "Stopped" }
            Details = if ($running) { "Process detected" } else { "No process detected" }
        }
    }

    $serviceStatus = foreach ($serviceConfig in $managedServices) {
        $record = Get-ServiceRecord -Name $serviceConfig.Name
        $task = Get-ResumeTask -TaskName $serviceConfig.TaskName
        $taskInfo = if ($task) { Get-ScheduledTaskInfo -TaskName $serviceConfig.TaskName } else { $null }

        if ($record) {
            [pscustomobject]@{
                Target = $serviceConfig.Name
                Kind = "Service"
                Status = $record.State
                Details = "Startup=$($record.StartMode); ResumeTask=" + $(if ($taskInfo) { $taskInfo.NextRunTime.ToString("yyyy-MM-dd HH:mm:ss") } else { "None" })
            }
        } else {
            [pscustomobject]@{
                Target = $serviceConfig.Name
                Kind = "Service"
                Status = "Missing"
                Details = "Service not installed"
            }
        }
    }

    $taskStatus = foreach ($serviceConfig in $managedServices) {
        $task = Get-ResumeTask -TaskName $serviceConfig.TaskName
        $taskInfo = if ($task) { Get-ScheduledTaskInfo -TaskName $serviceConfig.TaskName } else { $null }

        [pscustomobject]@{
            Target = $serviceConfig.TaskName
            Kind = "Task"
            Status = if ($task) { $task.State.ToString() } else { "Missing" }
            Details = if ($taskInfo) { "NextRun=$($taskInfo.NextRunTime.ToString("yyyy-MM-dd HH:mm:ss"))" } else { "No scheduled resume task" }
        }
    }

    $state = Load-State
    if ($state) {
        Write-Host "Saved state: enter at $($state.EnteredAt), resume at $($state.ResumeAt)" -ForegroundColor Cyan
    } else {
        Write-Host "Saved state: none" -ForegroundColor Yellow
    }

    Write-Host ""
    ($browserStatus + $serviceStatus + $taskStatus) | Format-Table -AutoSize
}

switch ($Mode) {
    "enter" { Enter-ExamMode }
    "restore" { Restore-ExamMode }
    "status" { Show-ExamModeStatus }
    default { throw "Unsupported mode: $Mode" }
}
