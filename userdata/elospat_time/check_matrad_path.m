% check_matrad_path.m
%
% RMS-ELO-SPAT — matRad Path Verification Utility
%
% PURPOSE:
%   Checks whether matRad is available on the current MATLAB path and
%   whether the required patient data file (PROSTATE.mat) can be found.
%   This script only WARNS if matRad is missing — it does NOT throw an
%   error so that infrastructure tasks can be validated independently.
%
% USAGE:
%   Run from the userdata/elospat_time/ directory:
%       check_matrad_path
%
% OUTPUTS:
%   Prints status messages to the command window.
%   Returns logical variables matradFound and prostateFound (if called
%   as a function-style script with output capture, use evalc or similar).
%
% =========================================================================
% RESEARCH-ONLY DISCLAIMER:
%   This script is part of RMS-ELO-SPAT, a research-only module.
%   Not for clinical use. matRad is used as a research sandbox only.
% =========================================================================

fprintf('\n=== RMS-ELO-SPAT: matRad Path Check ===\n\n');

% -------------------------------------------------------------------------
% CHECK 1: matRad runtime context (matRad_rc)
% -------------------------------------------------------------------------
% matRad_rc is the canonical indicator that matRad is on the MATLAB path.
% It is a script/function at the root of the matRad installation that sets
% up the matRad runtime environment and must be callable for matRad to work.

matradFound = (exist('matRad_rc', 'file') == 2);  % 2 = m-file on path

if matradFound
    fprintf('[OK]  matRad_rc found on MATLAB path.\n');
    fprintf('      matRad appears to be properly added to the path.\n');

    % Try to locate the matRad root directory for informational output
    matradRcPath = which('matRad_rc');
    if ~isempty(matradRcPath)
        [matradRoot, ~, ~] = fileparts(matradRcPath);
        fprintf('      matRad root: %s\n', matradRoot);
    end
else
    fprintf('[WARN] matRad_rc NOT found on MATLAB path.\n');
    fprintf('       matRad does not appear to be available.\n\n');
    fprintf('  --> To fix this, clone matRad from GitHub and add it to path:\n');
    fprintf('       1. Clone:  git clone https://github.com/e0404/matRad.git\n');
    fprintf('       2. In MATLAB, navigate to the matRad root directory.\n');
    fprintf('       3. Run:    matRad_rc\n');
    fprintf('          OR:     addpath(genpath(''path/to/matRad''))\n\n');
    fprintf('  NOTE: Tasks 1-4 (infrastructure) can proceed without matRad.\n');
    fprintf('        Tasks 5+ require matRad to be on the path.\n');
end

fprintf('\n');

% -------------------------------------------------------------------------
% CHECK 2: PROSTATE.mat patient data file
% -------------------------------------------------------------------------
% PROSTATE.mat is distributed with matRad in its data/ subdirectory.
% We search for it using MATLAB's built-in file locator, then also try
% common relative paths relative to matRad root (if found above).

prostateFound = false;
prostatePath  = '';

% Try MATLAB path search first
prostateOnPath = which('PROSTATE.mat');
if ~isempty(prostateOnPath)
    prostateFound = true;
    prostatePath  = prostateOnPath;
end

% If not found via which(), try searching relative to matRad root
if ~prostateFound && matradFound
    matradRcPath = which('matRad_rc');
    if ~isempty(matradRcPath)
        [matradRoot, ~, ~] = fileparts(matradRcPath);

        % Common locations within a matRad installation
        candidatePaths = {
            fullfile(matradRoot, 'data', 'PROSTATE.mat')
            fullfile(matradRoot, 'PROSTATE.mat')
            fullfile(matradRoot, 'examples', 'PROSTATE.mat')
        };

        for k = 1:numel(candidatePaths)
            if exist(candidatePaths{k}, 'file') == 2
                prostateFound = true;
                prostatePath  = candidatePaths{k};
                break;
            end
        end
    end
end

if prostateFound
    fprintf('[OK]  PROSTATE.mat found.\n');
    fprintf('      Path: %s\n', prostatePath);
else
    fprintf('[WARN] PROSTATE.mat NOT found.\n');
    if matradFound
        fprintf('       matRad is on the path but PROSTATE.mat was not located.\n');
        fprintf('  --> PROSTATE.mat should be distributed with matRad in its data/\n');
        fprintf('      subdirectory. Verify your matRad installation is complete.\n');
        fprintf('      Expected location: <matRad_root>/data/PROSTATE.mat\n');
    else
        fprintf('       (Cannot search matRad data paths since matRad is not on path.)\n');
        fprintf('  --> Once matRad is installed, PROSTATE.mat should be at:\n');
        fprintf('      <matRad_root>/data/PROSTATE.mat\n');
    end
end

fprintf('\n');

% -------------------------------------------------------------------------
% SUMMARY
% -------------------------------------------------------------------------
fprintf('--- Summary ---\n');
if matradFound
    matradStr = 'YES';
else
    matradStr = 'NO  (see warnings above)';
end
if prostateFound
    prostateStr = 'YES';
else
    prostateStr = 'NO  (see warnings above)';
end
fprintf('  matRad on path:    %s\n', matradStr);
fprintf('  PROSTATE.mat:      %s\n', prostateStr);

if matradFound && prostateFound
    fprintf('\n[READY] All requirements met. You can proceed with run_rms_elospat_main.m\n');
elseif matradFound && ~prostateFound
    fprintf('\n[PARTIAL] matRad found, but PROSTATE.mat is missing.\n');
    fprintf('          Tasks 1-7 (except full integration tests) may still run.\n');
    fprintf('          Locate PROSTATE.mat before running end-to-end tests.\n');
else
    fprintf('\n[NOT READY] Please install matRad before proceeding beyond Task 4.\n');
    fprintf('            Infrastructure tasks (Tasks 1-4) can still be completed.\n');
end

fprintf('==========================================\n\n');
