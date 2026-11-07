% test_task7_quick  Quick test for Task 7 with 30° arc spacing
%
% This script runs a faster end-to-end test of Task 7 using 30° arc spacing
% (12 control points) to verify the pipeline and artifact persistence.
%
% USAGE
%   Run from MATLAB command window:
%     cd('d:\Elo-SPAT\userdata\elospat_time\planning')
%     test_task7_quick

fprintf('\n========================================\n');
fprintf('QUICK TEST: Task 7 with 30° arc\n');
fprintf('========================================\n\n');

% -------------------------------------------------------------------------
% Setup
% -------------------------------------------------------------------------
fprintf('--- Setup ---\n');

if ~exist('matRad_rc', 'file')
    error('matRad not found on path. Run matRad_rc from matRad root directory.');
end

% Find PROSTATE.mat
possiblePaths = {
    'PROSTATE.mat',
    fullfile('..', '..', '..', 'matRad', 'data', 'PROSTATE.mat'),
    fullfile(fileparts(which('matRad_rc')), 'data', 'PROSTATE.mat')
};

patientFile = '';
for i = 1:numel(possiblePaths)
    if exist(possiblePaths{i}, 'file')
        patientFile = possiblePaths{i};
        break;
    end
end

if isempty(patientFile)
    error('PROSTATE.mat not found.');
end

fprintf('Loading patient data from: %s\n', patientFile);
load(patientFile, 'ct', 'cst');

% Load configuration
configDir = fullfile(fileparts(mfilename('fullpath')), '..', 'config');
addpath(configDir);

cfg = default_config();
cfg.arcSpacing = 30;  % 30-degree spacing → 12 control points (faster)

fprintf('[OK] Setup complete.\n');
fprintf('     Arc spacing: %d degrees (12 control points)\n\n', cfg.arcSpacing);

% -------------------------------------------------------------------------
% Step 1: Build discrete arc plan (Task 6)
% -------------------------------------------------------------------------
fprintf('========================================\n');
fprintf('STEP 1: Build discrete arc plan (30°)\n');
fprintf('========================================\n\n');

tic;
[pln, stf, dij] = build_discrete_arc_plan(ct, cst, cfg);
t_arc = toc;

fprintf('\n[OK] Arc plan built in %.1f seconds.\n', t_arc);
fprintf('     Control points: %d\n', numel(pln.propStf.gantryAngles));
fprintf('     Total bixels: %d\n\n', dij.totalNumOfBixels);

% -------------------------------------------------------------------------
% Step 2: Run Dense Arc baseline optimization (Task 7)
% -------------------------------------------------------------------------
fprintf('========================================\n');
fprintf('STEP 2: Dense Arc baseline (30°)\n');
fprintf('========================================\n\n');

tic;
[resultGUI, meta] = run_dense_arc_baseline(dij, cst, pln);
t_opt = toc;

fprintf('\n[OK] Optimization completed in %.1f seconds.\n', t_opt);
fprintf('     Converged: %s\n', mat2str(meta.converged));
fprintf('     Exit flag: %d\n', meta.exitStatus);
fprintf('     Status: %s\n\n', meta.statusMessage);

% -------------------------------------------------------------------------
% Step 3: Save artifacts
% -------------------------------------------------------------------------
fprintf('========================================\n');
fprintf('STEP 3: Save artifacts\n');
fprintf('========================================\n\n');

resultsDir = fullfile(fileparts(mfilename('fullpath')), '..', 'results');
if ~exist(resultsDir, 'dir')
    mkdir(resultsDir);
end

timestamp = datestr(now, 'yyyy-mm-dd_HHMMSS');
artifactFile = fullfile(resultsDir, sprintf('task7_quick_30deg_%s.mat', timestamp));

% Save all outputs
save(artifactFile, 'pln', 'stf', 'dij', 'resultGUI', 'meta', 'cfg', ...
    't_arc', 't_opt', '-v7.3');

fprintf('[OK] Artifacts saved to:\n');
fprintf('     %s\n', artifactFile);
fprintf('     File size: %.1f MB\n\n', dir(artifactFile).bytes / 1e6);

fprintf('Saving main visualization plot...\n');
try
    plotFiles = export_plan_visualizations( ...
        ct, cst, dij, resultGUI, meta, resultsDir, ...
        sprintf('task7_quick_30deg_%s', timestamp));
    fprintf('[OK] DVH plot saved: %s\n\n', plotFiles.dvh);
catch plotME
    fprintf('[INFO] Could not save DVH plot: %s\n\n', plotME.message);
end

% -------------------------------------------------------------------------
% Step 4: Verify outputs
% -------------------------------------------------------------------------
fprintf('========================================\n');
fprintf('STEP 4: Verify outputs\n');
fprintf('========================================\n\n');

% Check resultGUI
assert(isstruct(resultGUI), 'resultGUI must be a struct');
assert(isfield(resultGUI, 'w'), 'resultGUI.w missing');
assert(numel(resultGUI.w) == dij.totalNumOfBixels, 'Weight dimension mismatch');
fprintf('[OK] resultGUI validated.\n');

% Check meta
assert(isstruct(meta), 'meta must be a struct');
assert(strcmp(meta.methodName, 'Dense Arc baseline'), 'Wrong method name');
assert(meta.nControlPoints == 12, 'Expected 12 control points');
assert(meta.arcSpacing == 30, 'Expected 30° arc spacing');
fprintf('[OK] meta validated.\n');

% Check convergence
if meta.converged
    fprintf('[OK] Optimization converged successfully.\n');
else
    fprintf('[WARNING] Optimization did not converge (exit flag: %d).\n', meta.exitStatus);
    if meta.exitStatus == -1
        fprintf('          Max iterations reached - this is acceptable for testing.\n');
    end
end

% -------------------------------------------------------------------------
% Summary
% -------------------------------------------------------------------------
fprintf('\n========================================\n');
fprintf('QUICK TEST SUMMARY\n');
fprintf('========================================\n\n');

fprintf('Configuration:\n');
fprintf('  Arc spacing: %d degrees\n', cfg.arcSpacing);
fprintf('  Control points: %d\n', meta.nControlPoints);
fprintf('  Total bixels: %d\n', meta.nBixels);
fprintf('\n');

fprintf('Performance:\n');
fprintf('  Arc build time: %.1f seconds\n', t_arc);
fprintf('  Optimization time: %.1f seconds\n', t_opt);
fprintf('  Total time: %.1f seconds\n', t_arc + t_opt);
fprintf('\n');

fprintf('Results:\n');
fprintf('  Method: %s\n', meta.methodName);
fprintf('  Converged: %s\n', mat2str(meta.converged));
fprintf('  Exit flag: %d\n', meta.exitStatus);
fprintf('  Status: %s\n', meta.statusMessage);
fprintf('  Optimized weights: %d\n', numel(resultGUI.w));
fprintf('\n');

fprintf('Artifacts:\n');
fprintf('  Saved to: %s\n', artifactFile);
fprintf('  Size: %.1f MB\n', dir(artifactFile).bytes / 1e6);
if exist('plotFiles', 'var') && isstruct(plotFiles) && isfield(plotFiles, 'dvh') && ~isempty(plotFiles.dvh)
    fprintf('  DVH plot: %s\n', plotFiles.dvh);
end
fprintf('\n');

fprintf('Status: PASSED\n');
fprintf('========================================\n\n');

fprintf('Task 7 quick test (30°) completed successfully.\n');
fprintf('Ready to run full test with 20° arc spacing.\n\n');
