% run_task7_30deg_save  Run Task 7 with 30° and save intermediate results
%
% This script runs Task 7 with 30° arc spacing and saves results
% immediately after each major step to ensure artifacts are persisted
% even if optimization doesn't fully converge.

fprintf('\n========================================\n');
fprintf('Task 7: 30° Arc with Intermediate Saves\n');
fprintf('========================================\n\n');

% Setup
addpath('D:\matRad');
matRad_rc;

configDir = fullfile(fileparts(mfilename('fullpath')), '..', 'config');
addpath(configDir);

% Find PROSTATE.mat
matRadPath = fileparts(which('matRad_rc'));
patientFile = fullfile(matRadPath, 'data', 'PROSTATE.mat');
load(patientFile, 'ct', 'cst');

% Configure
cfg = default_config();
cfg.arcSpacing = 30;

resultsDir = fullfile(fileparts(mfilename('fullpath')), '..', 'results');
if ~exist(resultsDir, 'dir')
    mkdir(resultsDir);
end

timestamp = datestr(now, 'yyyy-mm-dd_HHMMSS');

% Step 1: Build arc
fprintf('Step 1: Building 30° arc...\n');
tic;
[pln, stf, dij] = build_discrete_arc_plan(ct, cst, cfg);
t_arc = toc;

% Save arc immediately
arcFile = fullfile(resultsDir, sprintf('task7_30deg_arc_%s.mat', timestamp));
save(arcFile, 'pln', 'stf', 'dij', 'cfg', 't_arc', '-v7.3');
fprintf('Arc saved to: %s (%.1f MB)\n\n', arcFile, dir(arcFile).bytes/1e6);

% Step 2: Optimize
fprintf('Step 2: Running Dense Arc optimization...\n');
fprintf('NOTE: This may take 10-20 minutes. Results will be saved when complete.\n\n');

tic;
[resultGUI, meta] = run_dense_arc_baseline(dij, cst, pln);
t_opt = toc;

% Save final result
finalFile = fullfile(resultsDir, sprintf('task7_30deg_final_%s.mat', timestamp));
save(finalFile, 'pln', 'stf', 'dij', 'resultGUI', 'meta', 'cfg', ...
    't_arc', 't_opt', '-v7.3');

plotFiles = struct('dvh', '');
try
    plotFiles = export_plan_visualizations( ...
        ct, cst, dij, resultGUI, meta, resultsDir, ...
        sprintf('task7_30deg_final_%s', timestamp));
    fprintf('DVH plot saved to: %s\n', plotFiles.dvh);
catch plotME
    fprintf('[INFO] Could not save DVH plot: %s\n', plotME.message);
end

fprintf('\n========================================\n');
fprintf('COMPLETE\n');
fprintf('========================================\n\n');
fprintf('Arc file: %s (%.1f MB)\n', arcFile, dir(arcFile).bytes/1e6);
fprintf('Final file: %s (%.1f MB)\n', finalFile, dir(finalFile).bytes/1e6);
if isfield(plotFiles, 'dvh') && ~isempty(plotFiles.dvh)
    fprintf('DVH plot: %s\n', plotFiles.dvh);
end
fprintf('Converged: %s\n', mat2str(meta.converged));
fprintf('Exit flag: %d\n', meta.exitStatus);
fprintf('Total time: %.1f minutes\n', (t_arc + t_opt)/60);
