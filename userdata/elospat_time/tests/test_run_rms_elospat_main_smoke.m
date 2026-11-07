% test_run_rms_elospat_main_smoke  Precomputed smoke test for Task 37.

clear; clc;

rootDir = fileparts(mfilename('fullpath'));
workspaceDir = fileparts(rootDir);
addpath(fullfile(rootDir, 'mocks'));
addpath(fullfile(workspaceDir, 'config'));
addpath(fullfile(workspaceDir, 'delivery'));
addpath(fullfile(workspaceDir, 'evaluation'));
addpath(fullfile(workspaceDir, 'mapping'));
addpath(fullfile(workspaceDir, 'objective'));
addpath(fullfile(workspaceDir, 'planning'));
addpath(fullfile(workspaceDir, 'reduction'));
addpath(workspaceDir);

fprintf('\n========================================\n');
fprintf('TEST: run_rms_elospat_main smoke (Task 37)\n');
fprintf('========================================\n');

cfg = default_config();
cfg.outputDir = fullfile(tempdir, 'rms_elospat_main_smoke');
cfg.captureConsoleOutput = false;
cfg.robust.enabled = false;

ct = struct('cubeDim', [6 1 1], 'resolution', struct('x', 1, 'y', 1, 'z', 1));
cst = cell(2, 6);
cst{1, 2} = 'CTV';
cst{1, 3} = 'TARGET';
cst{1, 4} = {int32([1; 2; 3])};
cst{1, 6} = {struct('parameters', {{5}}, 'robustness', 'none')};
cst{2, 2} = 'Rectum';
cst{2, 3} = 'OAR';
cst{2, 4} = {int32([4; 5; 6])};
cst{2, 6} = {struct('parameters', {{2}}, 'robustness', 'none')};

map = struct( ...
    'globalIx', num2cell((1:6)'), ...
    'beamNum', num2cell([1; 1; 1; 2; 2; 2]), ...
    'rayNum', num2cell([1; 1; 1; 1; 1; 1]), ...
    'bixelNum', num2cell([1; 2; 3; 1; 2; 3]), ...
    'gantryAngle', num2cell([0; 0; 0; 20; 20; 20]), ...
    'energy', num2cell([150; 150; 140; 150; 140; 140]), ...
    'cpIndex', num2cell([1; 1; 1; 2; 2; 2]), ...
    'energyIndex', num2cell([1; 1; 2; 1; 2; 2]));
[groups, cpGroups] = build_group_index(map);

dij = struct();
dij.totalNumOfBixels = 6;
dij.physicalDose = {diag([3, 2.5, 2, 3, 2.2, 1.8])};
dij.doseGrid = struct('dimensions', [6, 1, 1], 'numOfVoxels', 6);
dij.beamNum = (1:6)';
dij.rayNum = ones(6, 1);
dij.bixelNum = (1:6)';

pln = struct();
stf = struct();
wDense = [1.2; 0.8; 0.2; 1.0; 0.15; 0.05];
wImpt = [1.0; 0.6; 0.1; 0.8; 0.1; 0.02];

imptResult = struct('w', wImpt, 'info', struct('objective', 120, 'iter', 15, 'status', 0, 'cpu', 1.2), ...
    'usedOptimizer', 'IPOPT', 'converged', true, 'statusMessage', 'Solved');
imptMeta = struct('runtime', 1.3, 'converged', true, 'exitStatus', 0, 'optimizerName', 'IPOPT');

denseResult = struct('w', wDense, 'info', struct('objective', 100, 'iter', 20, 'status', 0, 'cpu', 1.5), ...
    'usedOptimizer', 'IPOPT', 'converged', true, 'statusMessage', 'Solved');
denseMeta = struct('runtime', 1.6, 'converged', true, 'exitStatus', 0, 'optimizerName', 'IPOPT');

cfg.precomputed = struct();
cfg.precomputed.ct = ct;
cfg.precomputed.cst = cst;
cfg.precomputed.impt = struct('pln', pln, 'stf', stf, 'dij', dij, 'resultGUI', imptResult, 'meta', imptMeta, ...
    'map', map, 'groups', groups, 'cpGroups', cpGroups);
cfg.precomputed.arc = struct('pln', pln, 'stf', stf, 'dij', dij, 'denseResult', denseResult, 'denseMeta', denseMeta, ...
    'map', map, 'groups', groups, 'cpGroups', cpGroups);

outputs = run_rms_elospat_main(cfg);

assert(istable(outputs.methodTable));
assert(height(outputs.methodTable) == 5);
assert(exist(outputs.files.matFile, 'file') == 2);
assert(exist(outputs.files.methodCsv, 'file') == 2);
assert(exist(outputs.files.latexFile, 'file') == 2);
assert(exist(outputs.files.summaryPlot, 'file') == 2);
assert(strcmp(outputs.robust.status, 'not_requested'));

fprintf('[PASS] Orchestrator completed and exported summary artifacts.\n');
