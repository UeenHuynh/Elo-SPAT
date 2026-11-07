% test_run_post_filter_baseline_smoke  Synthetic smoke test for Task 36.

clear; clc;

rootDir = fileparts(mfilename('fullpath'));
addpath(fullfile(rootDir, 'mocks'));
addpath(fullfile(rootDir, '..', 'config'));
addpath(fullfile(rootDir, '..', 'delivery'));
addpath(fullfile(rootDir, '..', 'evaluation'));
addpath(fullfile(rootDir, '..', 'mapping'));
addpath(fullfile(rootDir, '..', 'planning'));
addpath(fullfile(rootDir, '..', 'reduction'));

fprintf('\n========================================\n');
fprintf('TEST: run_post_filter_baseline smoke (Task 36)\n');
fprintf('========================================\n');

cfg = default_config();
cfg.spotThreshold = 0.05;
cfg.mmuThreshold = 0.10;
machine = default_machine_model();

map = struct( ...
    'globalIx', num2cell((1:6)'), ...
    'beamNum', num2cell([1; 1; 1; 2; 2; 2]), ...
    'rayNum', num2cell([1; 1; 1; 1; 1; 1]), ...
    'bixelNum', num2cell([1; 2; 3; 1; 2; 3]), ...
    'gantryAngle', num2cell([0; 0; 0; 20; 20; 20]), ...
    'energy', num2cell([150; 150; 140; 150; 140; 140]), ...
    'cpIndex', num2cell([1; 1; 1; 2; 2; 2]), ...
    'energyIndex', num2cell([1; 1; 2; 1; 2; 2]));

[groups, ~] = build_group_index(map);

dij = struct();
dij.totalNumOfBixels = 6;
dij.physicalDose = {diag([3, 2.5, 2, 3, 2.2, 1.8])};
dij.doseGrid = struct('dimensions', [6, 1, 1]);
dij.beamNum = (1:6)';
dij.rayNum = ones(6, 1);
dij.bixelNum = (1:6)';

cst = cell(1, 6);
cst{1, 2} = 'CTV';
cst{1, 3} = 'TARGET';
cst{1, 4} = {int32([1; 2; 3; 4; 5; 6])};
cst{1, 6} = {struct('parameters', {{5}}, 'robustness', 'none')};

pln = struct();
wInit = [1.2; 0.08; 0.2; 1.0; 0.15; 0.05];

[result, spotMask] = run_post_filter_baseline(dij, cst, pln, wInit, map, groups, cfg, machine);

assert(isfield(result, 'w') && numel(result.w) == 6);
assert(isfield(result, 'spotMask') && numel(result.spotMask) == 6);
assert(sum(spotMask) == 4);
assert(all(result.w(~spotMask) == 0));
assert(isfield(result, 'bdtMetrics') && isstruct(result.bdtMetrics));
assert(abs(result.thresholdUsed - 0.10) < 1e-12);

fprintf('[PASS] Post-filter baseline completed with expected pruning.\n');
