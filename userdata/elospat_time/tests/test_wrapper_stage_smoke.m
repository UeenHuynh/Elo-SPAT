% test_wrapper_stage_smoke  Synthetic smoke tests for Tasks 25-27.

clear; clc;

rootDir = fileparts(mfilename('fullpath'));
addpath(fullfile(rootDir, 'mocks'));
addpath(fullfile(rootDir, '..', 'config'));
addpath(fullfile(rootDir, '..', 'delivery'));
addpath(fullfile(rootDir, '..', 'evaluation'));
addpath(fullfile(rootDir, '..', 'mapping'));
addpath(fullfile(rootDir, '..', 'objective'));
addpath(fullfile(rootDir, '..', 'planning'));
addpath(fullfile(rootDir, '..', 'reduction'));
addpath(fullfile(rootDir, '..', 'robustness'));

fprintf('\n========================================\n');
fprintf('TEST: wrapper stage smoke tests (Tasks 25-27)\n');
fprintf('========================================\n\n');

cfg = default_config();
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

[groups, cpGroups] = build_group_index(map);

dij = struct();
dij.totalNumOfBixels = 6;
dij.physicalDose = {diag([3, 2.5, 2, 3, 2.2, 1.8])};
dij.doseGrid = struct('dimensions', [6, 1, 1]);
dij.beamNum = (1:6)';
dij.rayNum = ones(6, 1);
dij.bixelNum = (1:6)';

cst = cell(2, 6);
cst{1, 2} = 'CTV';
cst{1, 3} = 'TARGET';
cst{1, 4} = {int32([1; 2; 3])};
cst{1, 6} = {struct('parameters', {{5}}, 'robustness', 'none')};
cst{2, 2} = 'Rectum';
cst{2, 3} = 'OAR';
cst{2, 4} = {int32([4; 5; 6])};
cst{2, 6} = {struct('parameters', {{2}}, 'robustness', 'none')};

pln = struct();
wInit = [1.2; 0.8; 0.2; 1.0; 0.15; 0.05];

[resultLayer, layerMask] = run_elospat_like_stage(dij, cst, pln, wInit, groups, cpGroups, cfg, machine);
assert(isfield(resultLayer, 'w') && numel(resultLayer.w) == 6, ...
    'Layer stage should return embedded full-length weights.');
assert(numel(layerMask) == numel(groups), 'Layer mask length mismatch.');
fprintf('[OK] Task 25 smoke test passed.\n');

[resultSpot, spotMask] = run_spot_sparse_stage(dij, cst, pln, resultLayer.w, map, groups, cfg, machine);
assert(isfield(resultSpot, 'w') && numel(resultSpot.w) == 6, ...
    'Spot stage should return full-length weights.');
assert(numel(spotMask) == 6, 'Spot mask length mismatch.');
[spotReplay, ~] = run_spot_sparse_stage(dij, cst, pln, resultLayer.w, ...
    map, groups, cfg, machine, true);
assert(numel(spotReplay.supportSnapshots) == numel(spotReplay.history), ...
    'Snapshot capture should preserve every accepted history support.');
for k = 1:numel(spotReplay.supportSnapshots)
    snap = spotReplay.supportSnapshots(k);
    assert(nnz(snap.supportMask) == snap.nActiveSpots, ...
        'Snapshot spot count must match its support mask.');
    assert(all(snap.w(~snap.supportMask) == 0), ...
        'Snapshot weights must remain zero outside support.');
end
fprintf('[OK] Task 26 smoke test passed.\n');

cfg.robust.enabled = true;
supportMask = logical([1; 0; 1; 0; 1; 0]);
supportW = resultSpot.w;
supportW(~supportMask) = 0;
finalResult = run_final_reoptimization(dij, cst, pln, supportW, cfg, supportMask);
assert(isfield(finalResult, 'resultGUI') && isfield(finalResult.resultGUI, 'w'), ...
    'Final reoptimization should return a resultGUI with weights.');
assert(all(abs(finalResult.w(~supportMask)) < 1e-10), ...
    'Final reoptimization must preserve the supplied support.');
assert(strcmp(finalResult.robustMode, 'COWC'), 'Expected robust final mode.');
assert(isfield(finalResult, 'robustDeferred') && finalResult.robustDeferred, ...
    'Robust mode should be deferred for nominal-only final reoptimization.');
fprintf('[OK] Task 27 smoke test passed.\n');

fprintf('\n========================================\n');
fprintf('TEST PASSED\n');
fprintf('========================================\n\n');
