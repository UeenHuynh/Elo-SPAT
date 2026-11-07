% test_robustness_helpers  Focused unit tests for Tasks 28-29.

clear; clc;

rootDir = fileparts(mfilename('fullpath'));
addpath(fullfile(rootDir, '..', 'config'));
addpath(fullfile(rootDir, '..', 'evaluation'));
addpath(fullfile(rootDir, '..', 'robustness'));

fprintf('\n========================================\n');
fprintf('TEST: robustness helpers (Tasks 28-29)\n');
fprintf('========================================\n\n');

cfg = default_config();
cfg.robust.enabled = true;
cfg.robust.mode = 'COWC';
cfg.robust.setupMM = [3, 3, 3];
cfg.robust.rangePercent = 3.5;

pln = struct();
pln.radiationMode = "protons";
pln.machine = ['Generic'; 'Generic'];
pln = configure_robust_scenarios(pln, [], cfg);
assert(isfield(pln, 'multScen') && ~isempty(pln.multScen), ...
    'configure_robust_scenarios did not populate pln.multScen.');
assert(strcmp(pln.radiationMode, 'protons'), 'radiationMode normalization failed.');
assert(strcmp(pln.machine, 'Generic'), 'machine normalization failed.');
fprintf('[OK] Robust scenario metadata configured.\n');

cst = cell(2, 6);
cst{1, 2} = 'CTV';
cst{1, 3} = 'TARGET';
cst{1, 4} = {int32([1; 2])};
cst{1, 6} = {struct('parameters', {{5}}, 'robustness', 'none')};

cst{2, 2} = 'Rectum';
cst{2, 3} = 'OAR';
cst{2, 4} = {int32([3; 4])};
cst{2, 6} = {struct('parameters', {{2}}, 'robustness', 'none')};

cstRobust = configure_robust_objectives(cst, cfg);
assert(strcmp(cstRobust{1, 6}{1}.robustness, 'COWC'), ...
    'Objective robustness mode was not updated.');
fprintf('[OK] Robust objective mode configured.\n');

ct = struct('cubeDim', [4, 1, 1]);
dij = struct();
dij.physicalDose = {eye(4), diag([0.8, 0.9, 1.1, 1.2])};
dij.totalNumOfBixels = 4;
w = ones(4, 1);

robustEval = evaluate_robust_scenarios(ct, cstRobust, [], dij, pln, w, cfg);

assert(robustEval.nScenarios == 2, 'Expected two evaluated scenarios.');
assert(isfield(robustEval, 'worstCase') && isfinite(robustEval.worstCase.targetD95Min), ...
    'Worst-case metrics were not computed.');
assert(all(size(robustEval.stdDoseCube) == [4, 1, 1]), ...
    'Unexpected stdDoseCube dimensions.');
fprintf('[OK] Robust scenario evaluation completed.\n');

fprintf('\n========================================\n');
fprintf('TEST PASSED\n');
fprintf('========================================\n\n');
