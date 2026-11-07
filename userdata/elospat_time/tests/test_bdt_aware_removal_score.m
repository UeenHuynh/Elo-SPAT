% test_bdt_aware_removal_score  Focused unit test for Task 25 scoring helper.

clear; clc;

planningDir = fullfile(fileparts(mfilename('fullpath')), '..', 'planning');
addpath(planningDir);

fprintf('\n========================================\n');
fprintf('TEST: bdt-aware removal score (Task 25)\n');
fprintf('========================================\n\n');

doseImportance = [10; 10; 5];
bdtSaving = [1; 5; 5];
scores = compute_bdt_aware_removal_score(doseImportance, bdtSaving, 1e-12);

assert(scores(2) < scores(1), ...
    'Higher BDT saving at equal dose importance should reduce the score.');
assert(scores(3) < scores(2), ...
    'Lower dose importance at equal BDT saving should reduce the score.');
assert(all(isfinite(scores)), 'Scores should remain finite for positive savings.');

fprintf('[OK] Removal score ranks low-dose/high-saving candidates first.\n');

fprintf('\n========================================\n');
fprintf('TEST PASSED\n');
fprintf('========================================\n\n');
