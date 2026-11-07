% test_reoptimize_reduced_problem  Validation for Task 24 dimension safety.

clear; clc;

reductionDir = fullfile(fileparts(mfilename('fullpath')), '..', 'reduction');
addpath(reductionDir);

fprintf('\n========================================\n');
fprintf('TEST: reoptimize_reduced_problem (Task 24)\n');
fprintf('========================================\n\n');

dijReduced = struct('totalNumOfBixels', 3);
cst = cell(1, 6);
pln = struct('radiationMode', 'protons');
wInitWrong = [1; 2];

caughtExpectedError = false;

try
    reoptimize_reduced_problem(dijReduced, cst, pln, wInitWrong);
catch ME
    caughtExpectedError = strcmp(ME.identifier, ...
        'reoptimize_reduced_problem:dimensionMismatch');
    if ~caughtExpectedError
        rethrow(ME);
    end
end

assert(caughtExpectedError, ...
    'Expected a dimensionMismatch error for invalid warm-start length.');

fprintf('[OK] Warm-start dimension check validated.\n');

fprintf('\n========================================\n');
fprintf('TEST PASSED\n');
fprintf('========================================\n\n');
