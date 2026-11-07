% test_compute_dvh_metrics  Focused unit test for Task 30 DVH metrics.

clear; clc;

evaluationDir = fullfile(fileparts(mfilename('fullpath')), '..', 'evaluation');
addpath(evaluationDir);

fprintf('\n========================================\n');
fprintf('TEST: compute_dvh_metrics (Task 30)\n');
fprintf('========================================\n\n');

doseCube = reshape(1:8, [2, 2, 2]);
ct = struct('cubeDim', [2, 2, 2]);

cst = cell(2, 6);
cst{1, 2} = 'CTV';
cst{1, 3} = 'TARGET';
cst{1, 4} = {int32([1; 2; 3; 4])};
cst{1, 6} = {struct('parameters', {{5}})};

cst{2, 2} = 'Rectum';
cst{2, 3} = 'OAR';
cst{2, 4} = {int32([5; 6; 7; 8])};
cst{2, 6} = {struct('parameters', {{3}})};

cfg = struct();
metrics = compute_dvh_metrics(doseCube, ct, cst, cfg);

assert(numel(metrics.structures) == 2, 'Expected two structure metrics.');
assert(abs(metrics.structures(1).D95 - 1.15) < 1e-9, 'Unexpected CTV D95.');
assert(abs(metrics.structures(1).D98 - 1.06) < 1e-9, 'Unexpected CTV D98.');
assert(abs(metrics.structures(1).D5 - 3.85) < 1e-9, 'Unexpected CTV D5.');
assert(abs(metrics.structures(1).Dmean - 2.5) < 1e-9, 'Unexpected CTV Dmean.');
assert(abs(metrics.structures(1).Dmax - 4.0) < 1e-9, 'Unexpected CTV Dmax.');
assert(abs(metrics.structures(1).V95 - 0.0) < 1e-9, 'Unexpected CTV V95.');
assert(abs(metrics.summary.targetD95Min - metrics.structures(1).D95) < 1e-9, ...
    'Target summary did not propagate D95 correctly.');
assert(abs(metrics.summary.oarDmaxMax - 8.0) < 1e-9, 'Unexpected OAR Dmax summary.');

fprintf('[OK] DVH metrics computed as expected.\n');

fprintf('\n========================================\n');
fprintf('TEST PASSED\n');
fprintf('========================================\n\n');
