% test_barrier_gradient  Finite-difference validation for Task 12.

clear; clc;

objectiveDir = fullfile(fileparts(mfilename('fullpath')), '..', 'objective');
mappingDir = fullfile(fileparts(mfilename('fullpath')), '..', 'mapping');
addpath(objectiveDir);
addpath(mappingDir);

fprintf('\n========================================\n');
fprintf('TEST: control_point_barrier gradient (Task 12)\n');
fprintf('========================================\n\n');

rng(43);

w = [0.4; 0.2; 0.9; 0.3; 1.1; 0.5; 0.7];
cpGroups = repmat(struct('ix', []), 3, 1);
cpGroups(1).ix = [1; 2];
cpGroups(2).ix = [3; 4; 5];
cpGroups(3).ix = [6; 7];
epsB = 1e-6;

[val, grad] = control_point_barrier(w, cpGroups, epsB);

assert(isscalar(val) && isfinite(val), 'Barrier value must be a finite scalar.');
assert(isequal(size(grad), size(w)), 'Gradient size must match w.');
fprintf('[OK] Barrier and gradient shapes validated.\n');

h = 1e-6;
n = numel(w);
relErr = zeros(n, 1);
fdGrad = zeros(n, 1);

for i = 1:n
    e = zeros(n, 1);
    e(i) = 1;

    valPlus = control_point_barrier(w + h * e, cpGroups, epsB);
    valMinus = control_point_barrier(w - h * e, cpGroups, epsB);
    fdGrad(i) = (valPlus - valMinus) / (2 * h);

    denom = max(1, abs(fdGrad(i)), abs(grad(i)));
    relErr(i) = abs(fdGrad(i) - grad(i)) / denom;
end

passMask = relErr < 1e-4;
passRate = mean(passMask);

fprintf('[OK] Finite-difference check completed.\n');
fprintf('     Pass rate: %.2f%%\n', 100 * passRate);
fprintf('     Max relative error: %.3e\n', max(relErr));

assert(passRate >= 0.95, ...
    'Finite-difference pass rate %.2f%% is below required 95%%.', 100 * passRate);

fprintf('\n========================================\n');
fprintf('TEST PASSED\n');
fprintf('========================================\n\n');
