% test_l21half_gradient  Finite-difference validation for Task 11.

clear; clc;

objectiveDir = fullfile(fileparts(mfilename('fullpath')), '..', 'objective');
addpath(objectiveDir);

fprintf('\n========================================\n');
fprintf('TEST: l21half_group_penalty gradient (Task 11)\n');
fprintf('========================================\n\n');

rng(42);

w = [0.8; 0.1; 0.6; 0.3; 1.2; 0.5; 0.9];
groups = repmat(struct('ix', []), 3, 1);
groups(1).ix = [1; 2];
groups(2).ix = [3; 4; 5];
groups(3).ix = [6; 7];
epsEL = 1e-6;

[val, grad] = l21half_group_penalty(w, groups, epsEL);

assert(isscalar(val) && isfinite(val), 'Penalty value must be a finite scalar.');
assert(isequal(size(grad), size(w)), 'Gradient size must match w.');
fprintf('[OK] Penalty and gradient shapes validated.\n');

h = 1e-6;
n = numel(w);
relErr = zeros(n, 1);
fdGrad = zeros(n, 1);

for i = 1:n
    e = zeros(n, 1);
    e(i) = 1;

    valPlus = l21half_group_penalty(w + h * e, groups, epsEL);
    valMinus = l21half_group_penalty(w - h * e, groups, epsEL);
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
