% test_integrated_prox  Validate group-half proximal thresholding.

rootDir = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(fullfile(rootDir, 'integrated'));

groups = struct('ix', {[1; 2]});
step = 0.4;
lambda = 0.7;
spotWeights = [0.1; 0.2];
z = [2.0; 0.8];
w = prox_rms_elospat(z, step, groups, lambda, spotWeights);
assert(all(w >= 0), 'The proximal result must be non-negative.');

% For fixed direction after one-sided soft thresholding, compare the radial
% closed form against a dense one-dimensional objective scan.
v = max(0, z - step * spotWeights);
radius = norm(v);
grid = linspace(0, radius + 1, 200001)';
objective = 0.5 * (grid - radius) .^ 2 + step * lambda * sqrt(grid);
[~, ix] = min(objective);
expectedRadius = grid(ix);
assert(abs(norm(w) - expectedRadius) < 2e-4, ...
    'Group-half proximal radius disagrees with direct minimization.');

tau = step * lambda;
small = [1; 0] * (1.5 * tau ^ (2 / 3));
wZero = prox_rms_elospat(small, step, groups, lambda, zeros(2, 1));
assert(nnz(wZero) == 0, 'Threshold equality must select the zero minimizer.');
fprintf('[OK] Integrated group-half proximal map validated.\n');
