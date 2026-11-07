% test_integrated_optimizer_resume  Exact checkpoint continuation equivalence.

rootDir = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(fullfile(rootDir, 'integrated'));
addpath(fullfile(rootDir, 'config'));
addpath(fullfile(rootDir, 'mapping'));
addpath(fullfile(rootDir, 'objective'));

n = 12;
A = sparse(reshape(sin(1:(20*n)), 20, n));
A(abs(A) < 0.45) = 0;
xTrue = zeros(n, 1); xTrue([2, 3, 8]) = [1.0, 0.5, 0.8];
b = A * xTrue;
model = struct('nWeights', n, 'hasConstraints', false, ...
    'objectiveFunction', @(x) 0.5 * sum((A*x - b).^2), ...
    'objectiveGradient', @(x) A' * (A*x - b));

groups = repmat(struct('ix', [], 'cpIndex', 0, 'energy', 0), 4, 1);
cpGroups = repmat(struct('ix', [], 'cpIndex', 0, 'gantryAngle', 0), 2, 1);
for g = 1:4
    groups(g).ix = ((g-1)*3+1:g*3)';
    groups(g).cpIndex = 1 + (g > 2);
    groups(g).energy = 90 + 10*g;
end
cpGroups(1) = struct('ix', (1:6)', 'cpIndex', 1, 'gantryAngle', 0);
cpGroups(2) = struct('ix', (7:12)', 'cpIndex', 2, 'gantryAngle', 20);
map = repmat(struct('globalIx', 0), n, 1);
for i = 1:n, map(i).globalIx = i; end

cfg = default_config();
cfg.lambdaEL = 2e-2; cfg.lambdaS = 0; cfg.lambdaES = 0; cfg.lambdaT = 0;
cfg.integrated.initialStep = 1e-2;
cfg.integrated.maximumStep = 0.1;
cfg.integrated.minimumIterations = 1000;
cfg.integrated.verbose = false;
context = build_integrated_context( ...
    groups, cpGroups, map, default_machine_model(), cfg);
w0 = ones(n, 1) * 0.2;

cfgFull = cfg; cfgFull.integrated.maxIterations = 13;
fullRun = run_integrated_optimizer(model, context, w0, cfgFull);

cfgFirst = cfg; cfgFirst.integrated.maxIterations = 6;
firstRun = run_integrated_optimizer(model, context, w0, cfgFirst);
cfgSecond = cfg; cfgSecond.integrated.maxIterations = 7;
resumed = run_integrated_optimizer( ...
    model, context, firstRun.w, cfgSecond, firstRun);

assert(resumed.resumedFromIteration == 6 && resumed.segmentIterations == 7);
assert(resumed.iterations == fullRun.iterations);
assert(isequal(resumed.history.iteration, fullRun.history.iteration));
assert(max(abs(resumed.w - fullRun.w)) <= 1e-14, ...
    'Resumed weights differ from uninterrupted execution.');
assert(abs(resumed.objective - fullRun.objective) <= 1e-14, ...
    'Resumed objective differs from uninterrupted execution.');
assert(max(abs(resumed.history.objective - fullRun.history.objective)) <= 1e-14);
assert(max(abs(resumed.history.stepSize - fullRun.history.stepSize)) <= 1e-14);
assert(max(abs(resumed.history.normalizedStationarity - ...
    fullRun.history.normalizedStationarity)) <= 1e-14);
assert(resumed.resumeObjectiveRelativeError <= 1e-14);
fprintf('[OK] Exact integrated-optimizer checkpoint continuation validated.\n');

% Repeat across the middle of a frozen spot-MM block. This verifies that
% q and its anchor iteration, not just primal weights, survive a checkpoint.
cfgMM = cfg;
cfgMM.lambdaS = 3e-3;
cfgMM.epsS = 1e-3;
cfgMM.integrated.spotMMUpdatePeriod = 5;
cfgMM.integrated.maxIterations = 13;
fullMM = run_integrated_optimizer(model, context, w0, cfgMM);
cfgMMFirst = cfgMM; cfgMMFirst.integrated.maxIterations = 7;
firstMM = run_integrated_optimizer(model, context, w0, cfgMMFirst);
cfgMMSecond = cfgMM; cfgMMSecond.integrated.maxIterations = 6;
resumedMM = run_integrated_optimizer( ...
    model, context, firstMM.w, cfgMMSecond, firstMM);
assert(max(abs(resumedMM.w-fullMM.w)) <= 1e-14);
assert(abs(resumedMM.objective-fullMM.objective) <= 1e-14);
assert(isequal(resumedMM.history.spotMMUpdated,fullMM.history.spotMMUpdated));
assert(isequal(resumedMM.history.spotMMAnchorIteration, ...
    fullMM.history.spotMMAnchorIteration));
assert(max(abs(resumedMM.spotMMWeights-fullMM.spotMMWeights)) <= 1e-14);
assert(resumedMM.spotMMAnchorIteration == fullMM.spotMMAnchorIteration);
fprintf('[OK] Frozen spot-MM checkpoint continuation validated.\n');
