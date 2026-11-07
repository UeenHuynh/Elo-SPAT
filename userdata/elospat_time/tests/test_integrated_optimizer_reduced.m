% test_integrated_optimizer_reduced  Small sparse end-to-end engine test.

rootDir = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(fullfile(rootDir, 'integrated'));
addpath(fullfile(rootDir, 'config'));

n = 12;
A = sparse(reshape(sin(1:(20*n)), 20, n));
A(abs(A) < 0.45) = 0;
xTrue = zeros(n, 1); xTrue([2, 3, 8]) = [1.0, 0.5, 0.8];
b = A * xTrue;
model = struct();
model.nWeights = n;
model.hasConstraints = false;
model.objectiveFunction = @(x) 0.5 * sum((A*x - b).^2);
model.objectiveGradient = @(x) A' * (A*x - b);

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
cfg.lambdaEL = 2e-2;
cfg.lambdaS = 1e-3;
cfg.lambdaES = 0;
cfg.lambdaT = 0;
cfg.integrated.maxIterations = 150;
cfg.integrated.initialStep = 1e-2;
cfg.integrated.maximumStep = 0.1;
cfg.integrated.verbose = false;
context = build_integrated_context( ...
    groups, cpGroups, map, default_machine_model(), cfg);
w0 = ones(n, 1) * 0.2;
initialObjective = evaluate_integrated_objective(w0, model, context, cfg);
result = run_integrated_optimizer(model, context, w0, cfg);
assert(result.iterations >= cfg.integrated.minimumIterations || ...
        ~strcmp(result.status, 'converged'), ...
    'Stationarity convergence was declared before minimumIterations.');
assert(result.objective <= initialObjective + 1e-8, ...
    'Integrated optimizer did not reduce its objective.');
assert(all(diff(result.history.objective) <= 1e-8), ...
    'Integrated objective must be monotone under MM/backtracking.');
assert(all(result.w >= 0) && all(isfinite(result.w)), ...
    'Integrated optimizer returned invalid weights.');
assert(isfield(result.history, 'normalizedStationarity') && ...
        all(isfinite(result.history.normalizedStationarity)) && ...
        all(result.history.normalizedStationarity >= 0), ...
    'Stationarity diagnostics were not recorded correctly.');
assert(isfield(result, 'stoppingReason') && ~isempty(result.stoppingReason), ...
    'Optimizer did not report an explicit stopping reason.');
support = extract_integrated_support(result.w, w0, context, cfg);
assert(all(result.w(~support.mask) < max(cfg.spotThreshold, cfg.mmuThreshold)), ...
    'Support extraction retained inconsistent excluded weights.');
for bIx = 1:numel(cpGroups)
    assert(any(support.mask(cpGroups(bIx).ix)), ...
        'Control-point support coverage was not preserved.');
end
fprintf('[OK] Reduced integrated optimizer and support extraction validated.\n');
