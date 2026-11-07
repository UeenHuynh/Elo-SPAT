function result = run_stage2_admm_variant_h(problem, settings)
%RUN_STAGE2_ADMM_VARIANT_H Reduced MATLAB-native hierarchical ADMM.
%
% The represented Stage-1 layers are allowed, not mandatory. This routine
% never applies a layer-coverage repair and cannot represent excluded spots.

narginchk(2, 2);
required = {'doseModel','w0','minimumWeight','groups'};
for i = 1:numel(required)
    assert(isfield(problem, required{i}), 'Missing problem.%s.', required{i});
end
settings = defaults(settings);
w = max(0, reshape(double(problem.w0), [], 1));
n = numel(w);
m = expand(problem.minimumWeight, n);
z = project_mmu(w, m);
u = zeros(n, 1);
groups = problem.groups;
history = repmat(empty_row(), settings.maxIterations, 1);
supportHistory = cell(settings.maxIterations, 1);
totalCounts = struct('objective', 0, 'gradient', 0);
supportEvents = repmat(empty_support_event(), 0, 1);
timer = tic; stopReason = 'maximum_iterations';
rho = settings.rho;

for k = 1:settings.maxIterations
    zOld = z;
    supportOld = zOld > 0;
    center = z - u;
    inner = solve_w(problem.doseModel, w, center, rho, settings);
    w = max(0, reshape(double(inner.w), [], 1));
    uBefore = u;
    a = w + uBefore;
    [z, proxInfo] = prox_l0_mmu(a, settings.lambda0, rho, m);
    u = u + w - z;
    r = w - z;
    s = rho * (z - zOld);
    epsPri = sqrt(n) * settings.absoluteTolerance + ...
        settings.relativeTolerance * max(norm(w), norm(z));
    epsDual = sqrt(n) * settings.absoluteTolerance + ...
        settings.relativeTolerance * norm(rho * u);

    [fw, ~, dw] = evaluate_matrad_dose_model( ...
        problem.doseModel, w, false);
    [fz, ~, dz] = evaluate_matrad_dose_model( ...
        problem.doseModel, z, false);
    totalCounts.objective = totalCounts.objective + ...
        inner.callbackCounts.objective + 2;
    totalCounts.gradient = totalCounts.gradient + ...
        inner.callbackCounts.gradient;
    support = z > 0;
    supportHistory{k} = support;
    activeLayers = count_active_groups(support, groups);
    cyclingPeriod = detect_cycle(supportHistory, k);
    changed = find(support ~= supportOld);
    if ~isempty(changed)
        boundary = l0_mmu_switching_boundary(settings.lambda0, rho, m);
        supportEvents = [supportEvents; make_support_events( ...
            k, changed, supportOld, support, a, proxInfo, boundary, rho, ...
            problem)]; %#ok<AGROW>
    end
    bdt = NaN;
    if isfield(problem, 'bdtFunction') && ~isempty(problem.bdtFunction)
        bdt = double(problem.bdtFunction(z));
    end
    [rhoNext, uNext, rhoAction] = update_rho( ...
        rho, u, norm(r), norm(s), settings);
    history(k) = make_row(k, fw, fz, settings.lambda0 * nnz(support), ...
        norm(r), norm(s), epsPri, epsDual, support, activeLayers, ...
        proxInfo, dw, dz, inner, totalCounts, bdt, cyclingPeriod, ...
        rho, rhoNext, rhoAction);

    if settings.checkpointEvery > 0 && mod(k, settings.checkpointEvery) == 0 ...
            && isfield(settings, 'checkpointFunction') ...
            && ~isempty(settings.checkpointFunction)
        state = struct('w', w, 'z', z, 'u', u, 'iteration', k, ...
            'history', history(1:k), 'supportHistory', {supportHistory(1:k)}, ...
            'settings', settings, 'rho', rhoNext, ...
            'supportEvents', supportEvents);
        settings.checkpointFunction(state);
    end
    if norm(r) <= epsPri && norm(s) <= epsDual
        stopReason = 'admm_residual_converged';
        break;
    end
    if cyclingPeriod > 0 && k >= settings.minimumIterations
        stopReason = 'support_cycling';
        break;
    end
    rho = rhoNext;
    u = uNext;
end
history = history(1:k);
supportHistory = supportHistory(1:k);
result = struct('w', w, 'z', z, 'u', u, 'history', history, ...
    'supportHistory', {supportHistory}, 'stopReason', stopReason, ...
    'supportEvents', supportEvents, ...
    'iterations', k, 'wallSeconds', toc(timer), ...
    'callbackCounts', totalCounts, 'settings', settings, ...
    'mmuFeasible', all(z == 0 | z >= m), ...
    'activeSpots', nnz(z > 0), ...
    'activeLayers', count_active_groups(z > 0, groups));
end

function inner = solve_w(model, w, center, rho, settings)
switch lower(settings.wStrategy)
    case 'ipopt'
        inner = solve_admm_w_ipopt(model, w, center, rho, settings.ipopt);
    case 'pg'
        inner = solve_admm_w_first_order(model, w, center, rho, settings.inner);
    case 'lbfgs'
        inner = solve_admm_w_projected_lbfgs(model, w, center, rho, settings.inner);
    otherwise
        error('run_stage2_admm_variant_h:unknownStrategy', ...
            'Unknown wStrategy: %s.', settings.wStrategy);
end
end

function settings = defaults(settings)
base = struct('lambda0', 0, 'rho', 1, 'maxIterations', 10, ...
    'minimumIterations', 3, 'absoluteTolerance', 1e-5, ...
    'relativeTolerance', 1e-3, 'wStrategy', 'pg', ...
    'checkpointEvery', 0, 'inner', struct(), 'ipopt', struct(), ...
    'rhoStrategy', 'fixed', 'rhoBalance', struct('mu', 2, ...
    'increaseFactor', 2, 'decreaseFactor', 2, ...
    'minimum', 0.125, 'maximum', 8));
names = fieldnames(settings);
for i = 1:numel(names), base.(names{i}) = settings.(names{i}); end
settings = base;
assert(settings.lambda0 >= 0 && settings.rho > 0);
end

function [rhoNext, uNext, action] = update_rho(rho, u, r, s, settings)
switch lower(settings.rhoStrategy)
    case 'fixed'
        rhoNext = rho; uNext = u; action = 'fixed';
    case 'residual_balancing'
        [rhoNext, uNext, action] = balance_admm_rho( ...
            rho, u, r, s, settings.rhoBalance);
    otherwise
        error('run_stage2_admm_variant_h:unknownRhoStrategy', ...
            'Unknown rhoStrategy: %s.', settings.rhoStrategy);
end
end

function x = expand(x, n)
x = reshape(double(x), [], 1);
if isscalar(x), x = repmat(x, n, 1); end
assert(numel(x) == n && all(isfinite(x)) && all(x >= 0));
end

function z = project_mmu(a, m)
% Euclidean projection onto {0} union [m_i,Inf), with ties to zero.
z = prox_l0_mmu(a, 0, 1, m);
end

function nActive = count_active_groups(support, groups)
nActive = 0;
for g = 1:numel(groups)
    idx = group_indices(groups(g));
    nActive = nActive + any(support(idx));
end
end

function idx = group_indices(group)
if iscell(group), group = group{1}; end
if isstruct(group)
    if isfield(group, 'indices'), idx = group.indices;
    elseif isfield(group, 'spotIndices'), idx = group.spotIndices;
    elseif isfield(group, 'idx'), idx = group.idx;
    elseif isfield(group, 'ix'), idx = group.ix;
    else, error('Unknown group index field.'); end
else
    idx = group;
end
idx = reshape(double(idx), [], 1);
end

function period = detect_cycle(supportHistory, k)
period = 0;
for p = 2:min(4, floor(k/2))
    if isequal(supportHistory{k}, supportHistory{k-p})
        changedBetween = false;
        for j = k-p+1:k-1
            changedBetween = changedBetween || ...
                ~isequal(supportHistory{j}, supportHistory{k});
        end
        if changedBetween
            period = p;
            return;
        end
    end
end
end

function row = make_row(k, fw, fz, l0, r, s, ep, ed, support, layers, ...
    proxInfo, dw, dz, inner, counts, bdt, cycle, rho, rhoNext, rhoAction)
row = empty_row();
row.iteration = k; row.doseObjectiveW = fw; row.doseObjectiveZ = fz;
row.l0Term = l0; row.trueObjectiveZ = fz + l0;
row.primalResidual = r; row.dualResidual = s;
row.primalTolerance = ep; row.dualTolerance = ed;
row.activeSpots = nnz(support); row.exactZeroSpots = numel(support)-nnz(support);
row.activeLayers = layers; row.minimumPositiveZ = proxInfo.minimumPositive;
row.mmuFeasible = proxInfo.mmuFeasible;
row.constraintViolationW = dw.maxConstraintViolation;
row.constraintViolationZ = dz.maxConstraintViolation;
row.innerIterations = inner.iterations; row.innerWallSeconds = inner.wallSeconds;
row.innerReason = string(inner.reason);
row.objectiveEvaluations = counts.objective;
row.gradientEvaluations = counts.gradient;
row.provisionalBDT = bdt; row.supportCyclePeriod = cycle;
row.rho = rho; row.rhoNext = rhoNext; row.rhoAction = string(rhoAction);
end

function row = empty_row()
row = struct('iteration', NaN, 'doseObjectiveW', NaN, ...
    'doseObjectiveZ', NaN, 'l0Term', NaN, 'trueObjectiveZ', NaN, ...
    'primalResidual', NaN, 'dualResidual', NaN, ...
    'primalTolerance', NaN, 'dualTolerance', NaN, ...
    'activeSpots', NaN, 'exactZeroSpots', NaN, 'activeLayers', NaN, ...
    'minimumPositiveZ', NaN, 'mmuFeasible', false, ...
    'constraintViolationW', NaN, 'constraintViolationZ', NaN, ...
    'innerIterations', NaN, 'innerWallSeconds', NaN, 'innerReason', "", ...
    'objectiveEvaluations', NaN, 'gradientEvaluations', NaN, ...
    'provisionalBDT', NaN, 'supportCyclePeriod', NaN, ...
    'rho', NaN, 'rhoNext', NaN, 'rhoAction', "");
end

function events = make_support_events(k, changed, oldSupport, newSupport, ...
    a, proxInfo, boundary, rho, problem)
events = repmat(empty_support_event(), numel(changed), 1);
fullIndex = changed;
if isfield(problem, 'reducedToFull') && ~isempty(problem.reducedToFull)
    fullIndex = problem.reducedToFull(changed);
end
for j = 1:numel(changed)
    i = changed(j);
    events(j).iteration = k;
    events(j).reducedIndex = i;
    events(j).fullIndex = fullIndex(j);
    events(j).wasActive = oldSupport(i);
    events(j).isActive = newSupport(i);
    events(j).a = a(i);
    events(j).offCost = proxInfo.offCost(i);
    events(j).onCost = proxInfo.onCost(i);
    events(j).costGapOnMinusOff = proxInfo.deltaOnMinusOff(i);
    events(j).switchingBoundary = boundary(i);
    events(j).signedBoundaryDistance = a(i) - boundary(i);
    events(j).absoluteBoundaryDistance = abs(a(i) - boundary(i));
    events(j).relativeBoundaryDistance = abs(a(i)-boundary(i)) / ...
        max([1, abs(boundary(i))]);
    events(j).rho = rho;
end
end

function event = empty_support_event()
event = struct('iteration', NaN, 'reducedIndex', NaN, 'fullIndex', NaN, ...
    'wasActive', false, 'isActive', false, 'a', NaN, 'offCost', NaN, ...
    'onCost', NaN, 'costGapOnMinusOff', NaN, ...
    'switchingBoundary', NaN, 'signedBoundaryDistance', NaN, ...
    'absoluteBoundaryDistance', NaN, 'relativeBoundaryDistance', NaN, ...
    'rho', NaN);
end
