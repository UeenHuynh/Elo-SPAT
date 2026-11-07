function result = run_integrated_optimizer(doseModel, context, wInit, cfg, resumeState)
% run_integrated_optimizer  MATLAB-native nonconvex proximal-gradient solver.
%
% Layer sparsity uses the exact group L_{2,1/2} proximal map. Smooth-L0
% spot activation is handled by majorization-minimization (reweighted L1).
% Backtracking is applied to the smooth matRad/ES/time/barrier component.

narginchk(4, 5);
if nargin < 5
    resumeState = [];
end
opts = cfg.integrated;
w = max(0, reshape(double(wInit), [], 1));
if numel(w) ~= context.nWeights || context.nWeights ~= doseModel.nWeights
    error('run_integrated_optimizer:dimensionMismatch', ...
        'wInit, mapping and matRad objective dimensions must agree.');
end
maxIter = opts.maxIterations;
history = initialize_history(maxIter);
previousHistory = initialize_history(0);
iterationOffset = 0;
resumeObjectiveRelativeError = 0;
if nargin < 5 || isempty(resumeState)
    step = opts.initialStep;
    [objective, parts] = evaluate_integrated_objective(w, doseModel, context, cfg);
else
    [step, objective, parts, previousHistory, iterationOffset, ...
        resumeObjectiveRelativeError] = initialize_resume_state( ...
        resumeState, w, doseModel, context, cfg);
end
[spotWeights, spotMMAnchorIteration] = initialize_spot_mm_state( ...
    w, cfg, resumeState, iterationOffset);
status = 'maximum_iterations';
stoppingReason = 'additional_iteration_limit_reached';
stationarityStoppingEnabled = isfield(opts, 'stationarityStoppingEnabled') && ...
    logical(opts.stationarityStoppingEnabled);

for k = 1:maxIter
    stepStartIteration = iterationOffset + k - 1;
    spotMMUpdated = false;
    if cfg.lambdaS > 0 && ...
            stepStartIteration - spotMMAnchorIteration >= spot_mm_period(cfg)
        spotWeights = smooth_l0_mm_weights(w, cfg);
        spotMMAnchorIteration = stepStartIteration;
        spotMMUpdated = true;
    elseif cfg.lambdaS > 0 && stepStartIteration == 0
        spotMMUpdated = true;
    end
    [smoothValue, gradient] = evaluate_integrated_smooth( ...
        w, doseModel, context, cfg, true);

    accepted = false;
    trialStep = min(step / opts.backtrackGrowth, opts.maximumStep);
    for bt = 1:opts.maxBacktracking
        candidate = prox_rms_elospat(w - trialStep * gradient, trialStep, ...
            context.groups, cfg.lambdaEL, spotWeights);
        delta = candidate - w;
        [candidateSmooth, ~, candidateParts] = evaluate_integrated_smooth( ...
            candidate, doseModel, context, cfg, false);
        upper = smoothValue + gradient' * delta + ...
            sum(delta .^ 2) / (2 * trialStep) + opts.sufficientDecreaseTolerance;
        if candidateSmooth <= upper
            [candidateObjective, candidateParts] = add_sparse_penalties( ...
                candidateSmooth, candidateParts, candidate, context, cfg);
            objectiveTolerance = opts.sufficientDecreaseTolerance + ...
                1e-12 * max(1, abs(objective));
            if candidateObjective <= objective + objectiveTolerance
                accepted = true;
                break;
            end
        end
        trialStep = trialStep * opts.backtrackFactor;
        if trialStep < opts.minimumStep
            break;
        end
    end
    if ~accepted
        status = 'line_search_failed';
        stoppingReason = 'backtracking_exhausted_without_descent';
        break;
    end

    relStep = norm(delta) / max(1, norm(w));
    relObjective = abs(candidateObjective - objective) / max(1, abs(objective));
    objectiveDecrease = objective - candidateObjective;
    stationarity = evaluate_stationarity( ...
        candidate, trialStep, doseModel, context, cfg);
    [nExactZeroLayers, nExactActiveLayers, candidateZeroMask] = ...
        exact_layer_counts(candidate, context.groups);
    [~, ~, previousZeroMask] = exact_layer_counts(w, context.groups);
    globalIteration = iterationOffset + k;
    history.iteration(k) = globalIteration;
    history.objective(k) = candidateObjective;
    history.doseObjective(k) = candidateParts.dose;
    history.layerPenalty(k) = candidateParts.layerL21half;
    history.spotPenalty(k) = candidateParts.spotSmoothL0;
    history.energySequencing(k) = candidateParts.energySequencing;
    history.timeSurrogate(k) = candidateParts.timeSurrogate;
    history.maxConstraintViolation(k) = candidateParts.maxConstraintViolation;
    history.relativeStep(k) = relStep;
    history.relativeObjective(k) = relObjective;
    history.stepSize(k) = trialStep;
    history.objectiveDecrease(k) = objectiveDecrease;
    history.proxDisplacementNorm(k) = stationarity.proxDisplacementNorm;
    history.proxGradientResidualNorm(k) = ...
        stationarity.proxGradientResidualNorm;
    history.normalizedStationarity(k) = stationarity.normalizedStationarity;
    history.gradientNorm(k) = stationarity.gradientNorm;
    history.exactZeroLayers(k) = nExactZeroLayers;
    history.exactActiveLayers(k) = nExactActiveLayers;
    history.activeSpots(k) = nnz(candidate > cfg.spotThreshold);
    history.layerSupportChanged(k) = any(candidateZeroMask ~= previousZeroMask);
    history.spotSupportChanged(k) = any((candidate == 0) ~= (w == 0));
    history.spotMMUpdated(k) = spotMMUpdated;
    history.spotMMAnchorIteration(k) = spotMMAnchorIteration;
    history.spotMMWeightNorm(k) = norm(spotWeights);
    history.backtrackingTrials(k) = bt;

    w = candidate;
    objective = candidateObjective;
    parts = candidateParts;
    step = trialStep;
    if opts.verbose && (k == 1 || mod(globalIteration, opts.printEvery) == 0)
        fprintf(['Integrated RMS-ELO-SPAT iter %d: J=%.8g, Fdose=%.8g, ' ...
            'step=%.3e, stationarity=%.3e, zeroLayers=%d, active=%d\n'], ...
            globalIteration, objective, parts.dose, trialStep, ...
            stationarity.normalizedStationarity, nExactZeroLayers, ...
            nnz(w > cfg.spotThreshold));
    end
    if stationarityStoppingEnabled && ...
            globalIteration >= opts.minimumIterations && ...
            stationarity.normalizedStationarity <= opts.stationarityTolerance
        status = 'stationary_converged';
        stoppingReason = 'calibrated_normalized_proximal_gradient_tolerance';
        break;
    end
    if stationarityStoppingEnabled && ...
            is_numerically_stagnant(previousHistory, history, k, opts) && ...
            stationarity.normalizedStationarity > opts.stationarityTolerance
        status = 'numerical_stagnation';
        stoppingReason = 'relative_changes_stagnant_above_stationarity_tolerance';
        break;
    end
end

nDone = find(history.iteration > 0, 1, 'last');
if isempty(nDone)
    nDone = 0;
end
segmentHistory = trim_history(history, nDone);
history = append_history(previousHistory, segmentHistory);
result = struct('w', w, 'objective', objective, 'parts', parts, ...
    'history', history, 'iterations', iterationOffset + nDone, ...
    'segmentIterations', nDone, 'resumedFromIteration', iterationOffset, ...
    'acceptedStep', step, 'resumeObjectiveRelativeError', ...
    resumeObjectiveRelativeError, 'status', status, ...
    'stoppingReason', stoppingReason, ...
    'stationarityStoppingEnabled', stationarityStoppingEnabled, ...
    'stationarityTolerance', opts.stationarityTolerance, ...
    'spotMMUpdatePeriod', spot_mm_period(cfg), ...
    'spotMMWeights', spotWeights, ...
    'spotMMAnchorIteration', spotMMAnchorIteration, ...
    'finalNormalizedStationarity', final_history_value( ...
        history, 'normalizedStationarity', Inf), ...
    'solver', 'MATLAB proximal-gradient with group-half prox and spot MM');
end

function [weights, anchorIteration] = initialize_spot_mm_state( ...
    w, cfg, resumeState, iterationOffset)
if cfg.lambdaS <= 0
    weights = zeros(size(w));
    anchorIteration = iterationOffset;
    return;
end
period = spot_mm_period(cfg);
if isempty(resumeState)
    weights = smooth_l0_mm_weights(w, cfg);
    anchorIteration = 0;
    return;
end
if isfield(resumeState, 'spotMMUpdatePeriod') && ...
        double(resumeState.spotMMUpdatePeriod) ~= period
    error('run_integrated_optimizer:spotMMPeriodMismatch', ...
        'Resume checkpoint MM period differs from current configuration.');
end
if isfield(resumeState, 'spotMMWeights') && ...
        isfield(resumeState, 'spotMMAnchorIteration')
    weights = reshape(double(resumeState.spotMMWeights), [], 1);
    anchorIteration = double(resumeState.spotMMAnchorIteration);
    if numel(weights) ~= numel(w) || any(~isfinite(weights)) || ...
            any(weights < 0) || anchorIteration > iterationOffset
        error('run_integrated_optimizer:invalidSpotMMResumeState', ...
            'Saved spot-MM state is invalid.');
    end
elseif period == 1
    % Historical checkpoints recomputed q at every iteration, so rebuilding
    % it from the saved w is an exact continuation for period one.
    weights = smooth_l0_mm_weights(w, cfg);
    anchorIteration = iterationOffset;
else
    error('run_integrated_optimizer:missingSpotMMResumeState', ...
        'Block-MM continuation requires saved weights and anchor iteration.');
end
end

function weights = smooth_l0_mm_weights(w, cfg)
weights = cfg.lambdaS * cfg.epsS ./ (w + cfg.epsS) .^ 2;
end

function period = spot_mm_period(cfg)
if isfield(cfg.integrated, 'spotMMUpdatePeriod')
    period = double(cfg.integrated.spotMMUpdatePeriod);
else
    period = 1;
end
if ~isscalar(period) || ~isfinite(period) || period < 1 || period ~= round(period)
    error('run_integrated_optimizer:invalidSpotMMUpdatePeriod', ...
        'spotMMUpdatePeriod must be a positive integer.');
end
end

function [step, objective, parts, previousHistory, iterationOffset, ...
    relativeError] = initialize_resume_state( ...
    state, w, doseModel, context, cfg)
required = {'w','objective','parts','history','iterations'};
for i = 1:numel(required)
    if ~isfield(state, required{i})
        error('run_integrated_optimizer:invalidResumeState', ...
            'resumeState is missing %s.', required{i});
    end
end
stateW = reshape(double(state.w), [], 1);
if numel(stateW) ~= numel(w) || any(stateW ~= w)
    error('run_integrated_optimizer:resumeWeightMismatch', ...
        'wInit must exactly equal resumeState.w.');
end
previousHistory = state.history;
iterationOffset = double(state.iterations);
if iterationOffset ~= numel(previousHistory.iteration) || ...
        (iterationOffset > 0 && previousHistory.iteration(end) ~= iterationOffset)
    error('run_integrated_optimizer:invalidResumeHistory', ...
        'resumeState iteration count and history are inconsistent.');
end
if isfield(state, 'acceptedStep') && isfinite(state.acceptedStep)
    step = double(state.acceptedStep);
elseif iterationOffset > 0
    step = double(previousHistory.stepSize(end));
else
    step = opts_initial_step(cfg);
end
if step < cfg.integrated.minimumStep || step > cfg.integrated.maximumStep
    error('run_integrated_optimizer:invalidResumeStep', ...
        'Saved accepted step is outside configured bounds.');
end
objective = double(state.objective);
parts = state.parts;
[checkObjective, ~] = evaluate_integrated_objective(w, doseModel, context, cfg);
relativeError = abs(checkObjective - objective) / max(1, abs(objective));
if relativeError > 1e-11
    error('run_integrated_optimizer:resumeObjectiveMismatch', ...
        'Saved and recomputed objectives differ by %.3e.', relativeError);
end
end

function value = opts_initial_step(cfg)
value = cfg.integrated.initialStep;
end

function diagnostic = evaluate_stationarity(w, step, doseModel, context, cfg)
% Fixed-point residual for the current nonconvex prox/MM subproblem.
[~, gradient] = evaluate_integrated_smooth( ...
    w, doseModel, context, cfg, true);
if cfg.lambdaS > 0
    spotWeights = cfg.lambdaS * cfg.epsS ./ (w + cfg.epsS) .^ 2;
else
    spotWeights = zeros(size(w));
end
proxPoint = prox_rms_elospat(w - step * gradient, step, ...
    context.groups, cfg.lambdaEL, spotWeights);
displacementNorm = norm(w - proxPoint);
residualNorm = displacementNorm / step;
diagnostic = struct('gradientNorm', norm(gradient), ...
    'proxDisplacementNorm', displacementNorm, ...
    'proxGradientResidualNorm', residualNorm, ...
    'normalizedStationarity', residualNorm / max(1, norm(w)));
end

function [nZero, nActive, zeroMask] = exact_layer_counts(w, groups)
nZero = 0;
zeroMask = false(numel(groups), 1);
for g = 1:numel(groups)
    zeroMask(g) = all(w(groups(g).ix) == 0);
    nZero = nZero + zeroMask(g);
end
nActive = numel(groups) - nZero;
end

function tf = is_numerically_stagnant(previousHistory, history, k, opts)
tf = false;
globalIteration = numel(previousHistory.iteration) + k;
if globalIteration < max(opts.minimumIterations, opts.stagnationWindow)
    return;
end
relativeStep = [previousHistory.relativeStep(:); history.relativeStep(1:k)];
relativeObjective = [previousHistory.relativeObjective(:); ...
    history.relativeObjective(1:k)];
ix = (numel(relativeStep) - opts.stagnationWindow + 1):numel(relativeStep);
tf = all(relativeStep(ix) <= ...
        opts.stagnationRelativeStepTolerance) && ...
    all(relativeObjective(ix) <= ...
        opts.stagnationRelativeObjectiveTolerance);
end

function value = final_history_value(history, field, fallback)
if isfield(history, field) && ~isempty(history.(field))
    value = history.(field)(end);
else
    value = fallback;
end
end

function [value, parts] = add_sparse_penalties(smoothValue, parts, w, context, cfg)
layerPenalty = 0;
for g = 1:numel(context.groups)
    layerPenalty = layerPenalty + sqrt(norm(w(context.groups(g).ix), 2));
end
spotPenalty = sum(w ./ (w + cfg.epsS));
value = smoothValue + cfg.lambdaEL * layerPenalty + cfg.lambdaS * spotPenalty;
parts.layerL21half = layerPenalty;
parts.spotSmoothL0 = spotPenalty;
parts.total = value;
end

function history = initialize_history(n)
fields = {'iteration', 'objective', 'doseObjective', 'layerPenalty', ...
    'spotPenalty', 'energySequencing', 'timeSurrogate', ...
    'maxConstraintViolation', 'relativeStep', 'relativeObjective', 'stepSize'};
fields = [fields, {'objectiveDecrease', 'proxDisplacementNorm', ...
    'proxGradientResidualNorm', 'normalizedStationarity', 'gradientNorm', ...
    'exactZeroLayers', 'exactActiveLayers', 'activeSpots', ...
    'layerSupportChanged', 'spotSupportChanged', 'spotMMUpdated', ...
    'spotMMAnchorIteration', 'spotMMWeightNorm'}];
fields{end+1} = 'backtrackingTrials';
for i = 1:numel(fields)
    history.(fields{i}) = zeros(n, 1);
end
end

function history = trim_history(history, n)
fields = fieldnames(history);
for i = 1:numel(fields)
    history.(fields{i}) = history.(fields{i})(1:n);
end
end

function combined = append_history(first, second)
combined = first;
fields = fieldnames(second);
for i = 1:numel(fields)
    if ~isfield(combined, fields{i})
        % Historical checkpoints predate this diagnostic. NaN deliberately
        % prevents those unobserved iterations from proving support stability.
        combined.(fields{i}) = NaN(numel(first.iteration), 1);
    end
    combined.(fields{i}) = [combined.(fields{i})(:); second.(fields{i})(:)];
end
end
