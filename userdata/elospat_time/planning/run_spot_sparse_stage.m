function [result, spotMask] = run_spot_sparse_stage(dij, cst, pln, wInit, map, groups, cfg, machine, captureSnapshots)
% run_spot_sparse_stage  Iterative spot sparsity and MMU-aware filtering.

narginchk(8, 9);
if nargin < 9
    captureSnapshots = false;
end
captureSnapshots = logical(captureSnapshots);

validate_stage_inputs(dij, cst, pln, wInit, map, groups, cfg, machine);

nBixels = double(dij.totalNumOfBixels);
wFull = reshape(double(wInit), [], 1);
D = get_nominal_dose_matrix(dij);
groupOfBixel = build_group_lookup(groups, nBixels);

spotRemovalPercentile = get_optional_wrapper_value(cfg, 'spotRemovalPercentile', 20);
relTol = double(cfg.wrapper.relTol);
maxIterations = max(1, round(double(cfg.wrapper.maxSpotIterations)));
toleranceFraction = double(cfg.evaluation.ctvD95TolerancePercent) / 100;

prevDoseMetrics = compute_dvh_metrics(D * wFull, [], cst, cfg);
prevCtvD95 = prevDoseMetrics.summary.targetD95Min;

bestBdt = estimate_bdt(extract_delivery_sequence(wFull, map, groups, cfg), machine);
bestScore = bestBdt.T_total;
bestW = wFull;
bestMask = wFull > cfg.spotThreshold;
history = repmat(struct( ...
    'iteration', 0, ...
    'nActiveSpots', 0, ...
    'bdtTotal', NaN, ...
    'ctvD95', NaN, ...
    'nRemovedSpots', 0, ...
    'preReoptBdtTotal', NaN, ...
    'preReoptCtvD95', NaN, ...
    'preReoptActiveSpots', 0), 0, 1);
snapshots = repmat(struct('iteration', 0, 'nActiveSpots', 0, ...
    'w', [], 'supportMask', []), 0, 1);

prevBdtTotal = bestBdt.T_total;

for iter = 1:maxIterations
    activeMask = wFull > cfg.spotThreshold;
    activeIdx = find(activeMask);
    if isempty(activeIdx)
        break;
    end

    spotDoseImportance = zeros(nBixels, 1);
    spotBdtSaving = zeros(nBixels, 1);
    spotRemovalScore = inf(nBixels, 1);

    for k = 1:numel(activeIdx)
        idx = activeIdx(k);
        spotDoseImportance(idx) = norm(D(:, idx) * wFull(idx), 2);
        spotBdtSaving(idx) = machine.spot.switchTime + ...
            wFull(idx) / max(machine.irradiation.muRate, eps);

        groupIx = groupOfBixel(idx);
        if groupIx > 0
            groupWeights = wFull(groups(groupIx).ix);
            if sum(groupWeights > cfg.spotThreshold) <= 1
                spotBdtSaving(idx) = spotBdtSaving(idx) + estimate_layer_switch_saving(groups(groupIx), machine);
            end
        end

        spotRemovalScore(idx) = compute_bdt_aware_removal_score( ...
            spotDoseImportance(idx), spotBdtSaving(idx), 1e-12);
    end

    activeScores = spotRemovalScore(activeIdx);
    cutoff = compute_percentile(activeScores, spotRemovalPercentile);
    subMmuMask = (wFull > 0) & (wFull < cfg.mmuThreshold);
    lowScoreMask = activeMask & (spotRemovalScore <= cutoff);
    candidateMask = subMmuMask | lowScoreMask;
    candidateMask = candidateMask & activeMask;

    if ~any(candidateMask)
        break;
    end

    spotMaskCandidate = activeMask & ~candidateMask;
    if ~any(spotMaskCandidate)
        break;
    end

    % Capture the exact candidate-support state before the reduced solver.
    % This is needed to diagnose whether BDT/dose changes come from masking
    % or from re-optimization itself.
    wPreReopt = wFull;
    wPreReopt(~spotMaskCandidate) = 0;
    seqPreReopt = extract_delivery_sequence(wPreReopt, map, groups, cfg);
    if isempty(seqPreReopt)
        bdtPreReopt = struct('T_total', NaN);
    else
        bdtPreReopt = estimate_bdt(seqPreReopt, machine);
    end
    doseMetricsPreReopt = compute_dvh_metrics(D * wPreReopt, [], cst, cfg);
    ctvD95PreReopt = doseMetricsPreReopt.summary.targetD95Min;

    [dijReduced, ~, wReduced] = cut_dij_by_mask(dij, spotMaskCandidate, [], wFull);
    [resultReduced, optimizerReduced] = reoptimize_reduced_problem(dijReduced, cst, pln, wReduced); %#ok<NASGU>

    reducedWeights = extract_result_weights(resultReduced, numel(wReduced));
    wEmbedded = zeros(nBixels, 1);
    wEmbedded(spotMaskCandidate) = reducedWeights;

    doseMetricsNew = compute_dvh_metrics(D * wEmbedded, [], cst, cfg);
    newCtvD95 = doseMetricsNew.summary.targetD95Min;
    if isfinite(prevCtvD95) && prevCtvD95 > 0 && isfinite(newCtvD95)
        degradation = (prevCtvD95 - newCtvD95) / prevCtvD95;
        if degradation > toleranceFraction
            break;
        end
    end

    seqNew = extract_delivery_sequence(wEmbedded, map, groups, cfg);
    if isempty(seqNew)
        break;
    end
    bdtNew = estimate_bdt(seqNew, machine);

    history(end + 1, 1) = struct( ... %#ok<AGROW>
        'iteration', iter, ...
        'nActiveSpots', sum(spotMaskCandidate), ...
        'bdtTotal', bdtNew.T_total, ...
        'ctvD95', newCtvD95, ...
        'nRemovedSpots', sum(candidateMask), ...
        'preReoptBdtTotal', bdtPreReopt.T_total, ...
        'preReoptCtvD95', ctvD95PreReopt, ...
        'preReoptActiveSpots', sum(spotMaskCandidate));

    if captureSnapshots
        snapshots(end + 1, 1) = struct( ... %#ok<AGROW>
            'iteration', iter, ...
            'nActiveSpots', sum(spotMaskCandidate), ...
            'w', wEmbedded, ...
            'supportMask', spotMaskCandidate);
    end

    if bdtNew.T_total < bestScore
        bestScore = bdtNew.T_total;
        bestW = wEmbedded;
        bestMask = spotMaskCandidate;
    end

    relChange = abs(bdtNew.T_total - prevBdtTotal) / max(abs(prevBdtTotal), eps);
    wFull = wEmbedded;
    prevBdtTotal = bdtNew.T_total;
    prevCtvD95 = newCtvD95;

    if relChange < relTol
        break;
    end
end

spotMask = reshape(logical(bestMask), [], 1);
result = struct();
result.w = reshape(double(bestW), [], 1);
result.spotMask = spotMask;
result.bdtMetrics = estimate_bdt(extract_delivery_sequence(bestW, map, groups, cfg), machine);
result.history = history;
if captureSnapshots
    result.supportSnapshots = snapshots;
end
end

function validate_stage_inputs(dij, cst, pln, wInit, map, groups, cfg, machine)
if ~isstruct(dij) || ~isfield(dij, 'totalNumOfBixels')
    error('run_spot_sparse_stage:invalidInput', ...
        'dij.totalNumOfBixels is required.');
end
if ~iscell(cst) || size(cst, 2) < 6
    error('run_spot_sparse_stage:invalidInput', ...
        'cst must be a valid matRad cst cell array.');
end
if ~isstruct(pln)
    error('run_spot_sparse_stage:invalidInput', ...
        'pln must be a struct.');
end
if ~isnumeric(wInit) || ~isvector(wInit) || numel(wInit) ~= double(dij.totalNumOfBixels)
    error('run_spot_sparse_stage:invalidInput', ...
        'wInit must match dij.totalNumOfBixels.');
end
if ~isstruct(map) || numel(map) ~= double(dij.totalNumOfBixels)
    error('run_spot_sparse_stage:invalidInput', ...
        'map must match dij.totalNumOfBixels.');
end
if ~isstruct(groups) || isempty(groups) || ~all(isfield(groups, {'ix', 'cpIndex'}))
    error('run_spot_sparse_stage:invalidInput', ...
        'groups must expose ix/cpIndex.');
end
if ~isstruct(cfg) || ~isfield(cfg, 'wrapper') || ~isfield(cfg, 'spotThreshold') || ...
        ~isfield(cfg, 'mmuThreshold') || ~isfield(cfg, 'evaluation')
    error('run_spot_sparse_stage:invalidInput', ...
        'cfg.wrapper, cfg.spotThreshold, cfg.mmuThreshold, and cfg.evaluation are required.');
end
if ~isstruct(machine)
    error('run_spot_sparse_stage:invalidInput', ...
        'machine must be a machine timing struct.');
end
end

function D = get_nominal_dose_matrix(dij)
if iscell(dij.physicalDose)
    D = dij.physicalDose{1};
else
    D = dij.physicalDose;
end
end

function lookup = build_group_lookup(groups, nBixels)
lookup = zeros(nBixels, 1);
for g = 1:numel(groups)
    lookup(groups(g).ix) = g;
end
end

function value = estimate_layer_switch_saving(group, machine)
energyBase = machine.energy.downBase;
if isfield(group, 'energy') && isfinite(group.energy)
    energyBase = energyBase + abs(double(group.energy)) * machine.energy.perMeVDown * 0;
end
value = energyBase + machine.buffer.perLayer;
end

function value = get_optional_wrapper_value(cfg, fieldName, defaultValue)
value = defaultValue;
if isfield(cfg, 'wrapper') && isstruct(cfg.wrapper) && isfield(cfg.wrapper, fieldName)
    candidate = cfg.wrapper.(fieldName);
    if isnumeric(candidate) && isscalar(candidate) && isfinite(candidate) && candidate > 0
        value = double(candidate);
    end
end
end

function p = compute_percentile(values, percentileValue)
values = reshape(double(values), [], 1);
values = values(isfinite(values));
if isempty(values)
    p = NaN;
    return;
end
values = sort(values, 'ascend');
n = numel(values);
if n == 1
    p = values;
    return;
end
position = 1 + (n - 1) * (percentileValue / 100);
lowerIx = floor(position);
upperIx = ceil(position);
if lowerIx == upperIx
    p = values(lowerIx);
else
    alpha = position - lowerIx;
    p = (1 - alpha) * values(lowerIx) + alpha * values(upperIx);
end
end

function wResult = extract_result_weights(resultReduced, expectedLength)
if ~isstruct(resultReduced) || ~isfield(resultReduced, 'w') || isempty(resultReduced.w)
    error('run_spot_sparse_stage:missingWeights', ...
        'Reduced optimization result did not return result.w.');
end
wResult = reshape(double(resultReduced.w), [], 1);
if numel(wResult) ~= expectedLength
    error('run_spot_sparse_stage:dimensionMismatch', ...
        'Reduced result.w length %d does not match expected length %d.', ...
        numel(wResult), expectedLength);
end
end
