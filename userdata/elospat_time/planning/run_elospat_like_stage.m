function [result, layerMask] = run_elospat_like_stage(dij, cst, pln, wInit, groups, cpGroups, cfg, machine)
% run_elospat_like_stage  ELO-SPAT-inspired wrapper layer filtering stage.
%
% This is a research-wrapper approximation around matRad's nominal solver.

narginchk(8, 8);

validate_stage_inputs(dij, cst, pln, wInit, groups, cpGroups, cfg, machine);

nBixels = double(dij.totalNumOfBixels);
wFull = reshape(double(wInit), [], 1);
D = get_nominal_dose_matrix(dij);
maxIterations = max(1, round(double(cfg.wrapper.maxLayerIterations)));
relTol = double(cfg.wrapper.relTol);
removalFraction = get_optional_wrapper_value(cfg, 'removalFraction', 0.10);

bestScore = inf;
bestResult = struct();
bestMask = true(numel(groups), 1);
bestBixelMask = true(nBixels, 1);
prevScore = NaN;
history = repmat(struct( ...
    'iteration', 0, ...
    'nActiveLayers', 0, ...
    'score', NaN, ...
    'bdtTotal', NaN, ...
    'removedGroups', []), 0, 1);

for iter = 1:maxIterations
    activeLayerMask = build_layer_mask(wFull, groups, cfg.layerThreshold);
    activeGroups = find(activeLayerMask);
    if isempty(activeGroups)
        break;
    end

    currentSeq = build_group_delivery_sequence(wFull, groups, cpGroups, cfg);
    if isempty(currentSeq)
        break;
    end

    currentBdt = estimate_bdt(currentSeq, machine);
    [currentScore, currentComponents] = composite_score(wFull, groups, cpGroups, cfg, currentBdt); %#ok<ASGLU>

    if currentScore < bestScore
        bestScore = currentScore;
        bestMask = activeLayerMask;
        bestBixelMask = expand_layer_mask_to_bixel_mask(bestMask, groups, nBixels);
        bestResult = build_result_struct(bestResult, wFull, bestMask, bestBixelMask, iter, currentScore, currentBdt, history);
    end

    doseContribution = zeros(numel(groups), 1);
    bdtSaving = zeros(numel(groups), 1);
    removalScore = inf(numel(groups), 1);

    for k = 1:numel(activeGroups)
        g = activeGroups(k);
        ix = groups(g).ix;
        doseContribution(g) = norm(D(:, ix) * wFull(ix), 2);

        wWithout = wFull;
        wWithout(ix) = 0;
        seqWithout = build_group_delivery_sequence(wWithout, groups, cpGroups, cfg);
        if isempty(seqWithout)
            bdtSaving(g) = currentBdt.T_total;
        else
            bdtWithout = estimate_bdt(seqWithout, machine);
            bdtSaving(g) = max(currentBdt.T_total - bdtWithout.T_total, 0);
        end

        removalScore(g) = compute_bdt_aware_removal_score(doseContribution(g), bdtSaving(g), 1e-12);
    end

    candidateGroups = choose_groups_for_removal(activeGroups, removalScore, groups, removalFraction);
    if isempty(candidateGroups)
        break;
    end

    layerMaskCandidate = activeLayerMask;
    layerMaskCandidate(candidateGroups) = false;
    bixelMask = expand_layer_mask_to_bixel_mask(layerMaskCandidate, groups, nBixels);
    [dijReduced, ~, wReduced] = cut_dij_by_mask(dij, bixelMask, [], wFull);
    [resultReduced, optimizerReduced] = reoptimize_reduced_problem(dijReduced, cst, pln, wReduced); %#ok<NASGU>

    reducedWeights = extract_result_weights(resultReduced, numel(wReduced));
    wEmbedded = zeros(nBixels, 1);
    wEmbedded(bixelMask) = reducedWeights;
    wFull = wEmbedded;

    seqReduced = build_group_delivery_sequence(wFull, groups, cpGroups, cfg);
    if isempty(seqReduced)
        reducedBdt = currentBdt;
    else
        reducedBdt = estimate_bdt(seqReduced, machine);
    end

    [reducedScore, ~] = composite_score(wFull, groups, cpGroups, cfg, reducedBdt);

    history(end + 1, 1) = struct( ... %#ok<AGROW>
        'iteration', iter, ...
        'nActiveLayers', sum(layerMaskCandidate), ...
        'score', reducedScore, ...
        'bdtTotal', reducedBdt.T_total, ...
        'removedGroups', reshape(candidateGroups, [], 1));

    if reducedScore < bestScore
        bestScore = reducedScore;
        bestMask = layerMaskCandidate;
        bestBixelMask = bixelMask;
        bestResult = build_result_struct(resultReduced, wFull, bestMask, bestBixelMask, iter, reducedScore, reducedBdt, history);
    end

    if iter >= 3 && isfinite(prevScore) && prevScore ~= 0
        relChange = abs(reducedScore - prevScore) / abs(prevScore);
        if relChange < relTol
            break;
        end
    end
    prevScore = reducedScore;
end

if isempty(fieldnames(bestResult))
    seqFinal = build_group_delivery_sequence(wFull, groups, cpGroups, cfg);
    if isempty(seqFinal)
        finalBdt = struct('T_total', 0);
        finalScore = inf;
    else
        finalBdt = estimate_bdt(seqFinal, machine);
        [finalScore, ~] = composite_score(wFull, groups, cpGroups, cfg, finalBdt);
    end
    layerMask = build_layer_mask(wFull, groups, cfg.layerThreshold);
    bestBixelMask = expand_layer_mask_to_bixel_mask(layerMask, groups, nBixels);
    result = build_result_struct(struct(), wFull, layerMask, bestBixelMask, 0, finalScore, finalBdt, history);
else
    result = bestResult;
    layerMask = bestMask;
end

end

function validate_stage_inputs(dij, cst, pln, wInit, groups, cpGroups, cfg, machine)
if ~isstruct(dij) || ~isfield(dij, 'totalNumOfBixels')
    error('run_elospat_like_stage:invalidInput', ...
        'dij.totalNumOfBixels is required.');
end
if ~iscell(cst) || size(cst, 2) < 6
    error('run_elospat_like_stage:invalidInput', ...
        'cst must be a valid matRad cst cell array.');
end
if ~isstruct(pln)
    error('run_elospat_like_stage:invalidInput', ...
        'pln must be a struct.');
end
if ~isnumeric(wInit) || ~isvector(wInit) || numel(wInit) ~= double(dij.totalNumOfBixels)
    error('run_elospat_like_stage:invalidInput', ...
        'wInit must match dij.totalNumOfBixels.');
end
if ~isstruct(groups) || isempty(groups) || ~all(isfield(groups, {'ix', 'cpIndex', 'energy'}))
    error('run_elospat_like_stage:invalidInput', ...
        'groups must expose ix/cpIndex/energy.');
end
if ~isstruct(cpGroups) || isempty(cpGroups) || ~all(isfield(cpGroups, {'cpIndex', 'gantryAngle'}))
    error('run_elospat_like_stage:invalidInput', ...
        'cpGroups must expose cpIndex/gantryAngle.');
end
if ~isstruct(cfg) || ~isfield(cfg, 'wrapper') || ~isfield(cfg, 'layerThreshold')
    error('run_elospat_like_stage:invalidInput', ...
        'cfg.wrapper and cfg.layerThreshold are required.');
end
if ~isstruct(machine)
    error('run_elospat_like_stage:invalidInput', ...
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

function value = get_optional_wrapper_value(cfg, fieldName, defaultValue)
value = defaultValue;
if isfield(cfg, 'wrapper') && isstruct(cfg.wrapper) && isfield(cfg.wrapper, fieldName)
    candidate = cfg.wrapper.(fieldName);
    if isnumeric(candidate) && isscalar(candidate) && isfinite(candidate) && candidate > 0
        value = double(candidate);
    end
end
end

function candidateGroups = choose_groups_for_removal(activeGroups, removalScore, groups, removalFraction)
nActive = numel(activeGroups);
cpList = unique([groups(activeGroups).cpIndex]', 'stable');
maxRemovable = nActive - numel(cpList);
if maxRemovable <= 0
    candidateGroups = [];
    return;
end

nToRemove = max(1, ceil(nActive * removalFraction));
nToRemove = min(nToRemove, maxRemovable);
[~, order] = sort(removalScore(activeGroups), 'ascend');
rankedGroups = activeGroups(order);

candidateGroups = zeros(0, 1);
keptPerCp = count_active_layers_per_cp(activeGroups, groups, cpList);

for i = 1:numel(rankedGroups)
    g = rankedGroups(i);
    cpIndex = groups(g).cpIndex;
    cpLoc = find(cpList == cpIndex, 1, 'first');
    if isempty(cpLoc) || keptPerCp(cpLoc) <= 1
        continue;
    end
    candidateGroups(end + 1, 1) = g; %#ok<AGROW>
    keptPerCp(cpLoc) = keptPerCp(cpLoc) - 1;
    if numel(candidateGroups) >= nToRemove
        break;
    end
end
end

function counts = count_active_layers_per_cp(activeGroups, groups, cpList)
counts = zeros(numel(cpList), 1);
for i = 1:numel(activeGroups)
    cpLoc = find(cpList == groups(activeGroups(i)).cpIndex, 1, 'first');
    if ~isempty(cpLoc)
        counts(cpLoc) = counts(cpLoc) + 1;
    end
end
end

function seq = build_group_delivery_sequence(w, groups, cpGroups, cfg)
layerThreshold = double(cfg.layerThreshold);
spotThreshold = double(cfg.spotThreshold);

seq = repmat(struct( ...
    'cpIndex', 0, ...
    'gantryAngle', 0, ...
    'energy', 0, ...
    'energyIndex', 0, ...
    'spotIndices', [], ...
    'spotWeights', [], ...
    'numActiveSpots', 0, ...
    'layerWeightSum', 0, ...
    'prevEnergy', [], ...
    'nextEnergy', []), 0, 1);

cpLookup = containers.Map('KeyType', 'double', 'ValueType', 'double');
for i = 1:numel(cpGroups)
    cpLookup(double(cpGroups(i).cpIndex)) = double(cpGroups(i).gantryAngle);
end

activeEntries = repmat(seq, 0, 1);
for g = 1:numel(groups)
    ix = reshape(double(groups(g).ix), [], 1);
    groupWeights = reshape(double(w(ix)), [], 1);
    if norm(groupWeights, 2) <= layerThreshold
        continue;
    end

    activeMask = groupWeights > spotThreshold;
    activeSpots = ix(activeMask);
    activeSpotWeights = groupWeights(activeMask);
    if ~isKey(cpLookup, double(groups(g).cpIndex))
        continue;
    end

    activeEntries(end + 1, 1) = struct( ... %#ok<AGROW>
        'cpIndex', double(groups(g).cpIndex), ...
        'gantryAngle', cpLookup(double(groups(g).cpIndex)), ...
        'energy', double(groups(g).energy), ...
        'energyIndex', double(groups(g).energyIndex), ...
        'spotIndices', activeSpots, ...
        'spotWeights', activeSpotWeights, ...
        'numActiveSpots', numel(activeSpots), ...
        'layerWeightSum', sum(activeSpotWeights), ...
        'prevEnergy', [], ...
        'nextEnergy', []);
    end

if isempty(activeEntries)
    return;
end

sortMatrix = [[activeEntries.gantryAngle]', -[activeEntries.energy]'];
[~, order] = sortrows(sortMatrix, [1 2]);
seq = activeEntries(order);

for i = 1:numel(seq)
    if i > 1
        seq(i).prevEnergy = seq(i - 1).energy;
    end
    if i < numel(seq)
        seq(i).nextEnergy = seq(i + 1).energy;
    end
end
end

function wResult = extract_result_weights(resultReduced, expectedLength)
if ~isstruct(resultReduced) || ~isfield(resultReduced, 'w') || isempty(resultReduced.w)
    error('run_elospat_like_stage:missingWeights', ...
        'Reduced optimization result did not return result.w.');
end
wResult = reshape(double(resultReduced.w), [], 1);
if numel(wResult) ~= expectedLength
    error('run_elospat_like_stage:dimensionMismatch', ...
        'Reduced result.w length %d does not match expected length %d.', ...
        numel(wResult), expectedLength);
end
end

function result = build_result_struct(resultBase, wFull, layerMask, bixelMask, iter, score, bdtMetrics, history)
result = resultBase;
result.w = reshape(double(wFull), [], 1);
result.layerMask = reshape(logical(layerMask), [], 1);
result.bixelMask = reshape(logical(bixelMask), [], 1);
result.wrapperIteration = iter;
result.compositeScore = score;
result.bdtMetrics = bdtMetrics;
result.history = history;
end
