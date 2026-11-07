function report = run_integrated_reduced_validation(dij, cst, pln, denseW, ...
    map, cfg, machine, maxSpotsPerControlPoint, runFinalIpopt)
% run_integrated_reduced_validation  Validate native callbacks on small dij.
%
% This routine reuses an already available Dense Arc dij. It keeps the
% strongest Dense Arc spots per control point, validates dose equivalence and
% the matRad gradient, runs a short integrated solve, and optionally performs
% a small support-constrained final IPOPT solve. It never densifies dij.

narginchk(7, 9);
if nargin < 8 || isempty(maxSpotsPerControlPoint)
    maxSpotsPerControlPoint = 20;
end
if nargin < 9
    runFinalIpopt = false;
end
denseW = reshape(double(denseW), [], 1);
n = double(dij.totalNumOfBixels);
if numel(denseW) ~= n || numel(map) ~= n
    error('run_integrated_reduced_validation:dimensionMismatch', ...
        'dij, denseW and map dimensions must agree.');
end
[validationMask, selectionDiagnostics] = build_layer_stratified_support( ...
    denseW, map, maxSpotsPerControlPoint);

minimumRequestedLayers = 10;
if isfield(cfg.integrated, 'minimumValidationLayersPerCP')
    minimumRequestedLayers = cfg.integrated.minimumValidationLayersPerCP;
end
for k = 1:numel(selectionDiagnostics)
    required = min(minimumRequestedLayers, ...
        selectionDiagnostics(k).availableLayers);
    if selectionDiagnostics(k).selectedLayers < required
        error('run_integrated_reduced_validation:insufficientLayerDiversity', ...
            'CP %d retained %d layers; at least %d are required.', ...
            selectionDiagnostics(k).cpIndex, ...
            selectionDiagnostics(k).selectedLayers, required);
    end
end

[dijReduced, reducedToFull, fullToReduced, mapReduced] = ...
    reduce_dij_to_support(dij, validationMask, map);
for i = 1:numel(mapReduced)
    mapReduced(i).originalGlobalIx = mapReduced(i).globalIx;
    mapReduced(i).globalIx = i;
end
[groupsReduced, cpGroupsReduced] = build_group_index(mapReduced);
w0 = denseW(reducedToFull);

% D(:,S)*w_S must equal the same full-space multiplication exactly up to
% floating-point accumulation. No full dense matrix is constructed.
Dfull = nominal_matrix(dij);
Dreduced = nominal_matrix(dijReduced);
doseReduced = Dreduced * w0;
doseFromFullColumns = Dfull(:, reducedToFull) * w0;
doseEquivalenceRelativeError = norm(doseReduced - doseFromFullColumns) / ...
    max(1, norm(doseFromFullColumns));

doseModel = build_matrad_dose_model( ...
    dijReduced, cst, pln, cfg.integrated.constraintPenalty);
context = build_integrated_context( ...
    groupsReduced, cpGroupsReduced, mapReduced, machine, cfg);
[initialObjective, initialParts] = evaluate_integrated_objective( ...
    w0, doseModel, context, cfg);
runGradientCheck = true;
if isfield(cfg.integrated, 'runGradientCheck')
    runGradientCheck = logical(cfg.integrated.runGradientCheck);
end
if runGradientCheck
    gradientReport = validate_matrad_objective_gradient( ...
        doseModel, w0, [], []);
else
    gradientReport = struct('maximumRelativeError', NaN, ...
        'skipped', true);
end

validationCfg = configure_integrated_stage(cfg);
validationCfg.integrated.maxIterations = cfg.integrated.maxIterations;
optimizerResult = run_integrated_optimizer( ...
    doseModel, context, w0, validationCfg);
objectiveTrajectory = [initialObjective; optimizerResult.history.objective(:)];
descentTolerance = 1e-10 * max(1, max(abs(objectiveTrajectory)));
objectiveDescentPass = all(diff(objectiveTrajectory) <= descentTolerance);
exactLayerStats = exact_layer_stats(optimizerResult.w, groupsReduced);
supportReduced = extract_integrated_support( ...
    optimizerResult.w, w0, context, validationCfg);

supportFull = false(n, 1);
supportFull(reducedToFull(supportReduced.mask)) = true;
wOptimizerFull = zeros(n, 1);
wOptimizerFull(reducedToFull) = optimizerResult.w;
wSupportFull = zeros(n, 1);
wSupportFull(reducedToFull) = supportReduced.w;
if any(wSupportFull(~supportFull) ~= 0)
    error('run_integrated_reduced_validation:supportLeakage', ...
        'Reduced-to-full support embedding is inconsistent.');
end
supportEmbeddingPreserved = ~any(wSupportFull(~supportFull) ~= 0);

finalResult = [];
if runFinalIpopt
    finalCfg = validationCfg;
    finalCfg.enforceSelectedMinimumWeight = ...
        validationCfg.integrated.enforceFinalMinimumWeight;
    finalResult = run_final_reoptimization( ...
        dij, cst, pln, wSupportFull, finalCfg, supportFull, ...
        supportReduced.layerMask);
end

report = struct();
report.validationMask = validationMask;
report.selectionDiagnostics = selectionDiagnostics;
report.minimumSelectedLayersPerCP = min([selectionDiagnostics.selectedLayers]);
report.medianSelectedLayersPerCP = median([selectionDiagnostics.selectedLayers]);
report.maximumSelectedLayersPerCP = max([selectionDiagnostics.selectedLayers]);
report.reducedToFullIdx = reducedToFull;
report.fullToReducedIdx = fullToReduced;
report.nReducedSpots = nnz(validationMask);
report.doseEquivalenceRelativeError = doseEquivalenceRelativeError;
report.gradient = gradientReport;
report.initialObjective = initialObjective;
report.initialParts = initialParts;
report.optimizer = optimizerResult;
report.objectiveTrajectory = objectiveTrajectory;
report.objectiveDescentPass = objectiveDescentPass;
report.exactLayerStats = exactLayerStats;
report.constraintPass = optimizerResult.parts.maxConstraintViolation <= 1e-4;
report.supportEmbeddingPreserved = supportEmbeddingPreserved;
report.wInitialFull = zeros(n, 1);
report.wInitialFull(reducedToFull) = w0;
report.wOptimizerFull = wOptimizerFull;
report.wSupportFull = wSupportFull;
report.supportReduced = supportReduced;
report.supportFull = supportFull;
report.final = finalResult;
report.finalIpoptRun = logical(runFinalIpopt);
end

function stats = exact_layer_stats(w, groups)
w = reshape(double(w), [], 1);
nGroups = numel(groups);
groupNorm = zeros(nGroups, 1);
exactZeroMask = false(nGroups, 1);
for g = 1:nGroups
    values = w(groups(g).ix);
    groupNorm(g) = norm(values, 2);
    exactZeroMask(g) = all(values == 0);
end
stats = struct('nGroups', nGroups, ...
    'nExactZeroLayers', nnz(exactZeroMask), ...
    'nExactActiveLayers', nGroups - nnz(exactZeroMask), ...
    'exactZeroMask', exactZeroMask, 'groupNorm', groupNorm, ...
    'groupIndex', (1:nGroups)', ...
    'cpIndex', reshape(double([groups.cpIndex]), [], 1), ...
    'energyIndex', reshape(double([groups.energyIndex]), [], 1), ...
    'energy', reshape(double([groups.energy]), [], 1));
end

function D = nominal_matrix(dij)
if iscell(dij.physicalDose)
    D = dij.physicalDose{1};
else
    D = dij.physicalDose;
end
end
