function robustEval = evaluate_robust_scenarios(ct, cst, stf, dij, pln, w, cfg)
% evaluate_robust_scenarios  Evaluate dose metrics across available scenarios.
%
% INPUTS
%   ct, cst, stf, dij, pln, w, cfg - standard RMS-ELO-SPAT objects
%
% OUTPUT
%   robustEval - struct with per-scenario metrics, worst-case summary,
%                scenario pass-rate, and dose standard deviation cube

narginchk(7, 7);

if ~iscell(cst) || size(cst, 2) < 6
    error('evaluate_robust_scenarios:invalidInput', ...
        'cst must be a matRad cst cell array.');
end

if ~isstruct(dij) || ~isfield(dij, 'physicalDose')
    error('evaluate_robust_scenarios:invalidInput', ...
        'dij.physicalDose is required.');
end

if ~isnumeric(w) || isempty(w) || ~isvector(w)
    error('evaluate_robust_scenarios:invalidInput', ...
        'w must be a non-empty numeric vector.');
end

scenarioDoseFields = get_scenario_dose_fields(dij);
nScenarios = numel(scenarioDoseFields);
doseDims = infer_dose_dimensions(dij, ct, scenarioDoseFields{1});
nVoxels = prod(double(doseDims));

scenarioMetrics = repmat(struct( ...
    'scenarioIndex', 0, ...
    'doseMetrics', struct(), ...
    'targetD95Min', NaN, ...
    'targetD98Min', NaN, ...
    'targetV95Min', NaN, ...
    'oarDmaxMax', NaN, ...
    'scenarioInfo', struct()), nScenarios, 1);

sumDose = zeros(doseDims);
sumSqDose = zeros(doseDims);

for s = 1:nScenarios
    doseVector = compute_scenario_dose_vector(scenarioDoseFields{s}, w);
    if numel(doseVector) ~= nVoxels
        error('evaluate_robust_scenarios:dimensionMismatch', ...
            'Scenario %d dose has %d voxels, expected %d.', ...
            s, numel(doseVector), nVoxels);
    end

    doseForMetrics = reshape(doseVector, doseDims);
    scenarioDoseMetrics = compute_dvh_metrics(doseForMetrics, ct, cst, cfg);

    scenarioMetrics(s).scenarioIndex = s;
    scenarioMetrics(s).doseMetrics = scenarioDoseMetrics;
    scenarioMetrics(s).targetD95Min = scenarioDoseMetrics.summary.targetD95Min;
    scenarioMetrics(s).targetD98Min = scenarioDoseMetrics.summary.targetD98Min;
    scenarioMetrics(s).targetV95Min = scenarioDoseMetrics.summary.targetV95Min;
    scenarioMetrics(s).oarDmaxMax = scenarioDoseMetrics.summary.oarDmaxMax;
    scenarioMetrics(s).scenarioInfo = describe_scenario(pln, s);

    sumDose = sumDose + doseForMetrics;
    sumSqDose = sumSqDose + doseForMetrics .^ 2;
end

worstCase = compute_worst_case_metrics(scenarioMetrics);
nominalTargetD95 = scenarioMetrics(1).targetD95Min;
toleranceFraction = get_robust_tolerance(cfg) / 100;
scenarioPass = true(nScenarios, 1);

if isfinite(nominalTargetD95) && nominalTargetD95 > 0
    for s = 1:nScenarios
        currentD95 = scenarioMetrics(s).targetD95Min;
        if isfinite(currentD95)
            degradation = (nominalTargetD95 - currentD95) / nominalTargetD95;
            scenarioPass(s) = degradation <= toleranceFraction;
        else
            scenarioPass(s) = false;
        end
    end
end

meanDoseCube = sumDose / nScenarios;
stdDoseCube = sqrt(max(sumSqDose / nScenarios - meanDoseCube .^ 2, 0));

robustEval = struct();
robustEval.scenarioMetrics = scenarioMetrics;
robustEval.worstCase = worstCase;
robustEval.passRate = mean(double(scenarioPass));
robustEval.scenarioPass = scenarioPass;
robustEval.nominalScenarioIndex = 1;
robustEval.stdDoseCube = stdDoseCube;
robustEval.meanDoseCube = meanDoseCube;
robustEval.nScenarios = nScenarios;
robustEval.usedStf = ~isempty(stf); %#ok<NASGU>
robustEval.mode = resolve_robust_eval_mode(pln, dij);

end

function fields = get_scenario_dose_fields(dij)
physicalDose = dij.physicalDose;
if iscell(physicalDose)
    fields = physicalDose(:);
else
    fields = {physicalDose};
end
end

function doseDims = infer_dose_dimensions(dij, ct, firstDoseField)
if isfield(dij, 'doseGrid') && isfield(dij.doseGrid, 'dimensions') && ...
        numel(dij.doseGrid.dimensions) == 3
    doseDims = reshape(double(dij.doseGrid.dimensions), 1, []);
    return;
end

if isstruct(ct) && isfield(ct, 'cubeDim') && numel(ct.cubeDim) == 3
    doseDims = reshape(double(ct.cubeDim), 1, []);
    return;
end

if isstruct(ct) && isfield(ct, 'cubeHU') && iscell(ct.cubeHU) && ~isempty(ct.cubeHU) && ~isempty(ct.cubeHU{1})
    doseDims = size(ct.cubeHU{1});
    return;
end

if isvector(firstDoseField)
    doseDims = [numel(firstDoseField), 1, 1];
else
    error('evaluate_robust_scenarios:missingDimensions', ...
        'Could not infer dose dimensions from dij or ct.');
end
end

function doseVector = compute_scenario_dose_vector(doseMatrix, w)
if isvector(doseMatrix)
    doseVector = reshape(double(doseMatrix), [], 1);
    return;
end

doseVector = reshape(full(doseMatrix * reshape(double(w), [], 1)), [], 1);
end

function scenarioInfo = describe_scenario(pln, scenarioIndex)
scenarioInfo = struct('scenarioIndex', scenarioIndex);

if ~isstruct(pln) || ~isfield(pln, 'multScen') || isempty(pln.multScen)
    return;
end

multScen = pln.multScen;
try
    if isobject(multScen) && isprop(multScen, 'linearMask') && scenarioIndex <= size(multScen.linearMask, 1)
        scenarioInfo.linearMask = multScen.linearMask(scenarioIndex, :);
    elseif isstruct(multScen) && isfield(multScen, 'linearMask') && scenarioIndex <= size(multScen.linearMask, 1)
        scenarioInfo.linearMask = multScen.linearMask(scenarioIndex, :);
    end
catch
    % Best-effort metadata only.
end
end

function tolerance = get_robust_tolerance(cfg)
tolerance = 2.0;
if isstruct(cfg) && isfield(cfg, 'evaluation') && isstruct(cfg.evaluation) && ...
        isfield(cfg.evaluation, 'robustD95TolerancePercent') && ...
        isnumeric(cfg.evaluation.robustD95TolerancePercent) && ...
        isscalar(cfg.evaluation.robustD95TolerancePercent) && ...
        isfinite(cfg.evaluation.robustD95TolerancePercent)
    tolerance = double(cfg.evaluation.robustD95TolerancePercent);
end
end

function mode = resolve_robust_eval_mode(pln, dij)
if isstruct(pln) && isfield(pln, 'multScen') && ~isempty(pln.multScen)
    mode = 'robust_or_multiscenario';
elseif isstruct(dij) && isfield(dij, 'physicalDose') && iscell(dij.physicalDose) && numel(dij.physicalDose) > 1
    mode = 'evaluation_only_multiscenario';
else
    mode = 'nominal_only';
end
end
