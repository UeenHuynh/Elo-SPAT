function report = evaluate_numerical_stabilization(history, bdtHistory, rule)
% evaluate_numerical_stabilization  Multi-criterion convergence assessment.
%
% history must contain iteration, objective, relativeObjective, stepSize,
% normalizedStationarity and layerSupportChanged. bdtHistory contains
% iteration and BDT values evaluated at saved checkpoints.

narginchk(2, 3);
if nargin < 3 || isempty(rule), rule = default_stabilization_rule(); end
required = {'iteration','objective','relativeObjective','stepSize', ...
    'normalizedStationarity','layerSupportChanged'};
for i = 1:numel(required)
    if ~isfield(history, required{i})
        error('evaluate_numerical_stabilization:missingHistory', ...
            'History is missing %s.', required{i});
    end
end

n = numel(history.iteration);
w = rule.trailingWindow;
report = empty_report(rule);
report.iterations = iff(n == 0, 0, history.iteration(end));
report.enoughIterations = n >= w && report.iterations >= rule.minimumIterations;
if n < w
    report.reason = 'insufficient_iteration_history';
    report = stationarity_label(report, history, rule);
    return;
end

ix = (n-w+1):n;
objective = double(history.objective(ix));
relativeObjective = double(history.relativeObjective(ix));
steps = double(history.stepSize(ix));
report.relativeObjectiveWindow = (objective(1)-objective(end)) / ...
    max(1, abs(objective(1)));
report.maxPerIterationRelativeObjective = max(relativeObjective);
report.maximumObjectiveIncrease = max([0; diff(objective)]);
report.objectiveStable = report.relativeObjectiveWindow <= ...
    rule.relativeObjectiveWindowTolerance && ...
    report.maxPerIterationRelativeObjective <= ...
    rule.maxPerIterationRelativeObjectiveTolerance && ...
    report.maximumObjectiveIncrease <= rule.objectiveIncreaseTolerance * ...
        max(1, max(abs(objective)));
report.supportStable = ~any(logical(history.layerSupportChanged(ix)));
if isfield(history,'spotSupportChanged')
    report.spotSupportEvidenceAvailable = true;
    report.spotSupportStable = ~any(logical(history.spotSupportChanged(ix)));
else
    report.spotSupportEvidenceAvailable = false;
    report.spotSupportStable = true;
end
positiveSteps = steps(isfinite(steps) & steps > 0);
if numel(positiveSteps) == w
    report.stepRatio = max(positiveSteps) / min(positiveSteps);
    report.stepStable = report.stepRatio <= rule.maximumStepRatio;
end

[report.bdtEvidenceAvailable, report.bdtRelativeRange, report.bdtStable] = ...
    evaluate_bdt(bdtHistory, history.iteration(end), rule);
report.numerically_stabilized = report.enoughIterations && report.objectiveStable && ...
    report.supportStable && report.spotSupportStable && ...
    report.stepStable && report.bdtStable;
if report.numerically_stabilized
    report.reason = 'objective_support_step_and_bdt_stabilized';
elseif ~report.enoughIterations
    report.reason = 'insufficient_iteration_history';
elseif ~report.objectiveStable
    report.reason = 'objective_not_stabilized';
elseif ~report.supportStable
    report.reason = 'layer_support_not_stabilized';
elseif ~report.spotSupportStable
    report.reason = 'spot_support_not_stabilized';
elseif ~report.stepStable
    report.reason = 'accepted_step_not_stabilized';
elseif ~report.bdtEvidenceAvailable
    report.reason = 'insufficient_bdt_history';
else
    report.reason = 'bdt_not_stabilized';
end
report = stationarity_label(report, history, rule);
end

function [available, relativeRange, stable] = evaluate_bdt(bdt, endpoint, rule)
available = false; relativeRange = NaN; stable = false;
if isempty(bdt) || ~isstruct(bdt) || ...
        ~isfield(bdt, 'iteration') || ~isfield(bdt, 'value')
    return;
end
valid = isfinite(bdt.iteration) & isfinite(bdt.value) & ...
    bdt.iteration <= endpoint;
iterations = double(bdt.iteration(valid));
values = double(bdt.value(valid));
[iterations, order] = sort(iterations); values = values(order); %#ok<ASGLU>
q = rule.minimumBDTCheckpoints;
if numel(values) < q, return; end
values = values(end-q+1:end);
available = true;
relativeRange = (max(values)-min(values)) / max(1, abs(values(end)));
stable = relativeRange <= rule.relativeBDTRangeTolerance;
end

function report = stationarity_label(report, history, rule)
report.finalNormalizedStationarity = NaN;
if ~isempty(history.normalizedStationarity)
    report.finalNormalizedStationarity = ...
        double(history.normalizedStationarity(end));
end
tol = rule.calibratedStationarityTolerance;
report.stationarityToleranceCalibrated = isfinite(tol) && tol > 0;
report.stationary_converged = report.stationarityToleranceCalibrated && ...
    isfinite(report.finalNormalizedStationarity) && ...
    report.finalNormalizedStationarity <= tol;
end

function report = empty_report(rule)
report = struct('iterations',0,'enoughIterations',false, ...
    'relativeObjectiveWindow',NaN, ...
    'maxPerIterationRelativeObjective',NaN, ...
    'maximumObjectiveIncrease',NaN,'objectiveStable',false, ...
    'supportStable',false,'spotSupportEvidenceAvailable',false, ...
    'spotSupportStable',true,'stepRatio',NaN,'stepStable',false, ...
    'bdtEvidenceAvailable',false,'bdtRelativeRange',NaN, ...
    'bdtStable',false,'numerically_stabilized',false, ...
    'finalNormalizedStationarity',NaN, ...
    'stationarityToleranceCalibrated',false, ...
    'stationary_converged',false,'reason','not_evaluated', ...
    'rule',rule);
end

function value = iff(condition, yes, no)
if condition, value = yes; else, value = no; end
end
