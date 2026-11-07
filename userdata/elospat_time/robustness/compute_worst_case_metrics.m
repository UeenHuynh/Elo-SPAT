function worstCase = compute_worst_case_metrics(scenarioMetrics)
% compute_worst_case_metrics  Aggregate worst-case summary metrics.
%
% INPUT
%   scenarioMetrics - struct array from evaluate_robust_scenarios
%
% OUTPUT
%   worstCase - struct summarising target/OAR worst cases

narginchk(1, 1);

if ~isstruct(scenarioMetrics) || isempty(scenarioMetrics)
    error('compute_worst_case_metrics:invalidInput', ...
        'scenarioMetrics must be a non-empty struct array.');
end

targetD95 = reshape([scenarioMetrics.targetD95Min], [], 1);
targetD98 = reshape([scenarioMetrics.targetD98Min], [], 1);
targetV95 = reshape([scenarioMetrics.targetV95Min], [], 1);
oarDmax = reshape([scenarioMetrics.oarDmaxMax], [], 1);

[worstTargetD95, ixTargetD95] = min_with_index(targetD95);
[worstTargetD98, ixTargetD98] = min_with_index(targetD98);
[worstTargetV95, ixTargetV95] = min_with_index(targetV95);
[worstOarDmax, ixOarDmax] = max_with_index(oarDmax);

worstCase = struct();
worstCase.targetD95Min = worstTargetD95;
worstCase.targetD98Min = worstTargetD98;
worstCase.targetV95Min = worstTargetV95;
worstCase.oarDmaxMax = worstOarDmax;
worstCase.worstTargetD95Scenario = ixTargetD95;
worstCase.worstTargetD98Scenario = ixTargetD98;
worstCase.worstTargetV95Scenario = ixTargetV95;
worstCase.worstOarDmaxScenario = ixOarDmax;
end

function [value, index] = min_with_index(values)
validMask = isfinite(values);
if ~any(validMask)
    value = NaN;
    index = NaN;
    return;
end

validValues = values(validMask);
validIndices = find(validMask);
[value, localIx] = min(validValues);
index = validIndices(localIx);
end

function [value, index] = max_with_index(values)
validMask = isfinite(values);
if ~any(validMask)
    value = NaN;
    index = NaN;
    return;
end

validValues = values(validMask);
validIndices = find(validMask);
[value, localIx] = max(validValues);
index = validIndices(localIx);
end
