function cstOut = configure_robust_objectives(cst, cfg)
% configure_robust_objectives  Set robustness mode on all cst objectives.
%
% INPUTS
%   cst - matRad cst cell array
%   cfg - configuration struct with cfg.robust.enabled / cfg.robust.mode
%
% OUTPUT
%   cstOut - cst with objective robustness modes updated

narginchk(2, 2);

if ~iscell(cst) || size(cst, 2) < 6
    error('configure_robust_objectives:invalidInput', ...
        'cst must be a matRad cst cell array with at least 6 columns.');
end

if ~isstruct(cfg) || ~isfield(cfg, 'robust') || ~isstruct(cfg.robust)
    error('configure_robust_objectives:invalidInput', ...
        'cfg.robust must be provided.');
end

if isfield(cfg.robust, 'enabled') && logical(cfg.robust.enabled)
    robustnessMode = char(string(cfg.robust.mode));
else
    robustnessMode = 'none';
end

cstOut = cst;
for i = 1:size(cstOut, 1)
    objectives = cstOut{i, 6};
    if isempty(objectives)
        continue;
    end
    for j = 1:numel(objectives)
        objectives{j} = set_objective_robustness(objectives{j}, robustnessMode);
    end
    cstOut{i, 6} = objectives;
end

end

function objective = set_objective_robustness(objective, robustnessMode)
if isstruct(objective)
    objective.robustness = robustnessMode;
elseif isobject(objective) && isprop(objective, 'robustness')
    objective.robustness = robustnessMode;
end
end
