function [score, components] = composite_score(w, groups, cpGroups, cfg, bdtInput)
% composite_score  Compute the wrapper-mode composite score.
%
% score = lambdaEL * R_EL + lambdaB * B + lambdaES * Phi_ES
%         + lambdaS * R_spot + lambdaT * T_estimated
%
% INPUTS
%   w        - spot-weight vector
%   groups   - energy-layer grouping struct array
%   cpGroups - control-point grouping struct array
%   cfg      - configuration struct with lambda* and eps* fields
%   bdtInput - optional scalar estimated BDT or struct with field T_total
%
% OUTPUTS
%   score      - total weighted score
%   components - struct with individual component values

narginchk(4, 5);

validate_cfg_fields(cfg);

if nargin < 5
    bdtInput = [];
end

[l21halfVal, ~] = l21half_group_penalty(w, groups, cfg.epsEL);
[barrierVal, ~] = control_point_barrier(w, cpGroups, cfg.epsB);
[sequencingVal, switchUpCount] = energy_sequencing_penalty(w, groups, cfg.epsZ);
[spotVal, ~] = spot_sparsity_penalty(w, cfg.epsS);
bdtVal = resolve_bdt_value(bdtInput, cfg.lambdaT);

components = struct( ...
    'l21half', l21halfVal, ...
    'barrier', barrierVal, ...
    'sequencing', sequencingVal, ...
    'spotSparsity', spotVal, ...
    'bdt', bdtVal, ...
    'switchUpCount', switchUpCount, ...
    'weightedL21half', cfg.lambdaEL * l21halfVal, ...
    'weightedBarrier', cfg.lambdaB * barrierVal, ...
    'weightedSequencing', cfg.lambdaES * sequencingVal, ...
    'weightedSpotSparsity', cfg.lambdaS * spotVal, ...
    'weightedBdt', cfg.lambdaT * bdtVal);

score = components.weightedL21half + ...
    components.weightedBarrier + ...
    components.weightedSequencing + ...
    components.weightedSpotSparsity + ...
    components.weightedBdt;

end

function validate_cfg_fields(cfg)
requiredFields = { ...
    'lambdaEL', 'lambdaB', 'lambdaES', 'lambdaS', 'lambdaT', ...
    'epsEL', 'epsB', 'epsZ', 'epsS'};

if ~isstruct(cfg)
    error('composite_score:invalidInput', ...
        'cfg must be a struct.');
end

for i = 1:numel(requiredFields)
    fieldName = requiredFields{i};
    if ~isfield(cfg, fieldName)
        error('composite_score:missingField', ...
            'cfg.%s is required.', fieldName);
    end
end
end

function bdtVal = resolve_bdt_value(bdtInput, lambdaT)
if isempty(bdtInput)
    if lambdaT ~= 0
        error('composite_score:missingBdt', ...
            'A BDT input is required when cfg.lambdaT is non-zero.');
    end
    bdtVal = 0;
    return;
end

if isnumeric(bdtInput) && isscalar(bdtInput) && isfinite(bdtInput)
    bdtVal = double(bdtInput);
    return;
end

if isstruct(bdtInput) && isfield(bdtInput, 'T_total') && ...
        isscalar(bdtInput.T_total) && isnumeric(bdtInput.T_total) && ...
        isfinite(bdtInput.T_total)
    bdtVal = double(bdtInput.T_total);
    return;
end

error('composite_score:invalidBdt', ...
    'bdtInput must be empty, a finite scalar, or a struct with field T_total.');
end
