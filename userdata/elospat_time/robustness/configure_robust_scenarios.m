function pln = configure_robust_scenarios(pln, ct, cfg)
% configure_robust_scenarios  Attach robust scenario metadata to pln.multScen.
%
% INPUTS
%   pln - matRad plan struct
%   ct  - optional ct struct (used only if a scenario object is instantiated)
%   cfg - configuration struct with cfg.robust.*
%
% OUTPUT
%   pln - plan with pln.multScen configured for nominal or worst-case use

narginchk(2, 3);

if nargin < 3
    cfg = struct();
end

if ~isstruct(pln)
    error('configure_robust_scenarios:invalidInput', ...
        'pln must be a struct.');
end

pln = normalize_plan_machine_fields(pln);

if ~isstruct(cfg) || ~isfield(cfg, 'robust') || ~isstruct(cfg.robust)
    error('configure_robust_scenarios:invalidInput', ...
        'cfg.robust must be provided.');
end

robustCfg = cfg.robust;
enabled = isfield(robustCfg, 'enabled') && logical(robustCfg.enabled);

if ~enabled
    pln.multScen = 'nomScen';
    return;
end

setupMM = require_numeric_row(robustCfg, 'setupMM', 3);
rangePercent = require_numeric_scalar(robustCfg, 'rangePercent');

scenarioMetadata = build_scenario_metadata(setupMM, rangePercent);
pln.multScen = build_scenario_model(pln, ct, scenarioMetadata);

end

function multScen = build_scenario_model(pln, ct, scenarioMetadata)
try
    % matRad's class files are discoverable through its dynamic class path,
    % but exist(...,'class') can return 0 before the class has been loaded.
    % Construct directly and only fall back to metadata if construction is
    % genuinely unavailable.
    if isfield(pln, 'multScen') && isa(pln.multScen, 'matRad_WorstCaseScenarios')
        multScen = pln.multScen;
    elseif nargin >= 2 && ~isempty(ct)
        multScen = matRad_WorstCaseScenarios(ct);
    else
        multScen = matRad_WorstCaseScenarios();
    end

    multScen.shiftSD = sanitize_positive_shift_sd(scenarioMetadata.shiftSD);
    multScen.rangeRelSD = scenarioMetadata.rangeRelSD;
    % matRad 3.2.2 cannot update its scenario grid with rangeAbsSD=0
    % (the internal covariance normalization divides by a zero-width grid).
    % Keep the physical absolute range uncertainty at the documented
    % 1 mm default unless a positive value is supplied.
    if scenarioMetadata.rangeAbsSD > 0
        multScen.rangeAbsSD = scenarioMetadata.rangeAbsSD;
    end

    if isprop(multScen, 'wcSigma')
        multScen.wcSigma = scenarioMetadata.wcSigma;
    end
    if isprop(multScen, 'combinations')
        multScen.combinations = scenarioMetadata.combinations;
    end
    if isprop(multScen, 'combineRange')
        multScen.combineRange = scenarioMetadata.combineRange;
    end
catch
    multScen = scenarioMetadata;
end
end

function scenarioMetadata = build_scenario_metadata(setupMM, rangePercent)
scenarioMetadata = struct();
scenarioMetadata.model = 'wcScen';
scenarioMetadata.wcSigma = 1;
scenarioMetadata.shiftSD = abs(setupMM(:))';
scenarioMetadata.rangeAbsSD = 1;
scenarioMetadata.rangeRelSD = abs(rangePercent);
scenarioMetadata.combinations = 'none';
scenarioMetadata.combineRange = true;
end

function shiftSD = sanitize_positive_shift_sd(shiftSD)
shiftSD = reshape(double(shiftSD), 1, []);
shiftSD(shiftSD <= 0) = eps;
end

function value = require_numeric_row(s, fieldName, expectedLength)
if ~isfield(s, fieldName) || ~isnumeric(s.(fieldName)) || numel(s.(fieldName)) ~= expectedLength || ...
        any(~isfinite(s.(fieldName)))
    error('configure_robust_scenarios:invalidInput', ...
        'cfg.robust.%s must be a finite %d-element numeric vector.', ...
        fieldName, expectedLength);
end

value = reshape(double(s.(fieldName)), 1, []);
end

function value = require_numeric_scalar(s, fieldName)
if ~isfield(s, fieldName) || ~isnumeric(s.(fieldName)) || ~isscalar(s.(fieldName)) || ...
        ~isfinite(s.(fieldName))
    error('configure_robust_scenarios:invalidInput', ...
        'cfg.robust.%s must be a finite scalar.', fieldName);
end

value = double(s.(fieldName));
end
