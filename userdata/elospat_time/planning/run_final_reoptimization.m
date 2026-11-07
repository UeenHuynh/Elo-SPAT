function result = run_final_reoptimization(dij, cst, pln, wInit, cfg, supportMask, layerMask)
% run_final_reoptimization  Final nominal re-optimization stage.

narginchk(5, 7);

nFull = double(dij.totalNumOfBixels);
wInit = reshape(double(wInit), [], 1);
if numel(wInit) ~= nFull
    error('run_final_reoptimization:dimensionMismatch', ...
        'wInit length %d does not match full dij length %d.', numel(wInit), nFull);
end

% Legacy callers remain valid, but the orchestrator must pass an explicit
% support mask.  An omitted mask means no support restriction.
if nargin < 6 || isempty(supportMask)
    supportMask = true(nFull, 1);
else
    supportMask = reshape(logical(supportMask), [], 1);
    if numel(supportMask) ~= nFull
        error('run_final_reoptimization:dimensionMismatch', ...
            'supportMask length %d does not match full dij length %d.', ...
            numel(supportMask), nFull);
    end
end
if nargin < 7
    layerMask = [];
end
if any(abs(wInit(~supportMask)) > 0)
    error('run_final_reoptimization:supportLeakage', ...
        'Initial weights are nonzero outside the supplied support.');
end

if ~isstruct(cfg) || ~isfield(cfg, 'robust') || ~isstruct(cfg.robust)
    error('run_final_reoptimization:invalidInput', ...
        'cfg.robust must be provided.');
end

plnConfigured = normalize_plan_machine_fields(pln);
cstConfigured = cst;

if isfield(cfg.robust, 'enabled') && logical(cfg.robust.enabled)
    result.robustMode = char(string(cfg.robust.mode));
    result.robustDeferred = true;
    % Keep the final warm-start re-optimization nominal. The dense-arc dij
    % passed into this stage contains only the nominal scenario, so forcing a
    % multi-scenario pln here causes matRad to index past dij{1}. Robust
    % scenario expansion is handled later by run_task28_30_from_nominal()
    % after the nominal artifact has been saved.
else
    result.robustMode = 'none';
    result.robustDeferred = false;
end

[dijReduced, reducedToFullIdx, fullToReducedIdx] = reduce_dij_to_support(dij, supportMask);
wReduced = wInit(reducedToFullIdx);

% The integrated path may request an exact positive lower bound on every
% selected variable. This makes the selected support and minimum-weight
% definition consistent during final IPOPT reoptimization. Legacy wrapper
% callers do not set this opt-in field and therefore remain unchanged.
if isfield(cfg, 'enforceSelectedMinimumWeight') && ...
        logical(cfg.enforceSelectedMinimumWeight) && cfg.mmuThreshold > 0
    nReduced = double(dijReduced.totalNumOfBixels);
    if isfield(dijReduced, 'numParticlesPerMU') && ...
            ~isempty(dijReduced.numParticlesPerMU)
        particlesPerMU = reshape(double(dijReduced.numParticlesPerMU), [], 1);
        if isscalar(particlesPerMU)
            particlesPerMU = repmat(particlesPerMU, nReduced, 1);
        end
    else
        % Weight-space identity calibration: matRad computes
        % minimumW = numParticlesPerMU .* minMU / 1e6.
        particlesPerMU = 1e6 * ones(nReduced, 1);
        dijReduced.numParticlesPerMU = particlesPerMU;
    end
    if numel(particlesPerMU) ~= nReduced || any(particlesPerMU <= 0)
        error('run_final_reoptimization:invalidMUCalibration', ...
            'numParticlesPerMU must be positive and match reduced support.');
    end
    dijReduced.minMU = cfg.mmuThreshold * 1e6 ./ particlesPerMU;
    if ~isfield(plnConfigured, 'propOpt') || isempty(plnConfigured.propOpt)
        plnConfigured.propOpt = struct();
    end
    plnConfigured.propOpt.boundMU = true;
end
[resultGuiReduced, optimizer] = reoptimize_reduced_problem( ...
    dijReduced, cstConfigured, plnConfigured, wReduced);

wFinal = zeros(nFull, 1);
wFinal(reducedToFullIdx) = reshape(double(resultGuiReduced.w), [], 1);
tol = max(1e-10, 1e-8 * max(1, max(abs(wFinal))));
if any(abs(wFinal(~supportMask)) > tol) || nnz(abs(wFinal) > tol) > nnz(supportMask)
    error('run_final_reoptimization:supportLeakage', ...
        'Final optimization activated weights outside the allowed support.');
end
if isfield(cfg, 'enforceSelectedMinimumWeight') && ...
        logical(cfg.enforceSelectedMinimumWeight) && cfg.mmuThreshold > 0 && ...
        any(wFinal(supportMask) < cfg.mmuThreshold - tol)
    error('run_final_reoptimization:minimumWeightViolation', ...
        'Final optimization violated the selected-spot minimum weight.');
end
result = merge_result(result, resultGuiReduced, optimizer, plnConfigured, cstConfigured);
result.resultGUI.w = wFinal;
result.w = wFinal;
result.supportMask = supportMask;
result.reducedToFullIdx = reducedToFullIdx;
result.fullToReducedIdx = fullToReducedIdx;
result.layerMask = layerMask;
end

function result = merge_result(result, resultGui, optimizer, plnConfigured, cstConfigured)
result.w = resultGui.w;
result.resultGUI = resultGui;
result.optimizer = optimizer;
result.pln = plnConfigured;
result.cst = cstConfigured;
end
