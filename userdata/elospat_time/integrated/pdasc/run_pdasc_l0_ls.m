function result = run_pdasc_l0_ls(U, p0, settings)
%RUN_PDASC_L0_LS Published-form least-squares/L0 PDASC continuation.
%
% Solves the Zhao et al. (2023) support-selection model
%   min_w 0.5*||U*w-p0||_2^2 + lambda*||w||_0
% by geometric continuation and primal-dual active sets. Non-negativity is
% applied as the final projection reported in that paper; it is not inserted
% into the restricted least-squares update.

% Reference: Zhao et al., Acta Oncol. 2023;62:48-52,
% doi:10.1080/0284186X.2023.2172689, Algorithm 1.

narginchk(2, 3);
if nargin < 3 || isempty(settings), settings = struct(); end
settings = defaults(settings);
validateattributes(U, {'numeric'}, {'2d','real','finite'});
p0 = reshape(double(p0), [], 1);
if size(U, 1) ~= numel(p0) || isempty(U) || any(~isfinite(p0))
    error('run_pdasc_l0_ls:invalidInput', ...
        'U and p0 must be finite and have matching row dimensions.');
end
U = sparse(double(U));
n = size(U, 2);

% Positive diagonal column normalization is an equivalent change of
% variables x_i=||U_i||_2*w_i: it preserves residuals, L0 support, and signs.
columnNorm = full(sqrt(sum(U.^2, 1)))';
nonzeroColumn = columnNorm > settings.zeroColumnTolerance;
scale = columnNorm;
scale(~nonzeroColumn) = 1;
if settings.normalizeColumns
    Phi = U * spdiags(1 ./ scale, 0, n, n);
else
    Phi = U;
    scale(:) = 1;
end

x = zeros(n, 1);
d = reshape(double(Phi' * p0), [], 1);
lambdaMax = 0.5 * norm(d, Inf)^2;
if lambdaMax == 0
    result = empty_result(x, scale, nonzeroColumn, settings, lambdaMax);
    return;
end
if isempty(settings.lambdaMin)
    lambdaMin = lambdaMax * settings.lambdaMinRatio;
else
    lambdaMin = settings.lambdaMin;
end
if ~(isfinite(lambdaMin) && lambdaMin > 0 && lambdaMin <= lambdaMax)
    error('run_pdasc_l0_ls:invalidLambda', ...
        'lambdaMin must lie in (0, lambdaMax].');
end
lambdaGrid = exp(linspace(log(lambdaMax), log(lambdaMin), ...
    settings.maxContinuationLevels))';
lambdaGrid(1) = lambdaMax;
lambdaGrid(end) = lambdaMin;

history = repmat(empty_history(), settings.maxContinuationLevels, 1);
states = cell(settings.maxContinuationLevels, 1);
previousOuterSupport = false(n, 1);
timer = tic;
stopReason = 'maximum_continuation_levels';
usedLevels = settings.maxContinuationLevels;

for k = 1:settings.maxContinuationLevels
    lambda = lambdaGrid(k);
    threshold = settings.gamma * sqrt(2 * lambda);
    activePrevious = previousOuterSupport;
    innerStop = 'maximum_inner_iterations';
    innerIterations = 0;
    restrictedSolveSeconds = 0;
    restrictedResidual = NaN;
    gramConditionEstimate = NaN;

    for j = 1:settings.maxInnerIterations
        active = abs(x + settings.gamma * d) > threshold;
        active = active & nonzeroColumn;
        if isequal(active, activePrevious)
            innerStop = 'active_set_stable';
            break;
        end
        innerIterations = j;
        activePrevious = active;
        x(:) = 0;
        ix = find(active);
        solveTimer = tic;
        if ~isempty(ix)
            PhiA = Phi(:, ix);
            switch settings.restrictedSolver
                case 'normal_equations'
                    gram = PhiA' * PhiA;
                    rhs = PhiA' * p0;
                    if numel(ix) <= settings.maximumConditionEstimateSize
                        try, gramConditionEstimate = condest(gram); catch, end
                    end
                    x(ix) = gram \ rhs;
                case 'lsqr'
                    [x(ix), flag, relres] = lsqr(PhiA, p0, ...
                        settings.lsqrTolerance, settings.lsqrMaxIterations);
                    if flag > 1
                        warning('run_pdasc_l0_ls:lsqr', ...
                            'LSQR flag %d at continuation %d, inner %d.', flag, k, j);
                    end
                    restrictedResidual = relres;
                otherwise
                    error('run_pdasc_l0_ls:invalidSolver', ...
                        'Unknown restricted solver %s.', settings.restrictedSolver);
            end
        end
        restrictedSolveSeconds = restrictedSolveSeconds + toc(solveTimer);
        residual = p0 - Phi * x;
        d = reshape(double(Phi' * residual), [], 1);
    end

    rawActive = activePrevious;
    wRaw = x ./ scale;
    wProjected = max(wRaw, 0);
    projectedSupport = wProjected > settings.supportTolerance;
    residualRaw = U * wRaw - p0;
    residualProjected = U * wProjected - p0;
    rawObjective = 0.5 * sum(residualRaw.^2) + lambda * nnz(rawActive);
    projectedObjective = 0.5 * sum(residualProjected.^2) + ...
        lambda * nnz(projectedSupport);
    additions = nnz(rawActive & ~previousOuterSupport);
    removals = nnz(~rawActive & previousOuterSupport);

    h = empty_history();
    h.level = k; h.lambda = lambda; h.threshold = threshold;
    h.innerIterations = innerIterations; h.innerStop = innerStop;
    h.rawActiveSpots = nnz(rawActive);
    h.projectedActiveSpots = nnz(projectedSupport);
    h.negativeActiveCoefficients = nnz(rawActive & wRaw < 0);
    h.newActiveSpots = additions; h.removedActiveSpots = removals;
    h.rawLeastSquares = 0.5 * sum(residualRaw.^2);
    h.projectedLeastSquares = 0.5 * sum(residualProjected.^2);
    h.rawObjective = rawObjective; h.projectedObjective = projectedObjective;
    h.dualInfinityNorm = norm(d, Inf);
    h.restrictedResidual = restrictedResidual;
    h.gramConditionEstimate = gramConditionEstimate;
    h.restrictedSolveSeconds = restrictedSolveSeconds;
    h.elapsedSeconds = toc(timer);
    h.rawSupportHash = support_hash(rawActive);
    h.projectedSupportHash = support_hash(projectedSupport);
    history(k) = h;
    states{k} = struct('level', k, 'lambda', lambda, 'x', x, ...
        'd', d, 'wRaw', wRaw, 'w', wProjected, ...
        'rawSupport', rawActive, 'support', projectedSupport);
    previousOuterSupport = rawActive;

    if settings.verbose
        fprintf(['PDASC k=%d lambda=%.6e inner=%d raw/projected=%d/%d ' ...
            'LS=%.6e stop=%s\n'], k, lambda, innerIterations, ...
            nnz(rawActive), nnz(projectedSupport), ...
            h.projectedLeastSquares, innerStop);
    end
    if nnz(projectedSupport) >= settings.targetActiveSpots
        stopReason = 'target_active_spots_reached';
        usedLevels = k;
        break;
    end
end

history = history(1:usedLevels);
states = states(1:usedLevels);
last = states{end};
result = struct();
result.w = last.w; result.wRaw = last.wRaw; result.x = last.x;
result.d = last.d; result.support = last.support;
result.rawSupport = last.rawSupport; result.history = history;
result.states = states; result.lambdaGrid = lambdaGrid(1:usedLevels);
result.lambdaMax = lambdaMax; result.lambdaMin = lambdaMin;
result.columnNorm = columnNorm; result.changeOfVariablesScale = scale;
result.nonzeroColumn = nonzeroColumn; result.settings = settings;
result.levels = usedLevels; result.stopReason = stopReason;
result.wallSeconds = toc(timer);
result.formulation = 'published least-squares/L0 PDASC support selector';
end

function settings = defaults(override)
settings = struct('gamma', 1, 'lambdaMin', [], 'lambdaMinRatio', 1e-10, ...
    'maxContinuationLevels', 41, 'maxInnerIterations', 50, ...
    'targetActiveSpots', Inf, 'normalizeColumns', true, ...
    'zeroColumnTolerance', 0, 'supportTolerance', 0, ...
    'restrictedSolver', 'normal_equations', ...
    'maximumConditionEstimateSize', 2500, ...
    'lsqrTolerance', 1e-10, 'lsqrMaxIterations', 2000, ...
    'verbose', true);
names = fieldnames(override);
for i = 1:numel(names), settings.(names{i}) = override.(names{i}); end
validateattributes(settings.gamma, {'numeric'}, {'scalar','positive','finite'});
validateattributes(settings.lambdaMinRatio, {'numeric'}, ...
    {'scalar','positive','<=',1,'finite'});
validateattributes(settings.maxContinuationLevels, {'numeric'}, ...
    {'scalar','integer','>=',2});
validateattributes(settings.maxInnerIterations, {'numeric'}, ...
    {'scalar','integer','positive'});
settings.restrictedSolver = lower(char(settings.restrictedSolver));
end

function h = empty_history()
h = struct('level', NaN, 'lambda', NaN, 'threshold', NaN, ...
    'innerIterations', NaN, 'innerStop', '', 'rawActiveSpots', NaN, ...
    'projectedActiveSpots', NaN, 'negativeActiveCoefficients', NaN, ...
    'newActiveSpots', NaN, 'removedActiveSpots', NaN, ...
    'rawLeastSquares', NaN, 'projectedLeastSquares', NaN, ...
    'rawObjective', NaN, 'projectedObjective', NaN, ...
    'dualInfinityNorm', NaN, 'restrictedResidual', NaN, ...
    'gramConditionEstimate', NaN, 'restrictedSolveSeconds', NaN, ...
    'elapsedSeconds', NaN, 'rawSupportHash', '', 'projectedSupportHash', '');
end

function result = empty_result(x, scale, nonzeroColumn, settings, lambdaMax)
result = struct('w', x, 'wRaw', x, 'x', x, 'd', x, ...
    'support', false(size(x)), 'rawSupport', false(size(x)), ...
    'history', repmat(empty_history(), 0, 1), ...
    'lambdaGrid', zeros(0,1), 'lambdaMax', lambdaMax, 'lambdaMin', 0, ...
    'columnNorm', scale, 'changeOfVariablesScale', scale, ...
    'nonzeroColumn', nonzeroColumn, 'settings', settings, 'levels', 0, ...
    'stopReason', 'zero_correlation', 'wallSeconds', 0, ...
    'formulation', 'published least-squares/L0 PDASC support selector');
result.states = cell(0,1);
end

function value = support_hash(mask)
if exist('OCTAVE_VERSION', 'builtin')
    value = lower(hash('sha256', char(uint8(mask(:))' + uint8('0'))));
else
    md = java.security.MessageDigest.getInstance('SHA-256');
    md.update(uint8(mask(:))); bytes = typecast(md.digest(), 'uint8');
    value = lower(reshape(dec2hex(bytes, 2).', 1, []));
end
end
