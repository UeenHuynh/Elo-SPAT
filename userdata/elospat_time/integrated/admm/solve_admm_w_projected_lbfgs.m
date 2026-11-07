function result = solve_admm_w_projected_lbfgs(doseModel, w0, center, rho, settings)
%SOLVE_ADMM_W_PROJECTED_LBFGS Limited projected L-BFGS ADMM w-update.
%
% This is an inexact bound-constrained quasi-Newton benchmark, not a claim
% of equivalence to L-BFGS-B. A non-descent projected direction falls back
% to projected gradient.

narginchk(4, 5);
if nargin < 5, settings = struct(); end
settings = defaults(settings, struct('maxIterations', 5, 'memory', 5, ...
    'initialStep', 1, 'backtrackFactor', 0.5, 'armijo', 1e-4, ...
    'minimumStep', 1e-14, 'curvatureTolerance', 1e-12, ...
    'projectedGradientTolerance', 1e-6));
w = max(0, reshape(double(w0), [], 1));
center = reshape(double(center), [], 1);
S = zeros(numel(w), 0); Y = S; invCurvature = zeros(0, 1);
history = nan(settings.maxIterations, 6);
counts = struct('objective', 0, 'gradient', 0);
timer = tic; reason = 'maximum_iterations';
[f, g, diagnostics] = fg(w);

for k = 1:settings.maxIterations
    pg = w - max(0, w - g);
    pgNorm = norm(pg) / max(1, norm(w));
    if pgNorm <= settings.projectedGradientTolerance
        reason = 'projected_stationarity';
        history(k, :) = [k, f, pgNorm, 0, 0, size(S, 2)];
        break;
    end
    direction = -two_loop(g, S, Y, invCurvature);
    candidateDirection = max(0, w + direction) - w;
    if g' * candidateDirection >= -1e-12 * norm(g) * max(1, norm(candidateDirection))
        direction = -g;
    end
    step = settings.initialStep; accepted = false; trials = 0;
    while step >= settings.minimumStep
        candidate = max(0, w + step * direction);
        d = candidate - w;
        fCandidate = f_only(candidate);
        trials = trials + 1;
        if fCandidate <= f + settings.armijo * (g' * d)
            accepted = true;
            break;
        end
        step = step * settings.backtrackFactor;
    end
    if ~accepted
        reason = 'line_search_failed';
        history(k, :) = [k, f, pgNorm, step, trials, size(S, 2)];
        break;
    end
    [~, gNew, diagnostics] = fg(candidate);
    s = candidate - w; y = gNew - g; sy = s' * y;
    if sy > settings.curvatureTolerance * norm(s) * max(1, norm(y))
        if size(S, 2) == settings.memory
            S(:, 1) = []; Y(:, 1) = []; invCurvature(1) = [];
        end
        S(:, end+1) = s; %#ok<AGROW>
        Y(:, end+1) = y; %#ok<AGROW>
        invCurvature(end+1, 1) = 1 / sy; %#ok<AGROW>
    end
    w = candidate; f = fCandidate; g = gNew;
    history(k, :) = [k, f, pgNorm, step, trials, size(S, 2)];
end
last = find(isfinite(history(:, 1)), 1, 'last');
if isempty(last), last = 0; end
history = history(1:last, :);
result = struct('w', w, 'objective', f, 'gradient', g, ...
    'diagnostics', diagnostics, 'iterations', last, 'reason', reason, ...
    'history', history, 'callbackCounts', counts, ...
    'wallSeconds', toc(timer), 'strategy', 'limited_projected_lbfgs', ...
    'settings', settings);

    function [value, gradient, diag] = fg(x)
        counts.objective = counts.objective + 1;
        counts.gradient = counts.gradient + 1;
        [value, gradient, diag] = evaluate_admm_w_subproblem( ...
            doseModel, x, center, rho, true);
    end
    function value = f_only(x)
        counts.objective = counts.objective + 1;
        value = evaluate_admm_w_subproblem(doseModel, x, center, rho, false);
    end
end

function r = two_loop(g, S, Y, invCurvature)
q = g; m = size(S, 2); alpha = zeros(m, 1);
for j = m:-1:1
    alpha(j) = invCurvature(j) * (S(:, j)' * q);
    q = q - alpha(j) * Y(:, j);
end
if m > 0
    gamma = (S(:, end)' * Y(:, end)) / max(Y(:, end)' * Y(:, end), eps);
else
    gamma = 1;
end
r = gamma * q;
for j = 1:m
    beta = invCurvature(j) * (Y(:, j)' * r);
    r = r + S(:, j) * (alpha(j) - beta);
end
end

function out = defaults(in, base)
out = base;
names = fieldnames(in);
for i = 1:numel(names), out.(names{i}) = in.(names{i}); end
end
