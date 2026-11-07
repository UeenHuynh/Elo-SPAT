function result = solve_admm_w_first_order(doseModel, w0, center, rho, settings)
%SOLVE_ADMM_W_FIRST_ORDER Limited projected-gradient ADMM w-update.

narginchk(4, 5);
if nargin < 5, settings = struct(); end
settings = defaults(settings, struct('maxIterations', 5, ...
    'initialStep', 1e-3, 'backtrackFactor', 0.5, 'stepGrowth', 1.25, ...
    'armijo', 1e-4, 'minimumStep', 1e-14, 'maximumStep', 1, ...
    'projectedGradientTolerance', 1e-6));
w = max(0, reshape(double(w0), [], 1));
center = reshape(double(center), [], 1);
history = nan(settings.maxIterations, 5);
counts = struct('objective', 0, 'gradient', 0);
step = settings.initialStep;
timer = tic;
reason = 'maximum_iterations';

for k = 1:settings.maxIterations
    [f, g] = fg(w); %#ok<ASGLU>
    pg = w - max(0, w - g);
    pgNorm = norm(pg) / max(1, norm(w));
    if pgNorm <= settings.projectedGradientTolerance
        reason = 'projected_stationarity';
        history(k, :) = [k, f, pgNorm, step, 0];
        break;
    end
    trialStep = min(settings.maximumStep, step);
    accepted = false;
    trials = 0;
    while trialStep >= settings.minimumStep
        candidate = max(0, w - trialStep * g);
        d = candidate - w;
        fCandidate = f_only(candidate);
        trials = trials + 1;
        if fCandidate <= f + settings.armijo * (g' * d)
            accepted = true;
            break;
        end
        trialStep = trialStep * settings.backtrackFactor;
    end
    if ~accepted
        reason = 'line_search_failed';
        history(k, :) = [k, f, pgNorm, trialStep, trials];
        break;
    end
    w = candidate;
    step = min(settings.maximumStep, trialStep * settings.stepGrowth);
    history(k, :) = [k, fCandidate, pgNorm, trialStep, trials];
end
last = find(isfinite(history(:, 1)), 1, 'last');
if isempty(last), last = 0; end
history = history(1:last, :);
[finalObjective, finalGradient, diagnostics] = fg(w);
result = struct('w', w, 'objective', finalObjective, ...
    'gradient', finalGradient, 'diagnostics', diagnostics, ...
    'iterations', last, 'reason', reason, 'history', history, ...
    'callbackCounts', counts, 'wallSeconds', toc(timer), ...
    'strategy', 'limited_projected_gradient', 'settings', settings);

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

function out = defaults(in, base)
out = base;
names = fieldnames(in);
for i = 1:numel(names), out.(names{i}) = in.(names{i}); end
end
