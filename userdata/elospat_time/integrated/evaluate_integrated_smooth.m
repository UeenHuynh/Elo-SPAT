function [value, gradient, parts] = evaluate_integrated_smooth(w, doseModel, context, cfg, needGradient)
% evaluate_integrated_smooth  Smooth part of the integrated objective.

if nargin < 5
    needGradient = nargout >= 2;
end
[doseValue, doseGradient, doseParts] = evaluate_matrad_dose_model( ...
    doseModel, w, needGradient);
value = doseValue;
if needGradient
    gradient = doseGradient;
else
    gradient = [];
end

phiES = 0;
timeValue = 0;
timeParts = struct();
if cfg.lambdaES > 0 || cfg.lambdaT > 0
    [phiES, gradES, timeValue, gradTime, timeParts] = ...
        evaluate_delivery_surrogates(w, context);
    value = value + cfg.lambdaES * phiES + cfg.lambdaT * timeValue;
    if needGradient
        gradient = gradient + cfg.lambdaES * gradES + cfg.lambdaT * gradTime;
    end
end

barrierValue = 0;
if cfg.lambdaB > 0
    for b = 1:numel(context.cpGroups)
        ix = context.cpGroups(b).ix;
        q = sum(w(ix) .^ 2) + cfg.epsB;
        barrierValue = barrierValue - 0.5 * log(q);
        if needGradient
            gradient(ix) = gradient(ix) - cfg.lambdaB * w(ix) / q;
        end
    end
    value = value + cfg.lambdaB * barrierValue;
end
parts = struct('dose', doseParts.doseObjective, ...
    'constraintMerit', doseParts.constraintMerit, ...
    'maxConstraintViolation', doseParts.maxConstraintViolation, ...
    'energySequencing', phiES, 'timeSurrogate', timeValue, ...
    'timeParts', timeParts, 'controlPointBarrier', barrierValue);
end
