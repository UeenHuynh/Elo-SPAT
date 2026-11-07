function [value, gradient, diagnostics] = evaluate_matrad_dose_model(model, w, needGradient)
% evaluate_matrad_dose_model  Evaluate matRad objective and constraint merit.

narginchk(2, 3);
if nargin < 3
    needGradient = nargout >= 2;
end
w = reshape(double(w), [], 1);
if numel(w) ~= model.nWeights || any(~isfinite(w))
    error('evaluate_matrad_dose_model:invalidWeights', ...
        'w must be a finite vector with %d entries.', model.nWeights);
end

if isfield(model, 'objectiveFunction')
    doseValue = model.objectiveFunction(w);
    if needGradient
        doseGradient = model.objectiveGradient(w);
    else
        doseGradient = [];
    end
else
    doseValue = model.optiProb.matRad_objectiveFunction(w, model.dij, model.cst);
    if needGradient
        doseGradient = model.optiProb.matRad_objectiveGradient(w, model.dij, model.cst);
    else
        doseGradient = [];
    end
end

constraintValue = 0;
constraintGradient = zeros(size(w));
maxViolation = 0;
constraintValues = zeros(0, 1);
if model.hasConstraints
    constraintValues = reshape(double(model.optiProb.matRad_constraintFunctions( ...
        w, model.dij, model.cst)), [], 1);
    lowerResidual = min(constraintValues - model.constraintLower, 0);
    upperResidual = max(constraintValues - model.constraintUpper, 0);
    residual = lowerResidual + upperResidual;
    finiteResidual = residual(isfinite(residual));
    if isempty(finiteResidual)
        finiteResidual = 0;
    end
    maxViolation = max(abs(finiteResidual));
    constraintValue = 0.5 * model.constraintPenalty * sum(finiteResidual .^ 2);
    if needGradient && model.constraintPenalty > 0
        J = model.optiProb.matRad_constraintJacobian(w, model.dij, model.cst);
        constraintGradient = model.constraintPenalty * (J' * residual);
    end
end

value = double(doseValue) + constraintValue;
if needGradient
    gradient = reshape(double(doseGradient), [], 1) + constraintGradient;
else
    gradient = [];
end
diagnostics = struct( ...
    'doseObjective', double(doseValue), ...
    'constraintMerit', constraintValue, ...
    'maxConstraintViolation', maxViolation, ...
    'constraintValues', constraintValues);
end
