function [value, gradient, diagnostics] = evaluate_admm_w_subproblem( ...
    doseModel, w, center, rho, needGradient)
%EVALUATE_ADMM_W_SUBPROBLEM Exact matRad term plus ADMM quadratic.

narginchk(4, 5);
if nargin < 5
    needGradient = nargout >= 2;
end
w = reshape(double(w), [], 1);
center = reshape(double(center), [], 1);
if numel(w) ~= numel(center) || numel(w) ~= doseModel.nWeights
    error('evaluate_admm_w_subproblem:dimensionMismatch', ...
        'w and center must match doseModel.nWeights.');
end
[doseValue, doseGradient, doseDiagnostics] = ...
    evaluate_matrad_dose_model(doseModel, w, needGradient);
d = w - center;
quadratic = 0.5 * rho * (d' * d);
value = doseValue + quadratic;
if needGradient
    gradient = doseGradient + rho * d;
else
    gradient = [];
end
diagnostics = struct('doseObjective', doseValue, ...
    'quadratic', quadratic, ...
    'maxConstraintViolation', doseDiagnostics.maxConstraintViolation);
end
