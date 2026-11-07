function result = solve_admm_w_ipopt(doseModel, w0, center, rho, settings)
%SOLVE_ADMM_W_IPOPT Full IPOPT reference for the Variant-H w-subproblem.
%
% All benchmark strategies solve the same nonnegative objective
% F_matRad(w)+rho/2||w-center||^2. Clinical constraint semantics inside
% F_matRad (including any configured merit) are unchanged. Raw constraint
% violations are measured separately. Explicit IPOPT clinical constraints
% are deliberately not added here, because the inexact comparators do not
% solve constrained inner subproblems.

narginchk(4, 5);
if nargin < 5, settings = struct(); end
w0 = reshape(double(w0), [], 1);
center = reshape(double(center), [], 1);
referenceModel = struct();
referenceModel.nWeights = doseModel.nWeights;
referenceModel.objectiveFunction = @(w) objective_only( ...
    doseModel, w, center, rho);
referenceModel.objectiveGradient = @(w) gradient_only( ...
    doseModel, w, center, rho);
result = solve_weighted_matrad_ipopt( ...
    referenceModel, w0, zeros(size(w0)), settings);
result.strategy = 'full_ipopt_reference';
result.reason = result.statusMessage;
[~, ~, result.diagnostics] = evaluate_admm_w_subproblem( ...
    doseModel, result.w, center, rho, true);
end

function value = objective_only(model, w, center, rho)
value = evaluate_admm_w_subproblem(model, w, center, rho, false);
end

function gradient = gradient_only(model, w, center, rho)
[~, gradient] = evaluate_admm_w_subproblem(model, w, center, rho, true);
end
