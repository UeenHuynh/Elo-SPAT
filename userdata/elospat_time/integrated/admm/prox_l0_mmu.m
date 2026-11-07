function [z, info] = prox_l0_mmu(a, lambda0, rho, minimumWeight)
%PROX_L0_MMU Exact separable proximal operator for L0 plus spot MMU.
%
% Solves independently for each i
%   min_z lambda0_i * 1[z ~= 0] + rho/2 * (z-a_i)^2
%   s.t. z in {0} union [minimumWeight_i, Inf).
%
% Equality of the off/on costs is resolved to zero. No layer-coverage
% repair is performed: Stage-1 layers are allowed groups, not mandatory
% active groups.

narginchk(4, 4);
a = reshape(double(a), [], 1);
n = numel(a);
lambda0 = expand_parameter(lambda0, n, 'lambda0');
minimumWeight = expand_parameter(minimumWeight, n, 'minimumWeight');
rho = double(rho);
if ~isscalar(rho) || ~isfinite(rho) || rho <= 0 || ...
        any(~isfinite(a)) || any(~isfinite(lambda0)) || any(lambda0 < 0) || ...
        any(~isfinite(minimumWeight)) || any(minimumWeight < 0)
    error('prox_l0_mmu:invalidInput', ...
        'Inputs must be finite; rho>0, lambda0>=0, minimumWeight>=0.');
end

onValue = max(a, minimumWeight);
offCost = 0.5 * rho .* (a .^ 2);
onCost = lambda0 + 0.5 * rho .* ((onValue - a) .^ 2);
active = onCost < offCost;
z = zeros(n, 1);
z(active) = onValue(active);

scale = max(1, max(abs([offCost; onCost])));
tieTolerance = 16 * eps(scale);
tie = abs(onCost - offCost) <= tieTolerance;
% Deterministic sparse tie rule, including numerically indistinguishable
% costs at the exact analytic threshold.
z(tie) = 0;
active(tie) = false;

info = struct('offCost', offCost, 'onCost', onCost, ...
    'onValue', onValue, 'active', active, 'tie', tie, ...
    'deltaOnMinusOff', onCost - offCost, ...
    'nActive', nnz(active), 'nExactZero', n - nnz(active), ...
    'nTies', nnz(tie), 'minimumPositive', NaN, ...
    'mmuFeasible', all(z == 0 | z >= minimumWeight));
if any(active)
    info.minimumPositive = min(z(active));
end
end

function x = expand_parameter(x, n, name)
x = reshape(double(x), [], 1);
if isscalar(x)
    x = repmat(x, n, 1);
end
if numel(x) ~= n
    error('prox_l0_mmu:dimensionMismatch', ...
        '%s has %d entries; expected 1 or %d.', name, numel(x), n);
end
end
