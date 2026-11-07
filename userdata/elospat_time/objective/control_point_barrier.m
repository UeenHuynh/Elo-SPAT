function [val, grad] = control_point_barrier(w, cpGroups, epsB)
% control_point_barrier  Compute the smooth control-point non-emptiness barrier.
%
% INPUTS
%   w        - spot-weight vector (N x 1 or 1 x N)
%   cpGroups - struct array with field .ix defining control-point membership
%   epsB     - smoothing constant (must be positive)
%
% OUTPUTS
%   val  - scalar barrier value
%   grad - gradient vector with the same size as w

narginchk(3, 3);

if ~isnumeric(w) || isempty(w)
    error('control_point_barrier:invalidInput', ...
        'w must be a non-empty numeric vector.');
end

if ~isvector(w)
    error('control_point_barrier:invalidInput', ...
        'w must be a vector.');
end

if ~isstruct(cpGroups) || isempty(cpGroups) || ~isfield(cpGroups, 'ix')
    error('control_point_barrier:invalidInput', ...
        'cpGroups must be a non-empty struct array with field ix.');
end

if ~isscalar(epsB) || ~isnumeric(epsB) || ~isfinite(epsB) || epsB <= 0
    error('control_point_barrier:invalidInput', ...
        'epsB must be a positive finite scalar.');
end

originalSize = size(w);
w = reshape(double(w), [], 1);
n = numel(w);

val = 0;
grad = zeros(n, 1);

for b = 1:numel(cpGroups)
    ix = reshape(double(cpGroups(b).ix), [], 1);

    if isempty(ix)
        error('control_point_barrier:emptyGroup', ...
            'cpGroups(%d).ix must not be empty.', b);
    end

    if any(~isfinite(ix)) || any(ix < 1) || any(ix > n) || any(floor(ix) ~= ix)
        error('control_point_barrier:indexOutOfRange', ...
            'cpGroups(%d).ix contains invalid indices for w of length %d.', b, n);
    end

    groupWeights = w(ix);
    smoothTerm = sum(groupWeights .^ 2) + epsB;

    val = val - log(sqrt(smoothTerm));
    grad(ix) = grad(ix) - groupWeights ./ smoothTerm;
end

grad = reshape(grad, originalSize);

end
