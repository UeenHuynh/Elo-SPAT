function [val, grad] = l21half_group_penalty(w, groups, epsEL)
% l21half_group_penalty  Compute the ELO-SPAT L2,1/2 group sparsity penalty.
%
% INPUTS
%   w      - spot-weight vector (N x 1 or 1 x N)
%   groups - struct array with field .ix defining group membership
%   epsEL  - smoothing constant (must be positive)
%
% OUTPUTS
%   val  - scalar penalty value
%   grad - gradient vector with the same size as w

narginchk(3, 3);

if ~isnumeric(w) || isempty(w)
    error('l21half_group_penalty:invalidInput', ...
        'w must be a non-empty numeric vector.');
end

if ~isvector(w)
    error('l21half_group_penalty:invalidInput', ...
        'w must be a vector.');
end

if ~isstruct(groups) || isempty(groups) || ~isfield(groups, 'ix')
    error('l21half_group_penalty:invalidInput', ...
        'groups must be a non-empty struct array with field ix.');
end

if ~isscalar(epsEL) || ~isnumeric(epsEL) || ~isfinite(epsEL) || epsEL <= 0
    error('l21half_group_penalty:invalidInput', ...
        'epsEL must be a positive finite scalar.');
end

originalSize = size(w);
w = reshape(double(w), [], 1);
n = numel(w);

val = 0;
grad = zeros(n, 1);

for g = 1:numel(groups)
    ix = reshape(double(groups(g).ix), [], 1);

    if isempty(ix)
        error('l21half_group_penalty:emptyGroup', ...
            'groups(%d).ix must not be empty.', g);
    end

    if any(~isfinite(ix)) || any(ix < 1) || any(ix > n) || any(floor(ix) ~= ix)
        error('l21half_group_penalty:indexOutOfRange', ...
            'groups(%d).ix contains invalid indices for w of length %d.', g, n);
    end

    groupWeights = w(ix);
    groupNormSq = sum(groupWeights .^ 2);
    smoothTerm = groupNormSq + epsEL;

    val = val + smoothTerm .^ 0.25;
    grad(ix) = grad(ix) + 0.5 .* groupWeights .* smoothTerm .^ (-0.75);
end

grad = reshape(grad, originalSize);

end
