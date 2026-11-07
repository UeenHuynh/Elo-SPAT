function [val, grad] = spot_sparsity_penalty(w, epsS)
% spot_sparsity_penalty  Compute the smooth L0-like spot sparsity penalty.
%
% INPUTS
%   w    - spot-weight vector (N x 1 or 1 x N), assumed non-negative
%   epsS - smoothing constant (must be positive)
%
% OUTPUTS
%   val  - scalar penalty value
%   grad - gradient vector with the same size as w

narginchk(2, 2);

if ~isnumeric(w) || isempty(w)
    error('spot_sparsity_penalty:invalidInput', ...
        'w must be a non-empty numeric vector.');
end

if ~isvector(w)
    error('spot_sparsity_penalty:invalidInput', ...
        'w must be a vector.');
end

if any(~isfinite(w(:)))
    error('spot_sparsity_penalty:invalidInput', ...
        'w must contain only finite values.');
end

if any(w(:) < 0)
    error('spot_sparsity_penalty:invalidInput', ...
        'w must be non-negative for the smooth spot sparsity penalty.');
end

if ~isscalar(epsS) || ~isnumeric(epsS) || ~isfinite(epsS) || epsS <= 0
    error('spot_sparsity_penalty:invalidInput', ...
        'epsS must be a positive finite scalar.');
end

originalSize = size(w);
w = double(w);

val = sum(w(:) ./ (w(:) + epsS));
grad = epsS ./ (w + epsS) .^ 2;
grad = reshape(grad, originalSize);

end
