function spotMask = build_spot_mask(w, threshold)
% build_spot_mask  Identify active spots using a scalar weight threshold.
%
% INPUTS
%   w         - spot-weight vector
%   threshold - non-negative scalar activation threshold
%
% OUTPUT
%   spotMask  - logical vector (N x 1), true = keep spot

narginchk(2, 2);

if ~isnumeric(w) || isempty(w) || ~isvector(w) || any(~isfinite(w(:)))
    error('build_spot_mask:invalidInput', ...
        'w must be a non-empty finite numeric vector.');
end

if ~isscalar(threshold) || ~isnumeric(threshold) || ~isfinite(threshold) || threshold < 0
    error('build_spot_mask:invalidInput', ...
        'threshold must be a non-negative finite scalar.');
end

spotMask = reshape(double(w), [], 1) > threshold;

end
