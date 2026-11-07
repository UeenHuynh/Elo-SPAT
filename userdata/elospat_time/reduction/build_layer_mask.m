function layerMask = build_layer_mask(w, groups, threshold)
% build_layer_mask  Identify active energy layers using group L2 norm.
%
% INPUTS
%   w         - spot-weight vector
%   groups    - struct array with field .ix
%   threshold - non-negative scalar activation threshold
%
% OUTPUT
%   layerMask - logical vector (nGroups x 1), true = keep group

narginchk(3, 3);

if ~isnumeric(w) || isempty(w) || ~isvector(w) || any(~isfinite(w(:)))
    error('build_layer_mask:invalidInput', ...
        'w must be a non-empty finite numeric vector.');
end

if ~isstruct(groups) || isempty(groups) || ~isfield(groups, 'ix')
    error('build_layer_mask:invalidInput', ...
        'groups must be a non-empty struct array with field ix.');
end

if ~isscalar(threshold) || ~isnumeric(threshold) || ~isfinite(threshold) || threshold < 0
    error('build_layer_mask:invalidInput', ...
        'threshold must be a non-negative finite scalar.');
end

w = reshape(double(w), [], 1);
n = numel(w);
nGroups = numel(groups);
layerMask = false(nGroups, 1);

for g = 1:nGroups
    ix = reshape(double(groups(g).ix), [], 1);

    if isempty(ix)
        error('build_layer_mask:emptyGroup', ...
            'groups(%d).ix must not be empty.', g);
    end

    if any(~isfinite(ix)) || any(ix < 1) || any(ix > n) || any(floor(ix) ~= ix)
        error('build_layer_mask:indexOutOfRange', ...
            'groups(%d).ix contains invalid indices for w of length %d.', g, n);
    end

    layerMask(g) = norm(w(ix), 2) > threshold;
end

end
