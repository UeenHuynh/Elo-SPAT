function bixelMask = expand_layer_mask_to_bixel_mask(layerMask, groups, totalNumOfBixels)
% expand_layer_mask_to_bixel_mask  Convert group-level keep mask to bixel-level mask.
%
% INPUTS
%   layerMask        - logical or 0/1 vector of size nGroups x 1
%   groups           - struct array with field .ix
%   totalNumOfBixels - positive integer scalar
%
% OUTPUT
%   bixelMask - logical column vector of size totalNumOfBixels x 1

narginchk(3, 3);

if ~isstruct(groups) || isempty(groups) || ~isfield(groups, 'ix')
    error('expand_layer_mask_to_bixel_mask:invalidInput', ...
        'groups must be a non-empty struct array with field ix.');
end

nGroups = numel(groups);

if ~(islogical(layerMask) || isnumeric(layerMask)) || isempty(layerMask) || ~isvector(layerMask)
    error('expand_layer_mask_to_bixel_mask:invalidInput', ...
        'layerMask must be a non-empty logical or numeric vector.');
end

layerMask = reshape(layerMask, [], 1);
if numel(layerMask) ~= nGroups
    error('expand_layer_mask_to_bixel_mask:dimensionMismatch', ...
        'layerMask must have length %d (one entry per group), got %d.', ...
        nGroups, numel(layerMask));
end

if isnumeric(layerMask)
    if any(~isfinite(layerMask)) || any(layerMask ~= 0 & layerMask ~= 1)
        error('expand_layer_mask_to_bixel_mask:invalidInput', ...
            'Numeric layerMask must contain only 0/1 values.');
    end
    layerMask = logical(layerMask);
end

if ~isscalar(totalNumOfBixels) || ~isnumeric(totalNumOfBixels) || ...
        ~isfinite(totalNumOfBixels) || totalNumOfBixels < 1 || ...
        floor(totalNumOfBixels) ~= totalNumOfBixels
    error('expand_layer_mask_to_bixel_mask:invalidInput', ...
        'totalNumOfBixels must be a positive integer scalar.');
end

totalNumOfBixels = double(totalNumOfBixels);
bixelMask = false(totalNumOfBixels, 1);

for g = 1:nGroups
    ix = reshape(double(groups(g).ix), [], 1);

    if isempty(ix)
        error('expand_layer_mask_to_bixel_mask:emptyGroup', ...
            'groups(%d).ix must not be empty.', g);
    end

    if any(~isfinite(ix)) || any(ix < 1) || any(ix > totalNumOfBixels) || any(floor(ix) ~= ix)
        error('expand_layer_mask_to_bixel_mask:indexOutOfRange', ...
            'groups(%d).ix contains invalid indices for totalNumOfBixels=%d.', ...
            g, totalNumOfBixels);
    end

    if layerMask(g)
        bixelMask(ix) = true;
    end
end

end
