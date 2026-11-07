function [dijReduced, reducedToFullIdx, fullToReducedIdx, mapReduced] = reduce_dij_to_support(dij, supportMask, map)
% reduce_dij_to_support  Reduce all bixel-indexed dij fields consistently.
%
% The returned index maps make support-preserving re-optimization explicit:
%   wFull(supportMask) = wReduced;
%   wFull(~supportMask) = 0.

narginchk(2, 3);
if nargin < 3
    map = [];
end

nBixels = double(dij.totalNumOfBixels);
supportMask = reshape(logical(supportMask), [], 1);
if numel(supportMask) ~= nBixels
    error('reduce_dij_to_support:dimensionMismatch', ...
        'supportMask length %d does not match dij.totalNumOfBixels=%d.', ...
        numel(supportMask), nBixels);
end

reducedToFullIdx = find(supportMask);
fullToReducedIdx = zeros(nBixels, 1);
fullToReducedIdx(reducedToFullIdx) = (1:numel(reducedToFullIdx))';

[dijReduced, mapReduced] = cut_dij_by_mask(dij, supportMask, map);
dijReduced.reducedToFullIdx = reducedToFullIdx;
dijReduced.fullToReducedIdx = fullToReducedIdx;
end
