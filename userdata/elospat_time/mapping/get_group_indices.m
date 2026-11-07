function ix = get_group_indices(groups, groupId)
% get_group_indices  Return validated global indices for one energy-layer group.

narginchk(2, 2);

if ~isstruct(groups) || isempty(groups) || ~isfield(groups, 'ix')
    error('get_group_indices:invalidInput', ...
        'groups must be a non-empty struct array with field ix.');
end

if ~isscalar(groupId) || ~isfinite(groupId) || groupId < 1 || floor(groupId) ~= groupId
    error('get_group_indices:invalidInput', ...
        'groupId must be a positive integer scalar.');
end

if groupId > numel(groups)
    error('get_group_indices:indexOutOfRange', ...
        'groupId %d exceeds number of groups (%d).', groupId, numel(groups));
end

ix = reshape(double(groups(groupId).ix), [], 1);

if isempty(ix)
    error('get_group_indices:emptyGroup', ...
        'groups(%d).ix is empty.', groupId);
end

if any(~isfinite(ix)) || any(ix < 1) || any(floor(ix) ~= ix)
    error('get_group_indices:invalidIndices', ...
        'groups(%d).ix must contain positive integer indices.', groupId);
end
end
