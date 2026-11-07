function ix = get_control_point_indices(cpGroups, cpGroupId)
% get_control_point_indices  Return validated global indices for one control-point group.

narginchk(2, 2);

if ~isstruct(cpGroups) || isempty(cpGroups) || ~isfield(cpGroups, 'ix')
    error('get_control_point_indices:invalidInput', ...
        'cpGroups must be a non-empty struct array with field ix.');
end

if ~isscalar(cpGroupId) || ~isfinite(cpGroupId) || cpGroupId < 1 || floor(cpGroupId) ~= cpGroupId
    error('get_control_point_indices:invalidInput', ...
        'cpGroupId must be a positive integer scalar.');
end

if cpGroupId > numel(cpGroups)
    error('get_control_point_indices:indexOutOfRange', ...
        'cpGroupId %d exceeds number of control-point groups (%d).', cpGroupId, numel(cpGroups));
end

ix = reshape(double(cpGroups(cpGroupId).ix), [], 1);

if isempty(ix)
    error('get_control_point_indices:emptyGroup', ...
        'cpGroups(%d).ix is empty.', cpGroupId);
end

if any(~isfinite(ix)) || any(ix < 1) || any(floor(ix) ~= ix)
    error('get_control_point_indices:invalidIndices', ...
        'cpGroups(%d).ix must contain positive integer indices.', cpGroupId);
end
end
