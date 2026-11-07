function validate_mapping(map, dij, groups, cpGroups)
% validate_mapping  Validate consistency between map, dij, groups, and cpGroups.
%
% INPUTS
%   map      - struct array from map_bixel_to_cp_energy_spot()
%   dij      - matRad dose influence matrix
%   groups   - struct array from build_group_index()
%   cpGroups - struct array from build_group_index()

narginchk(4, 4);

validate_map_against_dij(map, dij);
validate_energy_groups(groups, dij.totalNumOfBixels);
validate_control_point_groups(cpGroups, dij.totalNumOfBixels);

mapGlobalIx = reshape(double([map.globalIx]), [], 1);
groupCoverage = collect_group_coverage(groups, 'groups');
cpCoverage = collect_group_coverage(cpGroups, 'cpGroups');

if ~isequal(groupCoverage, mapGlobalIx)
    missingIx = setdiff(mapGlobalIx, groupCoverage);
    extraIx = setdiff(groupCoverage, mapGlobalIx);
    duplicateIx = find_duplicate_indices(groupCoverage);
    error('validate_mapping:groupCoverageMismatch', ...
        ['Energy-layer groups do not cover map.globalIx exactly once. ', ...
         'Missing: %s | Extra: %s | Duplicates: %s'], ...
        format_index_list(missingIx), format_index_list(extraIx), format_index_list(duplicateIx));
end

if ~isequal(cpCoverage, mapGlobalIx)
    missingIx = setdiff(mapGlobalIx, cpCoverage);
    extraIx = setdiff(cpCoverage, mapGlobalIx);
    duplicateIx = find_duplicate_indices(cpCoverage);
    error('validate_mapping:controlPointCoverageMismatch', ...
        ['Control-point groups do not cover map.globalIx exactly once. ', ...
         'Missing: %s | Extra: %s | Duplicates: %s'], ...
        format_index_list(missingIx), format_index_list(extraIx), format_index_list(duplicateIx));
end

validate_group_metadata_consistency(map, groups, cpGroups);

end

function validate_map_against_dij(map, dij)
if ~isstruct(map) || isempty(map)
    error('validate_mapping:invalidInput', ...
        'map must be a non-empty struct array.');
end

if ~isstruct(dij) || ~isfield(dij, 'totalNumOfBixels')
    error('validate_mapping:invalidInput', ...
        'dij.totalNumOfBixels is required.');
end

requiredMapFields = {'globalIx', 'cpIndex', 'energyIndex', 'energy', 'gantryAngle'};
for i = 1:numel(requiredMapFields)
    if ~isfield(map, requiredMapFields{i})
        error('validate_mapping:missingField', ...
            'map.%s is required.', requiredMapFields{i});
    end
end

nMap = numel(map);
nDij = double(dij.totalNumOfBixels);
if nMap ~= nDij
    error('validate_mapping:lengthMismatch', ...
        'numel(map)=%d does not match dij.totalNumOfBixels=%d.', nMap, nDij);
end

globalIx = reshape(double([map.globalIx]), [], 1);
expectedIx = reshape((1:nDij)', [], 1);

if ~isequal(globalIx, expectedIx)
    missingIx = setdiff(expectedIx, globalIx);
    extraIx = setdiff(globalIx, expectedIx);
    duplicateIx = find_duplicate_indices(globalIx);
    error('validate_mapping:globalIxMismatch', ...
        ['map.globalIx must equal 1:N consecutively. ', ...
         'Missing: %s | Extra: %s | Duplicates: %s'], ...
        format_index_list(missingIx), format_index_list(extraIx), format_index_list(duplicateIx));
end
end

function validate_energy_groups(groups, nBixels)
if ~isstruct(groups) || isempty(groups)
    error('validate_mapping:invalidInput', ...
        'groups must be a non-empty struct array.');
end

requiredFields = {'ix', 'nSpots', 'energy', 'cpIndex', 'energyIndex'};
for i = 1:numel(requiredFields)
    if ~isfield(groups, requiredFields{i})
        error('validate_mapping:missingField', ...
            'groups.%s is required.', requiredFields{i});
    end
end

for i = 1:numel(groups)
    ix = reshape(double(groups(i).ix), [], 1);
    if isempty(ix)
        error('validate_mapping:emptyGroup', ...
            'groups(%d).ix must not be empty.', i);
    end
    if any(ix < 1) || any(ix > nBixels) || any(floor(ix) ~= ix)
        error('validate_mapping:indexOutOfRange', ...
            'groups(%d).ix contains indices outside 1:%d.', i, nBixels);
    end
    if groups(i).nSpots ~= numel(ix)
        error('validate_mapping:nSpotsMismatch', ...
            'groups(%d).nSpots=%d does not match numel(ix)=%d.', ...
            i, groups(i).nSpots, numel(ix));
    end
end
end

function validate_control_point_groups(cpGroups, nBixels)
if ~isstruct(cpGroups) || isempty(cpGroups)
    error('validate_mapping:invalidInput', ...
        'cpGroups must be a non-empty struct array.');
end

requiredFields = {'ix', 'gantryAngle', 'cpIndex'};
for i = 1:numel(requiredFields)
    if ~isfield(cpGroups, requiredFields{i})
        error('validate_mapping:missingField', ...
            'cpGroups.%s is required.', requiredFields{i});
    end
end

for i = 1:numel(cpGroups)
    ix = reshape(double(cpGroups(i).ix), [], 1);
    if isempty(ix)
        error('validate_mapping:emptyControlPointGroup', ...
            'cpGroups(%d).ix must not be empty.', i);
    end
    if any(ix < 1) || any(ix > nBixels) || any(floor(ix) ~= ix)
        error('validate_mapping:indexOutOfRange', ...
            'cpGroups(%d).ix contains indices outside 1:%d.', i, nBixels);
    end
end
end

function validate_group_metadata_consistency(map, groups, cpGroups)
mapCp = reshape(double([map.cpIndex]), [], 1);
mapEnergyIndex = reshape(double([map.energyIndex]), [], 1);
mapEnergy = reshape(double([map.energy]), [], 1);
mapAngle = reshape(double([map.gantryAngle]), [], 1);

for i = 1:numel(groups)
    ix = get_group_indices(groups, i);

    if any(mapCp(ix) ~= groups(i).cpIndex)
        error('validate_mapping:groupCpMismatch', ...
            'groups(%d).cpIndex does not match mapped cpIndex values.', i);
    end
    if any(mapEnergyIndex(ix) ~= groups(i).energyIndex)
        error('validate_mapping:groupEnergyIndexMismatch', ...
            'groups(%d).energyIndex does not match mapped energyIndex values.', i);
    end
    if any(abs(mapEnergy(ix) - groups(i).energy) > 1e-8)
        error('validate_mapping:groupEnergyMismatch', ...
            'groups(%d).energy does not match mapped energy values.', i);
    end
end

for i = 1:numel(cpGroups)
    ix = get_control_point_indices(cpGroups, i);

    if any(mapCp(ix) ~= cpGroups(i).cpIndex)
        error('validate_mapping:controlPointCpMismatch', ...
            'cpGroups(%d).cpIndex does not match mapped cpIndex values.', i);
    end
    if any(abs(mapAngle(ix) - cpGroups(i).gantryAngle) > 1e-8)
        error('validate_mapping:controlPointAngleMismatch', ...
            'cpGroups(%d).gantryAngle does not match mapped gantry angles.', i);
    end
end
end

function coverage = collect_group_coverage(groupStruct, label)
coverage = [];

for i = 1:numel(groupStruct)
    ix = reshape(double(groupStruct(i).ix), [], 1);
    coverage = [coverage; ix]; %#ok<AGROW>
end

coverage = sort(coverage);

if isempty(coverage)
    error('validate_mapping:emptyCoverage', ...
        '%s coverage is empty.', label);
end
end

function duplicateIx = find_duplicate_indices(ix)
if isempty(ix)
    duplicateIx = [];
    return;
end

uniqueIx = unique(ix);
counts = histc(ix, uniqueIx); %#ok<HISTC>
duplicateIx = uniqueIx(counts > 1);
end

function text = format_index_list(ix)
if isempty(ix)
    text = '[]';
    return;
end

ix = reshape(double(ix), 1, []);
maxItems = min(numel(ix), 10);
text = mat2str(ix(1:maxItems));
if numel(ix) > maxItems
    text = sprintf('%s ... (%d total)', text, numel(ix));
end
end
