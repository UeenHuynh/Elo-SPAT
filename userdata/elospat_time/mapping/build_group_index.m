function [groups, cpGroups] = build_group_index(map)
% build_group_index  Aggregate mapped bixels into energy-layer and control-point groups.
%
% INPUT
%   map - N x 1 struct array from map_bixel_to_cp_energy_spot()
%
% OUTPUTS
%   groups   - G x 1 struct array with fields:
%              .ix
%              .nSpots
%              .energy
%              .cpIndex
%              .energyIndex
%   cpGroups - B x 1 struct array with fields:
%              .ix
%              .gantryAngle
%              .cpIndex

narginchk(1, 1);

validate_map_input(map);

cpIndex = reshape(double([map.cpIndex]), [], 1);
energyIndex = reshape(double([map.energyIndex]), [], 1);
energy = reshape(double([map.energy]), [], 1);
gantryAngle = reshape(double([map.gantryAngle]), [], 1);
globalIx = reshape(double([map.globalIx]), [], 1);

pairMatrix = [cpIndex, energyIndex];
uniquePairs = unique(pairMatrix, 'rows', 'stable');
nGroups = size(uniquePairs, 1);

groups = repmat(struct( ...
    'ix', [], ...
    'nSpots', 0, ...
    'energy', 0, ...
    'cpIndex', 0, ...
    'energyIndex', 0), nGroups, 1);

for g = 1:nGroups
    pairMask = pairMatrix(:, 1) == uniquePairs(g, 1) & ...
        pairMatrix(:, 2) == uniquePairs(g, 2);
    pairIx = globalIx(pairMask);
    pairEnergy = unique(energy(pairMask));

    if numel(pairEnergy) ~= 1
        error('build_group_index:inconsistentGroup', ...
            'Group (%d, %d) contains multiple energies.', ...
            uniquePairs(g, 1), uniquePairs(g, 2));
    end

    groups(g).ix = reshape(pairIx, [], 1);
    groups(g).nSpots = numel(pairIx);
    groups(g).energy = pairEnergy;
    groups(g).cpIndex = uniquePairs(g, 1);
    groups(g).energyIndex = uniquePairs(g, 2);
end

uniqueCp = unique(cpIndex, 'stable');
nCpGroups = numel(uniqueCp);

cpGroups = repmat(struct( ...
    'ix', [], ...
    'gantryAngle', 0, ...
    'cpIndex', 0), nCpGroups, 1);

for b = 1:nCpGroups
    cpMask = cpIndex == uniqueCp(b);
    cpIx = globalIx(cpMask);
    cpAngle = unique(gantryAngle(cpMask));

    if numel(cpAngle) ~= 1
        error('build_group_index:inconsistentControlPoint', ...
            'Control point %d contains multiple gantry angles.', uniqueCp(b));
    end

    cpGroups(b).ix = reshape(cpIx, [], 1);
    cpGroups(b).gantryAngle = cpAngle;
    cpGroups(b).cpIndex = uniqueCp(b);
end

end

function validate_map_input(map)
requiredFields = {'globalIx', 'cpIndex', 'energyIndex', 'energy', 'gantryAngle'};

if ~isstruct(map) || isempty(map)
    error('build_group_index:invalidInput', ...
        'map must be a non-empty struct array.');
end

for i = 1:numel(requiredFields)
    fieldName = requiredFields{i};
    if ~isfield(map, fieldName)
        error('build_group_index:missingField', ...
            'map.%s is required.', fieldName);
    end
end

globalIx = reshape(double([map.globalIx]), [], 1);
cpIndex = reshape(double([map.cpIndex]), [], 1);
energyIndex = reshape(double([map.energyIndex]), [], 1);

if any(~isfinite(globalIx)) || any(globalIx < 1) || any(floor(globalIx) ~= globalIx)
    error('build_group_index:invalidInput', ...
        'map.globalIx must contain positive integer indices.');
end

if numel(unique(globalIx)) ~= numel(globalIx)
    error('build_group_index:invalidInput', ...
        'map.globalIx must be unique.');
end

if any(~isfinite(cpIndex)) || any(cpIndex < 1) || any(floor(cpIndex) ~= cpIndex)
    error('build_group_index:invalidInput', ...
        'map.cpIndex must contain positive integer indices.');
end

if any(~isfinite(energyIndex)) || any(energyIndex < 1) || any(floor(energyIndex) ~= energyIndex)
    error('build_group_index:invalidInput', ...
        'map.energyIndex must contain positive integer indices.');
end
end
