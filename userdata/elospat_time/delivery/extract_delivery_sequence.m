function seq = extract_delivery_sequence(w, map, groups, cfg)
% extract_delivery_sequence  Build a deterministic delivered-layer sequence.
%
% INPUTS
%   w      - spot-weight vector
%   map    - mapped bixel struct array from map_bixel_to_cp_energy_spot
%   groups - energy-layer grouping struct array from build_group_index
%   cfg    - configuration struct with layerThreshold and spotThreshold
%
% OUTPUT
%   seq - struct array with one entry per active delivered layer:
%         .cpIndex
%         .gantryAngle
%         .energy
%         .energyIndex
%         .spotIndices
%         .spotWeights
%         .numActiveSpots
%         .layerWeightSum
%         .prevEnergy
%         .nextEnergy

narginchk(4, 4);

validate_inputs(w, map, groups, cfg);

w = reshape(double(w), [], 1);
layerThreshold = double(cfg.layerThreshold);
spotThreshold = double(cfg.spotThreshold);

cpIndexMap = reshape(double([map.cpIndex]), [], 1);
gantryAngleMap = reshape(double([map.gantryAngle]), [], 1);
nGroups = numel(groups);

layerMeta = repmat(struct( ...
    'cpIndex', 0, ...
    'gantryAngle', 0, ...
    'energy', 0, ...
    'energyIndex', 0, ...
    'spotIndices', [], ...
    'spotWeights', [], ...
    'numActiveSpots', 0, ...
    'layerWeightSum', 0, ...
    'prevEnergy', [], ...
    'nextEnergy', []), 0, 1);

for g = 1:nGroups
    ix = reshape(double(groups(g).ix), [], 1);
    groupWeights = w(ix);
    layerNorm = norm(groupWeights, 2);

    if layerNorm <= layerThreshold
        continue;
    end

    activeSpotMask = groupWeights > spotThreshold;
    activeSpotIndices = ix(activeSpotMask);
    activeSpotWeights = groupWeights(activeSpotMask);

    cpIx = double(groups(g).cpIndex);
    cpMask = cpIndexMap == cpIx;
    cpAngles = unique(gantryAngleMap(cpMask));

    if numel(cpAngles) ~= 1
        error('extract_delivery_sequence:inconsistentControlPoint', ...
            'Control point %d must map to exactly one gantry angle.', cpIx);
    end

    layerMeta(end + 1, 1) = struct( ... %#ok<AGROW>
        'cpIndex', cpIx, ...
        'gantryAngle', cpAngles, ...
        'energy', double(groups(g).energy), ...
        'energyIndex', double(groups(g).energyIndex), ...
        'spotIndices', reshape(activeSpotIndices, [], 1), ...
        'spotWeights', reshape(activeSpotWeights, [], 1), ...
        'numActiveSpots', numel(activeSpotIndices), ...
        'layerWeightSum', sum(activeSpotWeights), ...
        'prevEnergy', [], ...
        'nextEnergy', []);
end

if isempty(layerMeta)
    seq = repmat(struct( ...
        'cpIndex', 0, ...
        'gantryAngle', 0, ...
        'energy', 0, ...
        'energyIndex', 0, ...
        'spotIndices', [], ...
        'spotWeights', [], ...
        'numActiveSpots', 0, ...
        'layerWeightSum', 0, ...
        'prevEnergy', [], ...
        'nextEnergy', []), 0, 1);
    return;
end

orderMode = resolve_order_mode(cfg);
orderedCp = unique([[layerMeta.gantryAngle]', [layerMeta.cpIndex]'], 'rows', 'sorted');

seq = repmat(struct( ...
    'cpIndex', 0, ...
    'gantryAngle', 0, ...
    'energy', 0, ...
    'energyIndex', 0, ...
    'spotIndices', [], ...
    'spotWeights', [], ...
    'numActiveSpots', 0, ...
    'layerWeightSum', 0, ...
    'prevEnergy', [], ...
    'nextEnergy', []), 0, 1);

for c = 1:size(orderedCp, 1)
    cpAngle = orderedCp(c, 1);
    cpIx = orderedCp(c, 2);
    cpMask = [layerMeta.cpIndex]' == cpIx & [layerMeta.gantryAngle]' == cpAngle;
    cpLayers = layerMeta(cpMask);
    cpLayers = order_layers(cpLayers, orderMode, cfg);
    seq = [seq; cpLayers(:)]; %#ok<AGROW>
end

for i = 1:numel(seq)
    if i == 1
        seq(i).prevEnergy = [];
    else
        seq(i).prevEnergy = seq(i - 1).energy;
    end

    if i == numel(seq)
        seq(i).nextEnergy = [];
    else
        seq(i).nextEnergy = seq(i + 1).energy;
    end
end

end

function validate_inputs(w, map, groups, cfg)
if ~isnumeric(w) || isempty(w) || ~isvector(w)
    error('extract_delivery_sequence:invalidInput', ...
        'w must be a non-empty numeric vector.');
end

if ~isstruct(map) || isempty(map) || ...
        ~all(isfield(map, {'cpIndex', 'gantryAngle'}))
    error('extract_delivery_sequence:invalidInput', ...
        'map must be a non-empty struct array with cpIndex and gantryAngle.');
end

if ~isstruct(groups) || isempty(groups) || ...
        ~all(isfield(groups, {'ix', 'cpIndex', 'energy', 'energyIndex'}))
    error('extract_delivery_sequence:invalidInput', ...
        'groups must be a non-empty struct array with ix/cpIndex/energy/energyIndex.');
end

if ~isstruct(cfg) || ~isfield(cfg, 'layerThreshold') || ~isfield(cfg, 'spotThreshold')
    error('extract_delivery_sequence:invalidInput', ...
        'cfg.layerThreshold and cfg.spotThreshold are required.');
end

if ~isscalar(cfg.layerThreshold) || ~isnumeric(cfg.layerThreshold) || ...
        ~isfinite(cfg.layerThreshold) || cfg.layerThreshold < 0
    error('extract_delivery_sequence:invalidInput', ...
        'cfg.layerThreshold must be a non-negative finite scalar.');
end

if ~isscalar(cfg.spotThreshold) || ~isnumeric(cfg.spotThreshold) || ...
        ~isfinite(cfg.spotThreshold) || cfg.spotThreshold < 0
    error('extract_delivery_sequence:invalidInput', ...
        'cfg.spotThreshold must be a non-negative finite scalar.');
end
end

function orderMode = resolve_order_mode(cfg)
orderMode = 'descending';

if isfield(cfg, 'layerOrderMode') && ~isempty(cfg.layerOrderMode)
    orderMode = char(cfg.layerOrderMode);
elseif isfield(cfg, 'delivery') && isstruct(cfg.delivery) && ...
        isfield(cfg.delivery, 'layerOrderMode') && ~isempty(cfg.delivery.layerOrderMode)
    orderMode = char(cfg.delivery.layerOrderMode);
end
end

function orderedLayers = order_layers(cpLayers, orderMode, cfg)
energies = [cpLayers.energy]';
customOrder = [];

if strcmpi(orderMode, 'custom')
    if isfield(cfg, 'customEnergyOrder') && ~isempty(cfg.customEnergyOrder)
        customOrder = cfg.customEnergyOrder;
    elseif isfield(cfg, 'delivery') && isstruct(cfg.delivery) && ...
            isfield(cfg.delivery, 'customEnergyOrder') && ...
            ~isempty(cfg.delivery.customEnergyOrder)
        customOrder = cfg.delivery.customEnergyOrder;
    else
        error('extract_delivery_sequence:missingCustomOrder', ...
            'A custom energy order must be supplied for order mode ''custom''.');
    end
end

orderedEnergies = order_layers_within_control_point(energies, orderMode, customOrder);
order = zeros(numel(energies), 1);
usedMask = false(numel(energies), 1);

for i = 1:numel(orderedEnergies)
    tolerance = max(1e-6, abs(orderedEnergies(i)) * 1e-8);
    match = find(~usedMask & abs(energies - orderedEnergies(i)) <= tolerance, 1, 'first');
    if isempty(match)
        error('extract_delivery_sequence:energyReorderingFailed', ...
            'Could not match reordered energy %.8g back to the CP layer list.', ...
            orderedEnergies(i));
    end
    order(i) = match;
    usedMask(match) = true;
end

orderedLayers = cpLayers(order);
end
