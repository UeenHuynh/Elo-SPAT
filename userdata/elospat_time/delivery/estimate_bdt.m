function metrics = estimate_bdt(seq, machine)
% estimate_bdt  Estimate machine-specific beam delivery time components.
%
% INPUTS
%   seq     - delivery sequence struct array from extract_delivery_sequence
%   machine - machine timing model struct from default_machine_model
%
% OUTPUT
%   metrics - struct with fields:
%             .T_energy, .T_spot, .T_irradiation, .T_gantry, .T_buffer
%             .T_total, .nActiveLayers, .nActiveSpots, .nSwitchUp
%             .totalSwitchMagnitude

narginchk(2, 2);

validate_seq(seq);
validate_machine(machine);

if isempty(seq)
    metrics = struct( ...
        'T_energy', 0, ...
        'T_spot', 0, ...
        'T_irradiation', 0, ...
        'T_gantry', 0, ...
        'T_buffer', 0, ...
        'T_total', 0, ...
        'nActiveLayers', 0, ...
        'nActiveSpots', 0, ...
        'nSwitchUp', 0, ...
        'totalSwitchMagnitude', 0);
    return;
end

nLayers = numel(seq);
cpIndex = [seq.cpIndex]';
gantryAngle = [seq.gantryAngle]';
energy = [seq.energy]';
numActiveSpots = [seq.numActiveSpots]';
layerWeightSum = [seq.layerWeightSum]';

T_energy = 0;
T_gantry = 0;
nSwitchUp = 0;
totalSwitchMagnitude = 0;

for i = 2:nLayers
    deltaEnergy = energy(i) - energy(i - 1);
    totalSwitchMagnitude = totalSwitchMagnitude + abs(deltaEnergy);

    if deltaEnergy > 0
        T_energy = T_energy + machine.energy.upBase + ...
            machine.energy.perMeVUp * deltaEnergy;
        nSwitchUp = nSwitchUp + 1;
    elseif deltaEnergy < 0
        T_energy = T_energy + machine.energy.downBase + ...
            machine.energy.perMeVDown * abs(deltaEnergy);
    end

    if cpIndex(i) ~= cpIndex(i - 1)
        theta = abs(gantryAngle(i) - gantryAngle(i - 1));
        T_gantry = T_gantry + estimate_gantry_time(theta, machine);
    end
end

T_spot = sum(max(numActiveSpots - 1, 0)) * machine.spot.switchTime;
T_irradiation = sum(layerWeightSum) / machine.irradiation.muRate;
nControlPoints = numel(unique(cpIndex, 'stable'));
T_buffer = nLayers * machine.buffer.perLayer + ...
    nControlPoints * machine.buffer.perControlPoint;

T_total = T_energy + T_spot + T_irradiation + T_gantry + T_buffer;

metrics = struct( ...
    'T_energy', T_energy, ...
    'T_spot', T_spot, ...
    'T_irradiation', T_irradiation, ...
    'T_gantry', T_gantry, ...
    'T_buffer', T_buffer, ...
    'T_total', T_total, ...
    'nActiveLayers', nLayers, ...
    'nActiveSpots', sum(numActiveSpots), ...
    'nSwitchUp', nSwitchUp, ...
    'totalSwitchMagnitude', totalSwitchMagnitude);

componentSum = metrics.T_energy + metrics.T_spot + metrics.T_irradiation + ...
    metrics.T_gantry + metrics.T_buffer;
if abs(metrics.T_total - componentSum) > 1e-10
    error('estimate_bdt:internalConsistency', ...
        'T_total must equal the sum of all BDT components.');
end

end

function validate_seq(seq)
requiredFields = {'cpIndex', 'gantryAngle', 'energy', 'numActiveSpots', 'layerWeightSum'};

if ~isstruct(seq)
    error('estimate_bdt:invalidInput', ...
        'seq must be a struct array.');
end

if isempty(seq)
    return;
end

for i = 1:numel(requiredFields)
    if ~isfield(seq, requiredFields{i})
        error('estimate_bdt:missingField', ...
            'seq.%s is required.', requiredFields{i});
    end
end

numericFields = {'cpIndex', 'gantryAngle', 'energy', 'numActiveSpots', 'layerWeightSum'};
for i = 1:numel(seq)
    for j = 1:numel(numericFields)
        value = seq(i).(numericFields{j});
        if isempty(value) || ~isscalar(value) || ~isnumeric(value) || ~isfinite(value)
            error('estimate_bdt:invalidInput', ...
                'seq(%d).%s must be a finite numeric scalar.', i, numericFields{j});
        end
    end
end
end

function validate_machine(machine)
requiredPaths = { ...
    {'energy', 'upBase'}, {'energy', 'downBase'}, ...
    {'energy', 'perMeVUp'}, {'energy', 'perMeVDown'}, ...
    {'spot', 'switchTime'}, ...
    {'irradiation', 'muRate'}, ...
    {'gantry', 'maxSpeed'}, {'gantry', 'maxAccel'}, ...
    {'buffer', 'perLayer'}, {'buffer', 'perControlPoint'}};

if ~isstruct(machine)
    error('estimate_bdt:invalidInput', ...
        'machine must be a struct.');
end

for i = 1:numel(requiredPaths)
    path = requiredPaths{i};
    if ~isfield(machine, path{1}) || ~isstruct(machine.(path{1})) || ...
            ~isfield(machine.(path{1}), path{2})
        error('estimate_bdt:missingField', ...
            'machine.%s.%s is required.', path{1}, path{2});
    end

    value = machine.(path{1}).(path{2});
    if ~isscalar(value) || ~isnumeric(value) || ~isfinite(value) || value < 0
        error('estimate_bdt:invalidInput', ...
            'machine.%s.%s must be a non-negative finite scalar.', path{1}, path{2});
    end
end

if machine.irradiation.muRate <= 0
    error('estimate_bdt:invalidInput', ...
        'machine.irradiation.muRate must be strictly positive.');
end

if machine.gantry.maxSpeed <= 0 || machine.gantry.maxAccel <= 0
    error('estimate_bdt:invalidInput', ...
        'machine.gantry.maxSpeed and maxAccel must be strictly positive.');
end
end
