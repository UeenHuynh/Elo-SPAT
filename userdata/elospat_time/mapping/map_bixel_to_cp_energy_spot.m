function map = map_bixel_to_cp_energy_spot(dij, stf, pln)
% map_bixel_to_cp_energy_spot  Map each dij column to control-point and energy-layer identity.
%
% INPUTS
%   dij - matRad dose influence matrix with beamNum, rayNum, bixelNum
%   stf - matRad steering information struct array
%   pln - matRad plan structure with pln.propStf.gantryAngles
%
% OUTPUT
%   map - N x 1 struct array with fields:
%         .globalIx
%         .beamNum
%         .rayNum
%         .bixelNum
%         .gantryAngle
%         .energy
%         .cpIndex
%         .energyIndex

narginchk(3, 3);

validate_dij_inputs(dij);
validate_stf_inputs(stf);
validate_pln_inputs(pln, numel(stf));

nBixels = double(dij.totalNumOfBixels);
beamNum = reshape(double(dij.beamNum), [], 1);
rayNum = reshape(double(dij.rayNum), [], 1);
bixelNum = reshape(double(dij.bixelNum), [], 1);
globalIx = reshape((1:nBixels)', [], 1);

gantryAngles = reshape(double(pln.propStf.gantryAngles), [], 1);
beamEnergyTables = build_beam_energy_tables(stf);

gantryAngleOut = zeros(nBixels, 1);
energyOut = zeros(nBixels, 1);
cpIndexOut = beamNum;
energyIndexOut = zeros(nBixels, 1);

for i = 1:nBixels
    b = beamNum(i);
    r = rayNum(i);
    s = bixelNum(i);

    validate_mapping_indices(stf, b, r, s, i);

    gantryAngleOut(i) = resolve_gantry_angle(stf, gantryAngles, b);
    energyOut(i) = double(stf(b).ray(r).energy(s));
    energyIndexOut(i) = resolve_energy_index(beamEnergyTables{b}, energyOut(i), b, i);
end

map = struct( ...
    'globalIx', num2cell(globalIx), ...
    'beamNum', num2cell(beamNum), ...
    'rayNum', num2cell(rayNum), ...
    'bixelNum', num2cell(bixelNum), ...
    'gantryAngle', num2cell(gantryAngleOut), ...
    'energy', num2cell(energyOut), ...
    'cpIndex', num2cell(cpIndexOut), ...
    'energyIndex', num2cell(energyIndexOut));

end

function validate_dij_inputs(dij)
requiredFields = {'totalNumOfBixels', 'beamNum', 'rayNum', 'bixelNum'};

if ~isstruct(dij)
    error('map_bixel_to_cp_energy_spot:invalidInput', ...
        'dij must be a struct.');
end

for i = 1:numel(requiredFields)
    fieldName = requiredFields{i};
    if ~isfield(dij, fieldName)
        error('map_bixel_to_cp_energy_spot:missingField', ...
            'dij.%s is required.', fieldName);
    end
end

nBixels = double(dij.totalNumOfBixels);
if ~isscalar(nBixels) || ~isfinite(nBixels) || nBixels < 0 || floor(nBixels) ~= nBixels
    error('map_bixel_to_cp_energy_spot:invalidInput', ...
        'dij.totalNumOfBixels must be a non-negative integer scalar.');
end

vectorFields = {'beamNum', 'rayNum', 'bixelNum'};
for i = 1:numel(vectorFields)
    fieldName = vectorFields{i};
    values = reshape(double(dij.(fieldName)), [], 1);

    if numel(values) ~= nBixels
        error('map_bixel_to_cp_energy_spot:invalidInput', ...
            'dij.%s must have %d entries (got %d).', fieldName, nBixels, numel(values));
    end

    if any(~isfinite(values)) || any(values < 1) || any(floor(values) ~= values)
        error('map_bixel_to_cp_energy_spot:invalidInput', ...
            'dij.%s must contain positive integer indices.', fieldName);
    end
end
end

function validate_stf_inputs(stf)
if ~isstruct(stf) || isempty(stf)
    error('map_bixel_to_cp_energy_spot:invalidInput', ...
        'stf must be a non-empty struct array.');
end

for b = 1:numel(stf)
    if ~isfield(stf(b), 'ray') || ~isstruct(stf(b).ray) || isempty(stf(b).ray)
        error('map_bixel_to_cp_energy_spot:invalidInput', ...
            'stf(%d).ray must be a non-empty struct array.', b);
    end

    for r = 1:numel(stf(b).ray)
        if ~isfield(stf(b).ray(r), 'energy') || isempty(stf(b).ray(r).energy)
            error('map_bixel_to_cp_energy_spot:invalidInput', ...
                'stf(%d).ray(%d).energy must be non-empty.', b, r);
        end
    end
end
end

function validate_pln_inputs(pln, nBeams)
if ~isstruct(pln) || ~isfield(pln, 'propStf') || ~isstruct(pln.propStf)
    error('map_bixel_to_cp_energy_spot:invalidInput', ...
        'pln.propStf is required.');
end

if ~isfield(pln.propStf, 'gantryAngles') || isempty(pln.propStf.gantryAngles)
    error('map_bixel_to_cp_energy_spot:missingField', ...
        'pln.propStf.gantryAngles is required.');
end

gantryAngles = reshape(double(pln.propStf.gantryAngles), [], 1);
if numel(gantryAngles) < nBeams
    error('map_bixel_to_cp_energy_spot:invalidInput', ...
        'pln.propStf.gantryAngles must provide at least %d beam angles.', nBeams);
end
end

function beamEnergyTables = build_beam_energy_tables(stf)
nBeams = numel(stf);
beamEnergyTables = cell(nBeams, 1);

for b = 1:nBeams
    allEnergies = [];
    for r = 1:numel(stf(b).ray)
        rayEnergies = reshape(double(stf(b).ray(r).energy), [], 1);
        allEnergies = [allEnergies; rayEnergies]; %#ok<AGROW>
    end

    if isempty(allEnergies)
        error('map_bixel_to_cp_energy_spot:invalidInput', ...
            'Beam %d has no energies in stf.', b);
    end

    beamEnergyTables{b} = sort(unique(allEnergies), 'descend');
end
end

function validate_mapping_indices(stf, beamIx, rayIx, bixelIx, globalIx)
if beamIx > numel(stf)
    error('map_bixel_to_cp_energy_spot:indexOutOfRange', ...
        'dij.beamNum(%d)=%d exceeds number of beams (%d).', ...
        globalIx, beamIx, numel(stf));
end

if rayIx > numel(stf(beamIx).ray)
    error('map_bixel_to_cp_energy_spot:indexOutOfRange', ...
        'dij.rayNum(%d)=%d exceeds number of rays in beam %d (%d).', ...
        globalIx, rayIx, beamIx, numel(stf(beamIx).ray));
end

if bixelIx > numel(stf(beamIx).ray(rayIx).energy)
    error('map_bixel_to_cp_energy_spot:indexOutOfRange', ...
        'dij.bixelNum(%d)=%d exceeds number of bixels in beam %d ray %d (%d).', ...
        globalIx, bixelIx, beamIx, rayIx, numel(stf(beamIx).ray(rayIx).energy));
end
end

function gantryAngle = resolve_gantry_angle(stf, gantryAngles, beamIx)
if isfield(stf(beamIx), 'gantryAngle') && ~isempty(stf(beamIx).gantryAngle)
    gantryAngle = double(stf(beamIx).gantryAngle);
else
    gantryAngle = gantryAngles(beamIx);
end
end

function energyIndex = resolve_energy_index(energyTable, energyValue, beamIx, globalIx)
tolerance = max(1e-6, abs(energyValue) * 1e-8);
energyIx = find(abs(energyTable - energyValue) <= tolerance, 1, 'first');

if isempty(energyIx)
    error('map_bixel_to_cp_energy_spot:energyLookupFailed', ...
        'Could not assign energyIndex for globalIx=%d in beam %d at energy %.8g.', ...
        globalIx, beamIx, energyValue);
end

energyIndex = energyIx;
end
