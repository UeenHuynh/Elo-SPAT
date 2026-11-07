% test_map_bixel_to_cp_energy_spot  Focused unit test for Task 8 mapping utility.

clear; clc;

mappingDir = fullfile(fileparts(mfilename('fullpath')), '..', 'mapping');
addpath(mappingDir);

fprintf('\n========================================\n');
fprintf('TEST: map_bixel_to_cp_energy_spot (Task 8)\n');
fprintf('========================================\n\n');

% -------------------------------------------------------------------------
% Synthetic inputs
% -------------------------------------------------------------------------
pln = struct();
pln.propStf = struct();
pln.propStf.gantryAngles = [0, 30];

stf = repmat(struct(), 2, 1);
stf(1).gantryAngle = 0;
stf(1).ray = repmat(struct(), 2, 1);
stf(1).ray(1).energy = [150, 140];
stf(1).ray(2).energy = [140];

stf(2).gantryAngle = 30;
stf(2).ray = struct();
stf(2).ray(1).energy = [160, 150, 150];

dij = struct();
dij.totalNumOfBixels = 6;
dij.beamNum = [1; 1; 1; 2; 2; 2];
dij.rayNum = [1; 1; 2; 1; 1; 1];
dij.bixelNum = [1; 2; 1; 1; 2; 3];

% -------------------------------------------------------------------------
% Run mapping
% -------------------------------------------------------------------------
map = map_bixel_to_cp_energy_spot(dij, stf, pln);

fprintf('[OK] Mapping function returned successfully.\n');
assert(isstruct(map), 'map must be a struct array.');
assert(numel(map) == dij.totalNumOfBixels, ...
    'map length must equal dij.totalNumOfBixels.');
fprintf('[OK] Output length matches dij.totalNumOfBixels (%d).\n', numel(map));

expectedEnergy = [150; 140; 140; 160; 150; 150];
expectedEnergyIndex = [1; 2; 2; 1; 2; 2];
expectedGantryAngle = [0; 0; 0; 30; 30; 30];

assert(isequal(reshape([map.globalIx], [], 1), (1:6)'), ...
    'globalIx must be consecutive from 1:N.');
assert(isequal(reshape([map.cpIndex], [], 1), dij.beamNum), ...
    'cpIndex must match beamNum for the discrete arc surrogate.');
assert(isequal(reshape([map.energy], [], 1), expectedEnergy), ...
    'Mapped energies do not match stf lookup.');
assert(isequal(reshape([map.energyIndex], [], 1), expectedEnergyIndex), ...
    'energyIndex assignment is incorrect.');
assert(isequal(reshape([map.gantryAngle], [], 1), expectedGantryAngle), ...
    'gantryAngle assignment is incorrect.');

fprintf('[OK] globalIx, cpIndex, energy, energyIndex, and gantryAngle validated.\n');
fprintf('[OK] Beam 1 energy layers sorted descending: [150, 140] -> indices [1, 2]\n');
fprintf('[OK] Beam 2 energy layers sorted descending: [160, 150] -> indices [1, 2]\n');

fprintf('\n========================================\n');
fprintf('TEST PASSED\n');
fprintf('========================================\n\n');
