% test_estimate_bdt_basic  Basic validation for Task 18 BDT decomposition.

clear; clc;

deliveryDir = fullfile(fileparts(mfilename('fullpath')), '..', 'delivery');
configDir = fullfile(fileparts(mfilename('fullpath')), '..', 'config');
addpath(deliveryDir);
addpath(configDir);

fprintf('\n========================================\n');
fprintf('TEST: estimate_bdt basic decomposition (Task 18)\n');
fprintf('========================================\n\n');

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
    'nextEnergy', []), 3, 1);

seq(1).cpIndex = 1; seq(1).gantryAngle = 0;  seq(1).energy = 200; seq(1).numActiveSpots = 3; seq(1).layerWeightSum = 2.0;
seq(2).cpIndex = 1; seq(2).gantryAngle = 0;  seq(2).energy = 180; seq(2).numActiveSpots = 2; seq(2).layerWeightSum = 1.0;
seq(3).cpIndex = 2; seq(3).gantryAngle = 20; seq(3).energy = 190; seq(3).numActiveSpots = 1; seq(3).layerWeightSum = 4.0;

machine = default_machine_model();
metrics = estimate_bdt(seq, machine);

componentSum = metrics.T_energy + metrics.T_spot + metrics.T_irradiation + ...
    metrics.T_gantry + metrics.T_buffer;

assert(abs(metrics.T_total - componentSum) < 1e-10, ...
    'T_total must equal the sum of all BDT components.');
assert(metrics.nActiveLayers == 3, 'Expected 3 active layers.');
assert(metrics.nActiveSpots == 6, 'Expected 6 active spots in total.');
assert(metrics.nSwitchUp == 1, 'Expected one low-to-high switch-up event.');
assert(metrics.totalSwitchMagnitude == 30, ...
    'Expected total absolute switch magnitude of 30 MeV.');

fprintf('[OK] BDT decomposition validated.\n');
fprintf('     T_total: %.6f s\n', metrics.T_total);
fprintf('     nSwitchUp: %d\n', metrics.nSwitchUp);

fprintf('\n========================================\n');
fprintf('TEST PASSED\n');
fprintf('========================================\n\n');
