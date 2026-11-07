% test_bdt_estimator  Validation for Tasks 18-19 BDT and gantry helpers.

clear; clc;

deliveryDir = fullfile(fileparts(mfilename('fullpath')), '..', 'delivery');
configDir = fullfile(fileparts(mfilename('fullpath')), '..', 'config');
addpath(deliveryDir);
addpath(configDir);

fprintf('\n========================================\n');
fprintf('TEST: BDT estimator and gantry helpers (Task 19)\n');
fprintf('========================================\n\n');

machine = default_machine_model();

tSmall = estimate_gantry_time(9, machine);
tLarge = estimate_gantry_time(40, machine);

expectedSmall = 2 * sqrt(9 / machine.gantry.maxAccel);
expectedLarge = 40 / machine.gantry.maxSpeed + machine.gantry.maxSpeed / machine.gantry.maxAccel;

assert(abs(tSmall - expectedSmall) < 1e-12, ...
    'Small-angle triangle-profile time is incorrect.');
assert(abs(tLarge - expectedLarge) < 1e-12, ...
    'Large-angle trapezoid-profile time is incorrect.');
fprintf('[OK] estimate_gantry_time validated on triangle/trapezoid cases.\n');

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

metrics = estimate_bdt(seq, machine);
[feasible, violations] = validate_gantry_constraints(seq, machine);

componentSum = metrics.T_energy + metrics.T_spot + metrics.T_irradiation + ...
    metrics.T_gantry + metrics.T_buffer;

assert(abs(metrics.T_total - componentSum) < 1e-10, ...
    'T_total must equal the sum of all BDT components.');
assert(feasible, 'Expected the nominal test sequence to satisfy gantry constraints.');
assert(isempty(violations), 'Expected no gantry constraint violations.');
fprintf('[OK] estimate_bdt and validate_gantry_constraints validated together.\n');

fprintf('\n========================================\n');
fprintf('TEST PASSED\n');
fprintf('========================================\n\n');
