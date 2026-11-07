% test_extract_delivery_sequence  Validation for Task 16 delivery ordering.

clear; clc;

deliveryDir = fullfile(fileparts(mfilename('fullpath')), '..', 'delivery');
addpath(deliveryDir);

fprintf('\n========================================\n');
fprintf('TEST: extract_delivery_sequence (Task 16)\n');
fprintf('========================================\n\n');

w = [0.90; 0.20; 0.80; 0.40; 0.70; 0.10];

map = repmat(struct('cpIndex', 0, 'gantryAngle', 0), 6, 1);
map(1).cpIndex = 2; map(1).gantryAngle = 40;
map(2).cpIndex = 2; map(2).gantryAngle = 40;
map(3).cpIndex = 2; map(3).gantryAngle = 40;
map(4).cpIndex = 1; map(4).gantryAngle = 20;
map(5).cpIndex = 1; map(5).gantryAngle = 20;
map(6).cpIndex = 1; map(6).gantryAngle = 20;

groups = repmat(struct('ix', [], 'cpIndex', 0, 'energy', 0, 'energyIndex', 0), 4, 1);
groups(1).ix = [1; 2]; groups(1).cpIndex = 2; groups(1).energy = 150; groups(1).energyIndex = 2;
groups(2).ix = 3;      groups(2).cpIndex = 2; groups(2).energy = 180; groups(2).energyIndex = 1;
groups(3).ix = 4;      groups(3).cpIndex = 1; groups(3).energy = 160; groups(3).energyIndex = 1;
groups(4).ix = [5; 6]; groups(4).cpIndex = 1; groups(4).energy = 120; groups(4).energyIndex = 2;

cfg = struct('layerThreshold', 0.15, 'spotThreshold', 0.15);

seq = extract_delivery_sequence(w, map, groups, cfg);

assert(numel(seq) == 4, 'Expected 4 active layers in the delivery sequence.');
assert(isequal([seq.cpIndex]', [1; 1; 2; 2]), ...
    'Control points must be ordered by increasing gantry angle.');
assert(isequal([seq.energy]', [160; 120; 180; 150]), ...
    'Layers within each CP must be ordered by descending energy by default.');

assert(isequal(seq(1).spotIndices, 4), 'Unexpected active spots for first layer.');
assert(isequal(seq(2).spotIndices, 5), ...
    'Spot thresholding should remove low-weight spots within a layer.');
assert(seq(2).numActiveSpots == 1, 'Expected one active spot after thresholding.');
assert(seq(3).prevEnergy == 120 && seq(3).nextEnergy == 150, ...
    'prevEnergy/nextEnergy must follow the final delivered sequence.');
assert(isempty(seq(1).prevEnergy) && isempty(seq(end).nextEnergy), ...
    'Boundary sequence entries must have empty prev/next energy.');

fprintf('[OK] Delivery sequence ordering validated.\n');
fprintf('     Sequence energies: [%s]\n', num2str([seq.energy]'));
fprintf('     Sequence CPs:      [%s]\n', num2str([seq.cpIndex]'));

fprintf('\n========================================\n');
fprintf('TEST PASSED\n');
fprintf('========================================\n\n');
