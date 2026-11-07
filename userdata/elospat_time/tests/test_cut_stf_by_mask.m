% test_cut_stf_by_mask  Validation for Task 23 STF reduction consistency.

clear; clc;

reductionDir = fullfile(fileparts(mfilename('fullpath')), '..', 'reduction');
addpath(reductionDir);

fprintf('\n========================================\n');
fprintf('TEST: cut_stf_by_mask (Task 23)\n');
fprintf('========================================\n\n');

stf = repmat(struct(), 2, 1);
stf(1).gantryAngle = 0;
stf(1).numOfRays = 2;
stf(1).ray = repmat(struct(), 2, 1);
stf(1).ray(1).energy = [150 140];
stf(1).ray(1).focusIx = [10 11];
stf(1).ray(2).energy = [130];
stf(1).ray(2).focusIx = [12];

stf(2).gantryAngle = 20;
stf(2).numOfRays = 1;
stf(2).ray = struct();
stf(2).ray(1).energy = [160 150 140];
stf(2).ray(1).focusIx = [21 22 23];

map = struct( ...
    'globalIx', num2cell((1:6)'), ...
    'beamNum', num2cell([1; 1; 1; 2; 2; 2]), ...
    'rayNum', num2cell([1; 1; 2; 1; 1; 1]), ...
    'bixelNum', num2cell([1; 2; 1; 1; 2; 3]));

mask = [true; false; false; true; false; true];
stfReduced = cut_stf_by_mask(stf, mask, map);

assert(numel(stfReduced) == 2, 'Expected both beams to remain.');
assert(stfReduced(1).numOfRays == 1, 'Beam 1 should retain one non-empty ray.');
assert(numel(stfReduced(1).ray) == 1, 'Beam 1 should have one ray after reduction.');
assert(isequal(stfReduced(1).ray(1).energy, 150), 'Beam 1 ray 1 energy is incorrect.');
assert(isequal(stfReduced(1).ray(1).focusIx, 10), 'Beam 1 ray 1 focusIx is incorrect.');

assert(stfReduced(2).numOfRays == 1, 'Beam 2 should retain one ray.');
assert(isequal(stfReduced(2).ray(1).energy, [160 140]), 'Beam 2 energies are incorrect.');
assert(isequal(stfReduced(2).ray(1).focusIx, [21 23]), 'Beam 2 focusIx is incorrect.');

fprintf('[OK] STF reduction validated.\n');

fprintf('\n========================================\n');
fprintf('TEST PASSED\n');
fprintf('========================================\n\n');
