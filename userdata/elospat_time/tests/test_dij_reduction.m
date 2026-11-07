% test_dij_reduction  Validation for Task 22 multi-scenario dij reduction.

clear; clc;

reductionDir = fullfile(fileparts(mfilename('fullpath')), '..', 'reduction');
addpath(reductionDir);

fprintf('\n========================================\n');
fprintf('TEST: cut_dij_by_mask (Task 22)\n');
fprintf('========================================\n\n');

dij = struct();
dij.totalNumOfBixels = 5;
dij.physicalDose = {sparse([1 2 3 4 5; 6 7 8 9 10]), sparse([10 20 30 40 50; 60 70 80 90 100])};
dij.mAlphaDose = {sparse([11 12 13 14 15; 16 17 18 19 20]), sparse([21 22 23 24 25; 26 27 28 29 30])};
dij.mSqrtBetaDose = {sparse([31 32 33 34 35; 36 37 38 39 40]), sparse([41 42 43 44 45; 46 47 48 49 50])};
dij.mLETDose = {sparse([51 52 53 54 55; 56 57 58 59 60]), sparse([61 62 63 64 65; 66 67 68 69 70])};
dij.beamNum = [1 1 2 2 3];
dij.rayNum = [1 2 1 2 1];
dij.bixelNum = [1 1 2 2 3];
dij.minMU = [0.1 0.2 0.3 0.4 0.5];
dij.maxMU = [1.1 1.2 1.3 1.4 1.5];
dij.numParticlesPerMU = [101 102 103 104 105];
dij.unrelatedField = 'keep_me';

map = struct('globalIx', num2cell((1:5)'));
w = [0.9; 0.8; 0.7; 0.6; 0.5];
mask = [true; false; true; false; true];

[dijReduced, mapReduced, wReduced] = cut_dij_by_mask(dij, mask, map, w);

assert(dijReduced.totalNumOfBixels == 3, 'Reduced dij must have 3 bixels.');
assert(isequal(size(dijReduced.physicalDose{1}), [2, 3]), ...
    'Reduced physicalDose has incorrect size.');
assert(isequal(full(dijReduced.physicalDose{1}), [1 3 5; 6 8 10]), ...
    'Reduced physicalDose scenario 1 is incorrect.');
assert(isequal(full(dijReduced.physicalDose{2}), [10 30 50; 60 80 100]), ...
    'Reduced physicalDose scenario 2 is incorrect.');
assert(isequal(dijReduced.beamNum, [1 2 3]), 'Reduced beamNum is incorrect.');
assert(isequal(dijReduced.rayNum, [1 1 1]), 'Reduced rayNum is incorrect.');
assert(isequal(dijReduced.bixelNum, [1 2 3]), 'Reduced bixelNum is incorrect.');
assert(isequal(dijReduced.minMU, [0.1 0.3 0.5]), 'Reduced minMU is incorrect.');
assert(isequal(dijReduced.maxMU, [1.1 1.3 1.5]), 'Reduced maxMU is incorrect.');
assert(isequal(dijReduced.numParticlesPerMU, [101 103 105]), ...
    'Reduced numParticlesPerMU is incorrect.');
assert(strcmp(dijReduced.unrelatedField, 'keep_me'), ...
    'Unrelated dij fields should be preserved unchanged.');
assert(isequal([mapReduced.globalIx]', [1; 3; 5]), ...
    'mapReduced must preserve original globalIx values.');
assert(isequal(wReduced, [0.9; 0.7; 0.5]), 'wReduced is incorrect.');
assert(isequal(dijReduced.keptBixelIndices, [1; 3; 5]), ...
    'keptBixelIndices must match the retained bixels.');

fprintf('[OK] Multi-scenario dij reduction validated.\n');

fprintf('\n========================================\n');
fprintf('TEST PASSED\n');
fprintf('========================================\n\n');
