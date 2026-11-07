% test_reduce_dij_to_support  Regression test for support-preserving reduction.

clear; clc;
rootDir = fileparts(mfilename('fullpath'));
addpath(fullfile(rootDir, '..', 'reduction'));

dij = struct();
dij.totalNumOfBixels = 5;
dij.physicalDose = {reshape(1:15, 3, 5)};
dij.mAlphaDose = {reshape(16:30, 3, 5)};
dij.beamNum = (11:15)';
dij.rayNum = (21:25)';
dij.bixelNum = (31:35)';

supportMask = logical([1; 0; 1; 0; 1]);
[reduced, toFull, toReduced] = reduce_dij_to_support(dij, supportMask);

assert(isequal(toFull, [1; 3; 5]));
assert(isequal(toReduced, [1; 0; 2; 0; 3]));
assert(isequal(reduced.physicalDose{1}, dij.physicalDose{1}(:, supportMask)));
assert(isequal(reduced.mAlphaDose{1}, dij.mAlphaDose{1}(:, supportMask)));
assert(isequal(reduced.beamNum, dij.beamNum(supportMask)));

xReduced = [0.7; 1.2; 0.4];
xFull = zeros(dij.totalNumOfBixels, 1);
xFull(toFull) = xReduced;
assert(norm(dij.physicalDose{1} * xFull - reduced.physicalDose{1} * xReduced) < 1e-12);

fprintf('[OK] reduce_dij_to_support regression test passed.\n');
