% test_build_masks  Validation for Task 20 layer and spot masks.

clear; clc;

reductionDir = fullfile(fileparts(mfilename('fullpath')), '..', 'reduction');
addpath(reductionDir);

fprintf('\n========================================\n');
fprintf('TEST: build_layer_mask / build_spot_mask (Task 20)\n');
fprintf('========================================\n\n');

w = [0.10; 0.20; 0.01; 0.00; 0.50; 0.02];
groups = repmat(struct('ix', []), 3, 1);
groups(1).ix = [1; 2];
groups(2).ix = [3; 4];
groups(3).ix = [5; 6];

layerThreshold = 0.15;
spotThreshold = 0.05;

layerMask = build_layer_mask(w, groups, layerThreshold);
spotMask = build_spot_mask(w, spotThreshold);

assert(islogical(layerMask) && isequal(size(layerMask), [3, 1]), ...
    'layerMask must be a logical column vector of size nGroups x 1.');
assert(islogical(spotMask) && isequal(size(spotMask), [6, 1]), ...
    'spotMask must be a logical column vector of size N x 1.');

expectedLayerMask = [true; false; true];
expectedSpotMask = [true; true; false; false; true; false];

assert(isequal(layerMask, expectedLayerMask), ...
    'Unexpected layerMask values.');
assert(isequal(spotMask, expectedSpotMask), ...
    'Unexpected spotMask values.');

fprintf('[OK] Layer and spot masks validated.\n');
fprintf('     layerMask = [%s]\n', num2str(layerMask'));
fprintf('     spotMask  = [%s]\n', num2str(spotMask'));

fprintf('\n========================================\n');
fprintf('TEST PASSED\n');
fprintf('========================================\n\n');
