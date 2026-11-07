% test_expand_layer_mask  Validation for Task 21 mask expansion.

clear; clc;

reductionDir = fullfile(fileparts(mfilename('fullpath')), '..', 'reduction');
addpath(reductionDir);

fprintf('\n========================================\n');
fprintf('TEST: expand_layer_mask_to_bixel_mask (Task 21)\n');
fprintf('========================================\n\n');

groups = repmat(struct('ix', []), 3, 1);
groups(1).ix = [1; 2];
groups(2).ix = [3; 5];
groups(3).ix = [4; 6; 7];

layerMask = [true; false; true];
bixelMask = expand_layer_mask_to_bixel_mask(layerMask, groups, 7);

expectedMask = [true; true; false; true; false; true; true];

assert(islogical(bixelMask) && isequal(size(bixelMask), [7, 1]), ...
    'bixelMask must be a logical column vector of length totalNumOfBixels.');
assert(isequal(bixelMask, expectedMask), ...
    'Expanded bixelMask does not match expected membership.');

fprintf('[OK] Layer-mask expansion validated.\n');
fprintf('     bixelMask = [%s]\n', num2str(bixelMask'));

fprintf('\n========================================\n');
fprintf('TEST PASSED\n');
fprintf('========================================\n\n');
