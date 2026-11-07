% test_mask_dimension_safety  Ensure layerMask cannot be mistaken for bixelMask.

clear; clc;

reductionDir = fullfile(fileparts(mfilename('fullpath')), '..', 'reduction');
addpath(reductionDir);

fprintf('\n========================================\n');
fprintf('TEST: mask dimension safety (Task 21)\n');
fprintf('========================================\n\n');

groups = repmat(struct('ix', []), 2, 1);
groups(1).ix = [1; 2];
groups(2).ix = [3; 4];

wrongMask = [true; false; true; false];
caughtExpectedError = false;

try
    expand_layer_mask_to_bixel_mask(wrongMask, groups, 4);
catch ME
    caughtExpectedError = strcmp(ME.identifier, ...
        'expand_layer_mask_to_bixel_mask:dimensionMismatch');
    if ~caughtExpectedError
        rethrow(ME);
    end
end

assert(caughtExpectedError, ...
    'Expected a dimensionMismatch error when passing a bixel-sized mask.');

fprintf('[OK] Dimension safety check validated.\n');

fprintf('\n========================================\n');
fprintf('TEST PASSED\n');
fprintf('========================================\n\n');
