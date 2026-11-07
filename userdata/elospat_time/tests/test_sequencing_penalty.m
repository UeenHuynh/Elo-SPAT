% test_sequencing_penalty  Validation for Task 13 energy sequencing penalty.

clear; clc;

objectiveDir = fullfile(fileparts(mfilename('fullpath')), '..', 'objective');
addpath(objectiveDir);

fprintf('\n========================================\n');
fprintf('TEST: energy_sequencing_penalty (Task 13)\n');
fprintf('========================================\n\n');

groups = repmat(struct('ix', [], 'cpIndex', 0, 'energy', 0), 6, 1);
groups(1).ix = 1; groups(1).cpIndex = 1; groups(1).energy = 200;
groups(2).ix = 2; groups(2).cpIndex = 1; groups(2).energy = 150;
groups(3).ix = 3; groups(3).cpIndex = 2; groups(3).energy = 180;
groups(4).ix = 4; groups(4).cpIndex = 2; groups(4).energy = 140;
groups(5).ix = 5; groups(5).cpIndex = 3; groups(5).energy = 160;
groups(6).ix = 6; groups(6).cpIndex = 3; groups(6).energy = 120;

epsZ = 1e-6;

% Descending active sequence: 200 -> 180 -> 160
wDescending = [1; 0; 1; 0; 1; 0];
[valDescending, countDescending] = energy_sequencing_penalty( ...
    wDescending, groups, epsZ);

assert(abs(valDescending) < 1e-12, ...
    'Descending active sequence should have zero smooth penalty.');
assert(countDescending == 0, ...
    'Descending active sequence should have zero hard switch-up count.');
fprintf('[OK] Monotonic descending sequence gives zero penalty.\n');

% Ascending active sequence: 150 -> 180 -> 160 (one switch-up then one down)
wAscending = [0; 1; 1; 0; 1; 0];
[valAscending, countAscending] = energy_sequencing_penalty( ...
    wAscending, groups, epsZ);

assert(valAscending > 0, ...
    'Ascending sequence should produce a positive smooth penalty.');
assert(countAscending >= 1, ...
    'Ascending sequence should produce at least one hard switch-up event.');
fprintf('[OK] Low-to-high transition gives positive penalty.\n');
fprintf('     Smooth penalty: %.6f\n', valAscending);
fprintf('     Hard switch-up count: %d\n', countAscending);

fprintf('\n========================================\n');
fprintf('TEST PASSED\n');
fprintf('========================================\n\n');
