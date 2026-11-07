% test_build_group_index  Focused unit test for Task 9 grouping utility.

clear; clc;

mappingDir = fullfile(fileparts(mfilename('fullpath')), '..', 'mapping');
addpath(mappingDir);

fprintf('\n========================================\n');
fprintf('TEST: build_group_index (Task 9)\n');
fprintf('========================================\n\n');

map = struct( ...
    'globalIx', num2cell((1:6)'), ...
    'beamNum', num2cell([1; 1; 1; 2; 2; 2]), ...
    'rayNum', num2cell([1; 1; 2; 1; 1; 1]), ...
    'bixelNum', num2cell([1; 2; 1; 1; 2; 3]), ...
    'gantryAngle', num2cell([0; 0; 0; 30; 30; 30]), ...
    'energy', num2cell([150; 140; 140; 160; 150; 150]), ...
    'cpIndex', num2cell([1; 1; 1; 2; 2; 2]), ...
    'energyIndex', num2cell([1; 2; 2; 1; 2; 2]));

[groups, cpGroups] = build_group_index(map);

fprintf('[OK] Grouping function returned successfully.\n');
assert(isstruct(groups) && numel(groups) == 4, ...
    'Expected 4 energy-layer groups.');
assert(isstruct(cpGroups) && numel(cpGroups) == 2, ...
    'Expected 2 control-point groups.');
fprintf('[OK] Correct number of groups and control-point groups.\n');

assert(isequal(groups(1).ix, [1]), 'Group 1 indices incorrect.');
assert(isequal(groups(2).ix, [2; 3]), 'Group 2 indices incorrect.');
assert(isequal(groups(3).ix, [4]), 'Group 3 indices incorrect.');
assert(isequal(groups(4).ix, [5; 6]), 'Group 4 indices incorrect.');

assert(isequal([groups.nSpots]', [1; 2; 1; 2]), 'Group spot counts incorrect.');
assert(isequal([groups.cpIndex]', [1; 1; 2; 2]), 'Group cpIndex values incorrect.');
assert(isequal([groups.energyIndex]', [1; 2; 1; 2]), 'Group energyIndex values incorrect.');
assert(isequal([groups.energy]', [150; 140; 160; 150]), 'Group energy values incorrect.');

assert(isequal(cpGroups(1).ix, [1; 2; 3]), 'cpGroups(1).ix incorrect.');
assert(isequal(cpGroups(2).ix, [4; 5; 6]), 'cpGroups(2).ix incorrect.');
assert(isequal([cpGroups.cpIndex]', [1; 2]), 'cpGroups cpIndex incorrect.');
assert(isequal([cpGroups.gantryAngle]', [0; 30]), 'cpGroups gantryAngle incorrect.');

fprintf('[OK] groups and cpGroups contents validated.\n');

fprintf('\n========================================\n');
fprintf('TEST PASSED\n');
fprintf('========================================\n\n');
