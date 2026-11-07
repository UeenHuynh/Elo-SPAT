% test_validate_mapping  Focused unit test for Task 10 mapping validation utilities.

clear; clc;

mappingDir = fullfile(fileparts(mfilename('fullpath')), '..', 'mapping');
addpath(mappingDir);

fprintf('\n========================================\n');
fprintf('TEST: validate_mapping (Task 10)\n');
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

dij = struct('totalNumOfBixels', 6);
[groups, cpGroups] = build_group_index(map);

validate_mapping(map, dij, groups, cpGroups);
fprintf('[OK] Happy-path validation passed.\n');

groupIx = get_group_indices(groups, 2);
cpIx = get_control_point_indices(cpGroups, 2);

assert(isequal(groupIx, [2; 3]), 'get_group_indices returned wrong indices.');
assert(isequal(cpIx, [4; 5; 6]), 'get_control_point_indices returned wrong indices.');
fprintf('[OK] Accessor helpers validated.\n');

badGroups = groups;
badGroups(2).ix = [2; 2];

didFail = false;
try
    validate_mapping(map, dij, badGroups, cpGroups);
catch ME
    didFail = strcmp(ME.identifier, 'validate_mapping:groupCoverageMismatch');
    if ~didFail
        rethrow(ME);
    end
end

assert(didFail, 'Expected duplicate group coverage to be detected.');
fprintf('[OK] Duplicate group coverage correctly rejected.\n');

badMap = map;
badMap(6).globalIx = 7;

didFail = false;
try
    validate_mapping(badMap, dij, groups, cpGroups);
catch ME
    didFail = strcmp(ME.identifier, 'validate_mapping:globalIxMismatch');
    if ~didFail
        rethrow(ME);
    end
end

assert(didFail, 'Expected non-consecutive globalIx to be detected.');
fprintf('[OK] Non-consecutive globalIx correctly rejected.\n');

fprintf('\n========================================\n');
fprintf('TEST PASSED\n');
fprintf('========================================\n\n');
