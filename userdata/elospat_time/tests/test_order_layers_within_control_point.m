% test_order_layers_within_control_point  Validation for Task 17 layer ordering.

clear; clc;

deliveryDir = fullfile(fileparts(mfilename('fullpath')), '..', 'delivery');
addpath(deliveryDir);

fprintf('\n========================================\n');
fprintf('TEST: order_layers_within_control_point (Task 17)\n');
fprintf('========================================\n\n');

energies = [150; 200; 120; 180];

orderedDescending = order_layers_within_control_point(energies, 'descending');
assert(isequal(orderedDescending, [200; 180; 150; 120]), ...
    'Descending order is incorrect.');
fprintf('[OK] Descending order validated.\n');

orderedAscending = order_layers_within_control_point(energies, 'ascending');
assert(isequal(orderedAscending, [120; 150; 180; 200]), ...
    'Ascending order is incorrect.');
fprintf('[OK] Ascending order validated.\n');

orderedCustom = order_layers_within_control_point(energies, 'custom', [180; 120; 200; 150]);
assert(isequal(orderedCustom, [180; 120; 200; 150]), ...
    'Custom order is incorrect.');
fprintf('[OK] Custom order validated.\n');

fprintf('\n========================================\n');
fprintf('TEST PASSED\n');
fprintf('========================================\n\n');
