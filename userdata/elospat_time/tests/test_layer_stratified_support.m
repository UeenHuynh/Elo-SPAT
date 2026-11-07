% test_layer_stratified_support  Layer coverage and nested-support checks.
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'integrated'));

weights = [9; 8; 1; 7; 6; 2; 5; 4; 3; 8; 7; 1];
cp = [ones(6,1); 2*ones(6,1)];
energyIndex = repmat([1;1;2;2;3;3], 2, 1);
map = struct('cpIndex', num2cell(cp), ...
    'energyIndex', num2cell(energyIndex));

[small, dSmall] = build_layer_stratified_support(weights, map, 3);
[large, dLarge] = build_layer_stratified_support(weights, map, 5);
assert(all(small <= large), 'Increasing the cap must produce a nested support.');
assert(nnz(small) == 6 && nnz(large) == 10, 'Unexpected support sizes.');
assert(all([dSmall.selectedLayers] == 3), ...
    'The small support must retain every available layer.');
assert(all([dLarge.selectedLayers] == 3), ...
    'The large support must retain every available layer.');
assert(all([dSmall.selectedSpots] == 3));
assert(all([dLarge.selectedSpots] == 5));

fprintf('[OK] Layer-stratified support is layer-complete and nested.\n');
