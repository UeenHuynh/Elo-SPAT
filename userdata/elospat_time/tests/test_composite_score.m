% test_composite_score  Validation for Task 15 composite_score.

clear; clc;

objectiveDir = fullfile(fileparts(mfilename('fullpath')), '..', 'objective');
addpath(objectiveDir);

fprintf('\n========================================\n');
fprintf('TEST: composite_score (Task 15)\n');
fprintf('========================================\n\n');

w = [0.8; 0.1; 0.6; 0.3; 1.2; 0.5];

groups = repmat(struct('ix', [], 'cpIndex', 0, 'energy', 0), 4, 1);
groups(1).ix = [1; 2]; groups(1).cpIndex = 1; groups(1).energy = 200;
groups(2).ix = 3;      groups(2).cpIndex = 1; groups(2).energy = 150;
groups(3).ix = [4; 5]; groups(3).cpIndex = 2; groups(3).energy = 210;
groups(4).ix = 6;      groups(4).cpIndex = 2; groups(4).energy = 140;

cpGroups = repmat(struct('ix', [], 'gantryAngle', 0, 'cpIndex', 0), 2, 1);
cpGroups(1).ix = [1; 2; 3]; cpGroups(1).gantryAngle = 0;  cpGroups(1).cpIndex = 1;
cpGroups(2).ix = [4; 5; 6]; cpGroups(2).gantryAngle = 20; cpGroups(2).cpIndex = 2;

cfg = struct( ...
    'lambdaEL', 0, ...
    'lambdaB', 0, ...
    'lambdaES', 0, ...
    'lambdaS', 0, ...
    'lambdaT', 0, ...
    'epsEL', 1e-8, ...
    'epsB', 1e-8, ...
    'epsZ', 1e-6, ...
    'epsS', 1e-3);

[zeroScore, zeroComponents] = composite_score(w, groups, cpGroups, cfg);
assert(abs(zeroScore) < 1e-12, 'All-zero lambdas must give score 0.');
assert(isstruct(zeroComponents), 'components must be returned as a struct.');
fprintf('[OK] All-zero lambda case validated.\n');

cfg.lambdaEL = 1.1;
cfg.lambdaB = 0.7;
cfg.lambdaES = 0.05;
cfg.lambdaS = 0.2;
cfg.lambdaT = 0.01;

bdtScalar = 12.5;
[score, components] = composite_score(w, groups, cpGroups, cfg, bdtScalar);

expectedScore = components.weightedL21half + ...
    components.weightedBarrier + ...
    components.weightedSequencing + ...
    components.weightedSpotSparsity + ...
    components.weightedBdt;

assert(abs(score - expectedScore) < 1e-12, ...
    'Composite score must equal the sum of weighted components.');
assert(abs(components.bdt - bdtScalar) < 1e-12, ...
    'BDT component must match the supplied scalar.');

fprintf('[OK] Weighted sum case validated.\n');
fprintf('     Composite score: %.6f\n', score);
fprintf('     Hard switch-up count: %d\n', components.switchUpCount);

fprintf('\n========================================\n');
fprintf('TEST PASSED\n');
fprintf('========================================\n\n');
