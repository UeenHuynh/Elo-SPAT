clear; clc;

rootDir = fileparts(mfilename('fullpath'));
workspaceDir = fileparts(rootDir);
addpath(fullfile(rootDir, 'mocks'));
addpath(fullfile(workspaceDir, 'planning'));
addpath(workspaceDir);

fprintf('\n========================================\n');
fprintf('TEST: run_task40_parameter_sweep smoke\n');
fprintf('========================================\n');

cfg = build_precomputed_smoke_fixture(fullfile(tempdir, 'task40_smoke'));
cfg.randomSeed = 1;

sweepOverride = struct();
sweepOverride.parameters = {'arcSpacing', 'wrapper.maxLayerIterations'};
sweepOverride.values = {[20; 30], [1; 2]};
sweepOverride.outputDir = fullfile(tempdir, 'task40_smoke_outputs');

outputs = run_task40_parameter_sweep(sweepOverride, cfg, sweepOverride.outputDir);

assert(isfield(outputs, 'summaryTable'));
assert(height(outputs.summaryTable) == 4);
assert(exist(outputs.files.summaryMatFile, 'file') == 2);
assert(exist(outputs.files.summaryCsvFile, 'file') == 2);
assert(exist(outputs.files.heatmapFile, 'file') == 2);

fprintf('[PASS] Task 40 sweep smoke completed and saved aggregate artifacts.\n');
