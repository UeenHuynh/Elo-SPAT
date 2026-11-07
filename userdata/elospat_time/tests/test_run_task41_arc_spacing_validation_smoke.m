clear; clc;

rootDir = fileparts(mfilename('fullpath'));
workspaceDir = fileparts(rootDir);
addpath(fullfile(rootDir, 'mocks'));
addpath(fullfile(workspaceDir, 'planning'));
addpath(workspaceDir);

fprintf('\n========================================\n');
fprintf('TEST: run_task41_arc_spacing_validation smoke\n');
fprintf('========================================\n');

cfg = build_precomputed_smoke_fixture(fullfile(tempdir, 'task41_smoke'));
outputs = run_task41_arc_spacing_validation([20, 10], cfg, fullfile(tempdir, 'task41_smoke_outputs'));

assert(height(outputs.summaryTable) == 2);
assert(all(strcmp(outputs.summaryTable.status, "completed")));
assert(exist(outputs.files.csvFile, 'file') == 2);
assert(exist(outputs.files.matFile, 'file') == 2);
assert(exist(outputs.files.plotFile, 'file') == 2);

fprintf('[PASS] Task 41 smoke completed and saved spacing summary artifacts.\n');
