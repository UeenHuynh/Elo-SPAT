clear; clc;

rootDir = fileparts(mfilename('fullpath'));
workspaceDir = fileparts(rootDir);
addpath(fullfile(rootDir, 'mocks'));
addpath(fullfile(workspaceDir, 'planning'));
addpath(workspaceDir);

fprintf('\n========================================\n');
fprintf('TEST: run_task38_integration smoke\n');
fprintf('========================================\n');

cfg = build_precomputed_smoke_fixture(fullfile(tempdir, 'task38_smoke'));
outputs = run_task38_integration(cfg, cfg.outputDir);

assert(outputs.report.verification.allMethodsPresent);
assert(outputs.report.verification.ctvD95WithinTolerance);
assert(exist(outputs.reportFile, 'file') == 2);
assert(isfield(outputs.report.verification, 'rmsBdtBelowDense'));

fprintf('[PASS] Task 38 smoke runner produced a verified report.\n');
