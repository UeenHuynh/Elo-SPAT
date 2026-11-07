clear; clc;

rootDir = fileparts(mfilename('fullpath'));
workspaceDir = fileparts(rootDir);
addpath(fullfile(workspaceDir, 'planning'));
addpath(workspaceDir);

fprintf('\n========================================\n');
fprintf('TEST: run_task42_documentation_check\n');
fprintf('========================================\n');

outputs = run_task42_documentation_check(fullfile(tempdir, 'task42_doc_check'));

assert(outputs.report.checks.readmeHasResearchOnly);
assert(outputs.report.checks.readmeHasProxyNote);
assert(outputs.report.checks.readmeHasUsageGuide);
assert(outputs.report.checks.readmeHasCitationPlaceholder);
assert(outputs.report.checks.latexHasDisclaimer);
assert(outputs.report.checks.latexHasMachineProxyNote);
assert(outputs.report.checks.latexHasCitationPlaceholder);
assert(exist(outputs.files.matFile, 'file') == 2);
assert(exist(outputs.files.txtFile, 'file') == 2);

fprintf('[PASS] Task 42 documentation check passed.\n');
