function outputs = run_task42_documentation_check(outputDir)
% run_task42_documentation_check  Verify Task 42 documentation/disclaimer items.

narginchk(0, 1);

planningDir = fileparts(mfilename('fullpath'));
workspaceDir = fileparts(planningDir);
evaluationDir = fullfile(workspaceDir, 'evaluation');
resultsDir = fullfile(workspaceDir, 'results');
readmeFile = fullfile(workspaceDir, 'README.md');

addpath(workspaceDir);
addpath(evaluationDir);

if nargin < 1 || isempty(outputDir)
    outputDir = resultsDir;
end
if exist(outputDir, 'dir') ~= 7
    mkdir(outputDir);
end

readmeText = fileread(readmeFile);
researchMeta = research_only_metadata();

checks = struct();
checks.readmeHasResearchOnly = contains(readmeText, 'Research-Only Disclaimer', 'IgnoreCase', true);
checks.readmeHasProxyNote = contains(readmeText, 'estimated configurable proxy', 'IgnoreCase', true) || ...
    contains(readmeText, 'estimated proxy', 'IgnoreCase', true);
checks.readmeHasUsageGuide = contains(readmeText, '## Usage', 'IgnoreCase', true);
checks.readmeHasCitationPlaceholder = contains(readmeText, 'placeholder', 'IgnoreCase', true) || ...
    contains(readmeText, 'TODO: add manuscript', 'IgnoreCase', true);

mockTable = table( ...
    ["Dense Arc"; "RMS-ELO-SPAT"], ...
    [70; 69], [68; 67], [74; 73], ...
    [18; 8], [100; 60], [5; 2], ...
    [10; 7], [12; 8], [6; 5], [28; 20], ...
    'VariableNames', {'Method', 'dose_targetD95Min', 'dose_targetD98Min', 'dose_oarDmaxMax', ...
    'delivery_nActiveLayers', 'delivery_nActiveSpots', 'delivery_nSwitchUp', ...
    'delivery_T_energy', 'delivery_T_spot', 'delivery_T_gantry', 'delivery_T_total'});
latexText = export_latex_table(mockTable);
checks.latexHasDisclaimer = contains(latexText, researchMeta.disclaimer, 'IgnoreCase', false);
checks.latexHasMachineProxyNote = contains(latexText, researchMeta.machineModelNote, 'IgnoreCase', false);
checks.latexHasCitationPlaceholder = contains(latexText, researchMeta.citationPlaceholder, 'IgnoreCase', false);

report = struct();
report.researchMeta = researchMeta;
report.checks = checks;
report.passed = all(cell2mat(struct2cell(checks)));

timestamp = datestr(now, 'yyyy-mm-dd_HHMMSS');
matFile = fullfile(outputDir, sprintf('task42_documentation_check_%s.mat', timestamp));
txtFile = fullfile(outputDir, sprintf('task42_documentation_check_%s.txt', timestamp));
save(matFile, 'report', '-v7.3');
write_report(txtFile, report);

outputs = struct();
outputs.report = report;
outputs.files = struct('matFile', matFile, 'txtFile', txtFile);
end

function write_report(txtFile, report)
fid = fopen(txtFile, 'w');
if fid == -1
    error('run_task42_documentation_check:fileOpenFailed', ...
        'Could not open %s for writing.', txtFile);
end
cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>

fprintf(fid, '%s\n', report.researchMeta.disclaimer);
fprintf(fid, 'Machine model note: %s\n', report.researchMeta.machineModelNote);
fprintf(fid, 'Citation placeholder: %s\n\n', report.researchMeta.citationPlaceholder);
fprintf(fid, 'Task 42 checks:\n');
checkNames = fieldnames(report.checks);
for i = 1:numel(checkNames)
    fprintf(fid, '  %s: %d\n', checkNames{i}, report.checks.(checkNames{i}));
end
fprintf(fid, '\nPassed: %d\n', report.passed);
end
