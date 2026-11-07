function outputs = run_task41_arc_spacing_validation(spacingList, cfgOverride, outputDir)
% run_task41_arc_spacing_validation  Run Task 41 spacing comparison.
%
% This helper runs the nominal orchestrator across several arc spacings,
% records estimated BDT reduction, and documents memory-related failures.

narginchk(0, 3);

planningDir = fileparts(mfilename('fullpath'));
workspaceDir = fileparts(planningDir);
resultsDir = fullfile(workspaceDir, 'results');

addpath(workspaceDir);

if nargin < 1 || isempty(spacingList)
    spacingList = [20, 10, 5];
end
spacingList = reshape(double(spacingList), 1, []);

cfgBase = struct();
cfgBase.captureConsoleOutput = false;
cfgBase.runtime = struct('headless', true, 'forceSoftwareOpenGL', false);
cfgBase.robust = struct('enabled', false);

if nargin >= 2 && ~isempty(cfgOverride)
    cfgBase = merge_override(cfgBase, cfgOverride);
end
if ~isfield(cfgBase, 'robust') || ~isstruct(cfgBase.robust)
    cfgBase.robust = struct();
end
cfgBase.robust.enabled = false;

if nargin >= 3 && ~isempty(outputDir)
    outputRoot = outputDir;
else
    outputRoot = fullfile(resultsDir, sprintf('task41_%s', datestr(now, 'yyyy-mm-dd_HHMMSS')));
end
outputRoot = resolve_output_root(outputRoot, workspaceDir);
if exist(outputRoot, 'dir') ~= 7
    mkdir(outputRoot);
end

fprintf('\n========================================\n');
fprintf('Task 41 Arc-Spacing Validation Runner\n');
fprintf('========================================\n');
fprintf('Output root: %s\n\n', outputRoot);

records = repmat(struct( ...
    'arcSpacing', NaN, ...
    'status', "", ...
    'errorMessage', "", ...
    'memoryInsufficient', false, ...
    'runtimeSeconds', NaN, ...
    'denseBDT', NaN, ...
    'rmsBDT', NaN, ...
    'bdtReductionPercent', NaN, ...
    'estimatedDijMemoryMB', NaN, ...
    'totalNumOfBixels', NaN, ...
    'resultMatFile', ""), numel(spacingList), 1);

for i = 1:numel(spacingList)
    arcSpacing = spacingList(i);
    cfgRun = cfgBase;
    cfgRun.arcSpacing = arcSpacing;
    cfgRun.outputDir = fullfile(outputRoot, sprintf('arc_%gdeg', arcSpacing));
    if exist(cfgRun.outputDir, 'dir') ~= 7
        mkdir(cfgRun.outputDir);
    end

    records(i).arcSpacing = arcSpacing;
    fprintf('Running arc spacing %.2f deg...\n', arcSpacing);
    startTime = tic;
    try
        runOutputs = run_rms_elospat_main(cfgRun);
        runtimeSeconds = toc(startTime);
        [denseRow, rmsRow] = resolve_dense_and_rms_rows(runOutputs.methodTable);
        records(i).status = "completed";
        records(i).runtimeSeconds = runtimeSeconds;
        records(i).denseBDT = double(denseRow.delivery_T_total);
        records(i).rmsBDT = double(rmsRow.delivery_T_total);
        records(i).bdtReductionPercent = 100 * (records(i).denseBDT - records(i).rmsBDT) / max(records(i).denseBDT, eps);
        records(i).estimatedDijMemoryMB = estimate_dij_memory_mb_local(runOutputs.stages.arc.dij);
        records(i).totalNumOfBixels = double(runOutputs.stages.arc.dij.totalNumOfBixels);
        records(i).resultMatFile = string(runOutputs.files.matFile);
    catch ME
        records(i).status = "failed";
        records(i).errorMessage = string(ME.message);
        records(i).memoryInsufficient = contains(ME.message, 'memory', 'IgnoreCase', true) || ...
            contains(ME.message, 'Out of memory', 'IgnoreCase', true);
        fprintf(2, '  FAILED: %s\n', ME.message);
    end
end

summaryTable = struct2table(records);
csvFile = fullfile(outputRoot, 'task41_arc_spacing_summary.csv');
matFile = fullfile(outputRoot, 'task41_arc_spacing_summary.mat');
plotFile = fullfile(outputRoot, 'task41_arc_spacing_bdt_curve.png');
save_spacing_plot(summaryTable, plotFile);

writetable(summaryTable, csvFile);
researchMeta = research_only_metadata();
save(matFile, 'summaryTable', 'records', 'researchMeta', '-v7.3');

outputs = struct();
outputs.summaryTable = summaryTable;
outputs.files = struct('csvFile', csvFile, 'matFile', matFile, 'plotFile', plotFile);

fprintf('\nTask 41 summary saved: %s\n', matFile);
fprintf('Task 41 CSV saved: %s\n', csvFile);
fprintf('Task 41 plot saved: %s\n\n', plotFile);
end

function merged = merge_override(base, override)
merged = base;
fields = fieldnames(override);
for i = 1:numel(fields)
    key = fields{i};
    if isfield(base, key) && isstruct(base.(key)) && isstruct(override.(key))
        merged.(key) = merge_override(base.(key), override.(key));
    else
        merged.(key) = override.(key);
    end
end
end

function outputRoot = resolve_output_root(outputRoot, workspaceDir)
if is_absolute_path(outputRoot)
    return;
end
outputRoot = fullfile(workspaceDir, outputRoot);
end

function tf = is_absolute_path(pathValue)
tf = false;
if ~ischar(pathValue) && ~(isstring(pathValue) && isscalar(pathValue))
    return;
end
pathValue = char(pathValue);
if isempty(pathValue)
    return;
end
if pathValue(1) == '/' || pathValue(1) == '\'
    tf = true;
    return;
end
if numel(pathValue) >= 3 && pathValue(2) == ':' && (pathValue(3) == '\' || pathValue(3) == '/')
    tf = true;
end
end

function [denseRow, rmsRow] = resolve_dense_and_rms_rows(methodTable)
denseMask = strcmp(string(methodTable.Method), "Dense Arc");
rmsMask = strcmp(string(methodTable.Method), "RMS-ELO-SPAT");
if ~any(denseMask) || ~any(rmsMask)
    error('run_task41_arc_spacing_validation:missingMethods', ...
        'Method table must contain Dense Arc and RMS-ELO-SPAT.');
end
denseRow = methodTable(find(denseMask, 1, 'first'), :);
rmsRow = methodTable(find(rmsMask, 1, 'first'), :);
end

function estimatedMemoryMB = estimate_dij_memory_mb_local(dij)
estimatedMemoryMB = NaN;
if ~isstruct(dij) || ~isfield(dij, 'doseGrid') || ~isfield(dij.doseGrid, 'numOfVoxels') || ...
        ~isfield(dij, 'totalNumOfBixels') || isempty(dij.doseGrid.numOfVoxels) || isempty(dij.totalNumOfBixels)
    return;
end

estimatedNonzeros = double(dij.doseGrid.numOfVoxels) * double(dij.totalNumOfBixels) * 0.01;
if isfield(dij, 'physicalDose') && iscell(dij.physicalDose)
    nScenarios = numel(dij.physicalDose);
else
    nScenarios = 1;
end
estimatedMemoryMB = (estimatedNonzeros * 8 * nScenarios) / (1024^2);
end

function save_spacing_plot(summaryTable, plotFile)
fig = figure('Visible', 'off', 'Color', 'w');
successMask = strcmp(summaryTable.status, "completed");

if any(successMask)
    plot(summaryTable.arcSpacing(successMask), summaryTable.bdtReductionPercent(successMask), ...
        'o-', 'LineWidth', 1.6, 'MarkerSize', 7);
    hold on;
    failedMask = ~successMask;
    if any(failedMask)
        scatter(summaryTable.arcSpacing(failedMask), zeros(nnz(failedMask), 1), 60, 'x', 'LineWidth', 1.5);
        legend({'Completed', 'Failed / memory-limited'}, 'Location', 'best');
    end
    ylabel('Estimated BDT reduction (%)');
else
    text(0.1, 0.5, 'No successful arc-spacing runs were completed.', 'FontSize', 11);
    axis off;
end

xlabel('Arc spacing (deg)');
title('Task 41 arc spacing vs estimated BDT reduction');
grid on;
set(fig, 'PaperPositionMode', 'auto');
print(fig, plotFile, '-dpng', '-r150');
close(fig);
end
