function outputs = run_task38_integration(cfgOverride, outputDir)
% run_task38_integration  Run and verify Task 38 nominal integration.
%
% This helper executes the full nominal orchestrator on PROSTATE.mat with a
% 20-degree arc, then verifies:
%   - all 5 methods are present
%   - RMS-ELO-SPAT estimated BDT is below Dense Arc estimated BDT
%   - RMS-ELO-SPAT CTV D95 stays within the configured tolerance
%   - comparison CSV and LaTeX outputs exist

narginchk(0, 2);

planningDir = fileparts(mfilename('fullpath'));
workspaceDir = fileparts(planningDir);
resultsDir = fullfile(workspaceDir, 'results');

addpath(workspaceDir);

cfg = struct();
cfg.arcSpacing = 20;
cfg.captureConsoleOutput = false;
cfg.runtime = struct('headless', true, 'forceSoftwareOpenGL', false);
cfg.robust = struct('enabled', false);

if nargin >= 1 && ~isempty(cfgOverride)
    cfg = merge_override(cfg, cfgOverride);
end

if nargin >= 2 && ~isempty(outputDir)
    cfg.outputDir = outputDir;
elseif ~isfield(cfg, 'outputDir') || isempty(cfg.outputDir)
    cfg.outputDir = fullfile(resultsDir, sprintf('task38_%s', datestr(now, 'yyyy-mm-dd_HHMMSS')));
end

fprintf('\n========================================\n');
fprintf('Task 38 Integration Runner\n');
fprintf('========================================\n');
fprintf('Output directory: %s\n\n', cfg.outputDir);

startTime = tic;
mainOutputs = run_rms_elospat_main(cfg);
runtimeSeconds = toc(startTime);

methodTable = mainOutputs.methodTable;
requiredMethods = ["IMPT", "Dense Arc", "ELO-SPAT-like", "Post-filter", "RMS-ELO-SPAT"];
presentMethods = string(methodTable.Method);
allMethodsPresent = all(ismember(requiredMethods, presentMethods));

[denseRow, rmsRow] = resolve_dense_and_rms_rows(methodTable);
tolerancePercent = double(mainOutputs.cfg.evaluation.ctvD95TolerancePercent);
d95DegradationPercent = 100 * (double(denseRow.dose_targetD95Min) - double(rmsRow.dose_targetD95Min)) / ...
    max(double(denseRow.dose_targetD95Min), eps);
bdtReductionPercent = 100 * (double(denseRow.delivery_T_total) - double(rmsRow.delivery_T_total)) / ...
    max(double(denseRow.delivery_T_total), eps);

verification = struct();
verification.allMethodsPresent = allMethodsPresent;
verification.rmsBdtBelowDense = double(rmsRow.delivery_T_total) < double(denseRow.delivery_T_total);
verification.ctvD95WithinTolerance = d95DegradationPercent <= tolerancePercent;
verification.comparisonCsvExists = exist(mainOutputs.files.methodCsv, 'file') == 2;
verification.latexExists = exist(mainOutputs.files.latexFile, 'file') == 2;
verification.summaryPlotExists = exist(mainOutputs.files.summaryPlot, 'file') == 2;
verification.runtimeSeconds = runtimeSeconds;
verification.d95DegradationPercent = d95DegradationPercent;
verification.bdtReductionPercent = bdtReductionPercent;
verification.tolerancePercent = tolerancePercent;

report = struct();
report.researchMeta = research_only_metadata();
report.taskId = 38;
report.methodTable = methodTable;
report.verification = verification;
report.runtimeSeconds = runtimeSeconds;
report.outputFiles = mainOutputs.files;

reportFile = fullfile(cfg.outputDir, 'task38_integration_report.mat');
save(reportFile, 'report', '-v7.3');

outputs = struct();
outputs.mainOutputs = mainOutputs;
outputs.report = report;
outputs.reportFile = reportFile;

fprintf('Task 38 report saved: %s\n', reportFile);
fprintf('All methods present: %d\n', verification.allMethodsPresent);
fprintf('RMS BDT below Dense Arc: %d\n', verification.rmsBdtBelowDense);
fprintf('CTV D95 degradation: %.3f%% (tol %.3f%%)\n', ...
    verification.d95DegradationPercent, verification.tolerancePercent);
fprintf('Runtime: %.2f minutes\n\n', runtimeSeconds / 60);
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

function [denseRow, rmsRow] = resolve_dense_and_rms_rows(methodTable)
denseMask = strcmp(string(methodTable.Method), "Dense Arc");
rmsMask = strcmp(string(methodTable.Method), "RMS-ELO-SPAT");
if ~any(denseMask) || ~any(rmsMask)
    error('run_task38_integration:missingMethods', ...
        'Method table must contain Dense Arc and RMS-ELO-SPAT.');
end
denseRow = methodTable(find(denseMask, 1, 'first'), :);
rmsRow = methodTable(find(rmsMask, 1, 'first'), :);
end
