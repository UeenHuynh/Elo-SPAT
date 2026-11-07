function outputs = run_task39_robust_integration(cfgOverride, outputDir)
% run_task39_robust_integration  Run and verify Task 39 robust integration.
%
% This helper executes the orchestrator with robust evaluation enabled and
% verifies that worst-case CTV D95 is reported and remains within the
% configured robust tolerance. It also reports scenario pass rate.

narginchk(0, 2);

planningDir = fileparts(mfilename('fullpath'));
workspaceDir = fileparts(planningDir);
resultsDir = fullfile(workspaceDir, 'results');

addpath(workspaceDir);

cfg = struct();
cfg.arcSpacing = 20;
cfg.captureConsoleOutput = false;
cfg.runtime = struct( ...
    'headless', true, ...
    'forceSoftwareOpenGL', false);
cfg.robust = struct( ...
    'enabled', true, ...
    'mode', 'COWC', ...
    'setupMM', [3 3 3], ...
    'rangePercent', 3.5);

if nargin >= 1 && ~isempty(cfgOverride)
    cfg = merge_override(cfg, cfgOverride);
end

if nargin >= 2 && ~isempty(outputDir)
    cfg.outputDir = outputDir;
elseif ~isfield(cfg, 'outputDir') || isempty(cfg.outputDir)
    cfg.outputDir = fullfile(resultsDir, sprintf('task39_%s', datestr(now, 'yyyy-mm-dd_HHMMSS')));
end

fprintf('\n========================================\n');
fprintf('Task 39 Robust Integration Runner\n');
fprintf('========================================\n');
fprintf('Output directory: %s\n\n', cfg.outputDir);

startTime = tic;
mainOutputs = run_rms_elospat_main(cfg);
runtimeSeconds = toc(startTime);

if ~isfield(mainOutputs, 'robust') || ~isstruct(mainOutputs.robust) || ...
        ~isfield(mainOutputs.robust, 'status') || ~strcmp(mainOutputs.robust.status, 'completed')
    error('run_task39_robust_integration:robustFailed', ...
        'Robust stage did not complete successfully.');
end

robustData = load(mainOutputs.robust.matFile, 'robustEval', 'nominalMetrics');
nominalD95 = robustData.nominalMetrics.summary.targetD95Min;
worstCaseD95 = robustData.robustEval.worstCase.targetD95Min;
tolerancePercent = double(mainOutputs.cfg.evaluation.robustD95TolerancePercent);
d95DegradationPercent = 100 * (double(nominalD95) - double(worstCaseD95)) / max(double(nominalD95), eps);

verification = struct();
verification.robustCompleted = true;
verification.worstCaseComputed = isfinite(worstCaseD95);
verification.ctvD95WithinTolerance = d95DegradationPercent <= tolerancePercent;
verification.passRate = double(robustData.robustEval.passRate);
verification.nScenarios = double(robustData.robustEval.nScenarios);
verification.runtimeSeconds = runtimeSeconds;
verification.d95DegradationPercent = d95DegradationPercent;
verification.tolerancePercent = tolerancePercent;
if isfield(mainOutputs.robust, 'warmStartSource')
    verification.warmStartSource = mainOutputs.robust.warmStartSource;
else
    verification.warmStartSource = 'unknown';
end

report = struct();
report.researchMeta = research_only_metadata();
report.taskId = 39;
report.verification = verification;
report.robustMatFile = mainOutputs.robust.matFile;
report.robustDvhPlot = mainOutputs.robust.dvhPlot;
report.robustSummaryPlot = mainOutputs.robust.summaryPlot;
report.runtimeSeconds = runtimeSeconds;

reportFile = fullfile(cfg.outputDir, 'task39_robust_integration_report.mat');
save(reportFile, 'report', '-v7.3');

outputs = struct();
outputs.mainOutputs = mainOutputs;
outputs.report = report;
outputs.reportFile = reportFile;

fprintf('Task 39 report saved: %s\n', reportFile);
fprintf('Worst-case CTV D95: %.6f\n', worstCaseD95);
fprintf('Robust pass rate: %.3f\n', verification.passRate);
fprintf('Robust evaluation source: %s\n', verification.warmStartSource);
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
