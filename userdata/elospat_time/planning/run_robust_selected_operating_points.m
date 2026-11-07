function summaryTable = run_robust_selected_operating_points(outputDir)
% Evaluate the three frozen operating points under the configured COWC model.
% Research-only; this is sequential evaluation, not robust re-optimization.

narginchk(0, 1);
rootDir = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(rootDir, 'config'));
addpath(fullfile(rootDir, 'planning'));
addpath(fullfile(rootDir, 'mapping'));
addpath(fullfile(rootDir, 'delivery'));
addpath(fullfile(rootDir, 'reduction'));
addpath(fullfile(rootDir, 'objective'));
addpath(fullfile(rootDir, 'evaluation'));
addpath(fullfile(rootDir, 'robustness'));

if nargin < 1 || isempty(outputDir)
    outputDir = fullfile(rootDir, 'results', 'robust_selected_operating_points');
end
if exist(outputDir, 'dir') ~= 7
    mkdir(outputDir);
end

specs = repmat(struct('spacing', 0, 'sourceArtifact', '', ...
    'selectedFile', '', 'supportSpots', 0), 3, 1);
specs(1).spacing = 20;
specs(1).sourceArtifact = fullfile(rootDir, 'results', 'support_fix_20deg', ...
    'arc_20deg', 'rms_elospat_main_2026-09-15_105936_nominal_artifact.mat');
specs(1).selectedFile = fullfile(rootDir, 'results', ...
    'support_quality_sweep_20deg', 'final_05544.mat');
specs(1).supportSpots = 5544;
specs(2).spacing = 10;
specs(2).sourceArtifact = fullfile(rootDir, 'results', ...
    'full_nominal_then_robust_2026-06-02_003113', 'nominal', 'task41', ...
    'arc_10deg', 'rms_elospat_main_2026-06-03_172706_nominal_artifact.mat');
specs(2).selectedFile = fullfile(rootDir, 'results', ...
    'support_quality_sweep_10deg', 'final_08888.mat');
specs(2).supportSpots = 8888;
specs(3).spacing = 5;
specs(3).sourceArtifact = fullfile(rootDir, 'results', ...
    'full_nominal_then_robust_2026-06-02_003113', 'nominal', 'task41', ...
    'arc_5deg', 'rms_elospat_main_2026-06-04_022851_nominal_artifact.mat');
specs(3).selectedFile = fullfile(rootDir, 'results', ...
    'support_quality_sweep_5deg', 'final_09097.mat');
specs(3).supportSpots = 9097;

records = repmat(struct( ...
    'arcSpacingDeg', NaN, 'supportSpots', NaN, 'nScenarios', NaN, ...
    'setupMM', "", 'rangePercent', NaN, 'mode', "COWC", ...
    'nominalD95', NaN, 'nominalD98', NaN, 'nominalOarDmax', NaN, ...
    'worstD95', NaN, 'worstD98', NaN, 'worstOarDmax', NaN, ...
    'worstD95Scenario', NaN, 'worstD98Scenario', NaN, ...
    'worstOarScenario', NaN, 'passRate', NaN, 'runtimeSeconds', NaN, ...
    'status', "pending", 'errorMessage', ""), numel(specs), 1);

for i = 1:numel(specs)
    spec = specs(i);
    fprintf('\n=== Robust operating point %g deg (%d spots) ===\n', ...
        spec.spacing, spec.supportSpots);
    records(i).arcSpacingDeg = spec.spacing;
    records(i).supportSpots = spec.supportSpots;
    try
        derivedDir = fullfile(outputDir, sprintf('arc_%gdeg', spec.spacing));
        if exist(derivedDir, 'dir') ~= 7
            mkdir(derivedDir);
        end
        derivedArtifact = fullfile(derivedDir, 'selected_nominal_artifact.mat');
        make_selected_artifact(spec.sourceArtifact, spec.selectedFile, derivedArtifact);

        cfgData = load(derivedArtifact, 'cfg');
        cfg = cfgData.cfg;
        cfg.robust.enabled = true;
        cfg.robust.mode = 'COWC';
        cfg.robust.setupMM = [3, 3, 3];
        cfg.robust.rangePercent = 3.5;
        save(derivedArtifact, 'cfg', '-append');

        startTime = tic;
        robust = run_task28_30_from_nominal(derivedArtifact, derivedDir, false);
        records(i).runtimeSeconds = toc(startTime);
        records(i).nScenarios = robust.nScenarios;
        records(i).passRate = robust.passRate;
        records(i).worstD95 = robust.worstCase.targetD95Min;
        records(i).worstD98 = robust.worstCase.targetD98Min;
        records(i).worstOarDmax = robust.worstCase.oarDmaxMax;
        records(i).worstD95Scenario = robust.worstCase.worstTargetD95Scenario;
        records(i).worstD98Scenario = robust.worstCase.worstTargetD98Scenario;
        records(i).worstOarScenario = robust.worstCase.worstOarDmaxScenario;
        result = load(robust.matFile, 'nominalMetrics', 'cfg');
        records(i).nominalD95 = result.nominalMetrics.summary.targetD95Min;
        records(i).nominalD98 = result.nominalMetrics.summary.targetD98Min;
        records(i).nominalOarDmax = result.nominalMetrics.summary.oarDmaxMax;
        records(i).setupMM = strjoin(string(cfg.robust.setupMM), ',');
        records(i).rangePercent = cfg.robust.rangePercent;
        records(i).status = "completed";
    catch ME
        records(i).status = "failed";
        records(i).errorMessage = string(ME.message);
        fprintf(2, 'Robust %g deg failed: %s\n', spec.spacing, ME.message);
    end
    writetable(struct2table(records), fullfile(outputDir, 'robust_selected_summary.csv'));
end

summaryTable = struct2table(records);
if any(strcmp(summaryTable.status, "failed"))
    error('run_robust_selected_operating_points:incomplete', ...
        'At least one robust operating point failed; inspect robust_selected_summary.csv.');
end
end

function make_selected_artifact(sourceArtifact, selectedFile, destination)
source = load(sourceArtifact);
selected = load(selectedFile, 'wFinal', 'supportMask');
if ~isfield(source, 'resultFinal') || ~isfield(source.resultFinal, 'w')
    error('run_robust_selected_operating_points:missingFinal', ...
        'Source artifact lacks resultFinal.w: %s', sourceArtifact);
end
if numel(source.resultFinal.w) ~= numel(selected.wFinal)
    error('run_robust_selected_operating_points:lengthMismatch', ...
        'Selected weight length does not match source artifact.');
end
source.resultFinal.w = reshape(double(selected.wFinal), [], 1);
if isfield(source.resultFinal, 'resultGUI')
    source.resultFinal.resultGUI.w = source.resultFinal.w;
end
source.resultFinal.supportMask = logical(selected.supportMask);
source.selectedSupportMask = logical(selected.supportMask);
source.selectedSupportSpots = nnz(source.selectedSupportMask);
save(destination, '-struct', 'source', '-v7.3');
end
