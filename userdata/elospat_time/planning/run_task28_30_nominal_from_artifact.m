function outputs = run_task28_30_nominal_from_artifact(artifactFile, patientFile, outputDir)
% run_task28_30_nominal_from_artifact  Compute nominal-only Task 28-30 style outputs from a saved 25-27 artifact.
%
% This is a nominal proxy for Tasks 28-30. It does not evaluate uncertainty
% scenarios and should not be interpreted as a robust result.

narginchk(0, 3);

planningDir = fileparts(mfilename('fullpath'));
workspaceDir = fileparts(planningDir);
configDir = fullfile(workspaceDir, 'config');
evaluationDir = fullfile(workspaceDir, 'evaluation');
resultsDir = fullfile(workspaceDir, 'results');

addpath(configDir);
addpath(evaluationDir);

if nargin < 1 || isempty(artifactFile)
    artifactFile = resolve_latest_artifact(resultsDir);
end

if nargin < 2 || isempty(patientFile)
    patientFile = resolve_patient_file();
end

if nargin < 3 || isempty(outputDir)
    outputDir = resultsDir;
end

if exist(outputDir, 'dir') ~= 7
    mkdir(outputDir);
end

fprintf('\n========================================\n');
fprintf('Task 28-30 Nominal Proxy\n');
fprintf('========================================\n');
fprintf('Artifact: %s\n', artifactFile);
fprintf('Patient:  %s\n', patientFile);
fprintf('Output:   %s\n\n', outputDir);

data = load(artifactFile, 'cfg', 'denseResult', 'resultFinal');
patientData = load(patientFile, 'ct', 'cst');
cfg = resolve_cfg(data);

denseDose = resolve_dose_cube(data.denseResult, 'denseResult');
finalDose = resolve_dose_cube(data.resultFinal.resultGUI, 'resultFinal.resultGUI');

denseMetrics = compute_dvh_metrics(denseDose, patientData.ct, patientData.cst, cfg);
finalMetrics = compute_dvh_metrics(finalDose, patientData.ct, patientData.cst, cfg);

[~, artifactBaseName] = fileparts(artifactFile);
plotPrefix = sprintf('%s_task28_30_nominal', artifactBaseName);

plotFiles = struct();
plotFiles.dvhOverlay = save_nominal_overlay_plot(patientData.cst, denseDose, finalDose, outputDir, plotPrefix);
plotFiles.metricBars = save_nominal_metric_plot(denseMetrics, finalMetrics, outputDir, plotPrefix);

summary = struct();
summary.mode = 'nominal_only_proxy';
summary.artifactFile = artifactFile;
summary.patientFile = patientFile;
summary.timestamp = datetime('now');
summary.dense = denseMetrics.summary;
summary.final = finalMetrics.summary;

matFile = fullfile(outputDir, sprintf('%s.mat', plotPrefix));
save(matFile, 'cfg', 'denseMetrics', 'finalMetrics', 'summary', 'plotFiles', '-v7.3');

outputs = struct();
outputs.matFile = matFile;
outputs.plotFiles = plotFiles;
outputs.summary = summary;

fprintf('Saved MAT: %s\n', matFile);
fprintf('Saved DVH overlay: %s\n', plotFiles.dvhOverlay);
fprintf('Saved metric bars: %s\n\n', plotFiles.metricBars);
end

function artifactFile = resolve_latest_artifact(resultsDir)
listing = dir(fullfile(resultsDir, 'task25_27_pipeline_*.mat'));
if isempty(listing)
    error('run_task28_30_nominal_from_artifact:missingArtifact', ...
        'No task25_27_pipeline_*.mat artifact found in %s.', resultsDir);
end

[~, ix] = max([listing.datenum]);
artifactFile = fullfile(listing(ix).folder, listing(ix).name);
end

function patientFile = resolve_patient_file()
patientFile = which('PROSTATE.mat');
if ~isempty(patientFile)
    return;
end

matRadRcPath = which('matRad_rc');
if isempty(matRadRcPath)
    error('run_task28_30_nominal_from_artifact:missingPatient', ...
        'PROSTATE.mat was not found and matRad_rc is not on the MATLAB path.');
end

matRadRoot = fileparts(matRadRcPath);
candidate = fullfile(matRadRoot, 'data', 'PROSTATE.mat');
if exist(candidate, 'file') ~= 2
    error('run_task28_30_nominal_from_artifact:missingPatient', ...
        'Could not locate PROSTATE.mat at %s.', candidate);
end

patientFile = candidate;
end

function cfg = resolve_cfg(data)
if isfield(data, 'cfg') && isstruct(data.cfg)
    cfg = data.cfg;
else
    cfg = default_config();
end
end

function doseCube = resolve_dose_cube(resultStruct, label)
if isfield(resultStruct, 'RBExDose') && ~isempty(resultStruct.RBExDose)
    doseCube = resultStruct.RBExDose;
elseif isfield(resultStruct, 'physicalDose') && ~isempty(resultStruct.physicalDose)
    doseCube = resultStruct.physicalDose;
else
    error('run_task28_30_nominal_from_artifact:missingDose', ...
        '%s does not contain RBExDose or physicalDose.', label);
end
end

function plotPath = save_nominal_overlay_plot(cst, denseDose, finalDose, outputDir, plotPrefix)
dvhDense = matRad_calcDVH(cst, denseDose, 'cum');
dvhFinal = matRad_calcDVH(cst, finalDose, 'cum');

fig = figure('Visible', 'off', 'Color', 'w');
hold on;

for i = 1:numel(dvhDense)
    plot(dvhDense(i).doseGrid, dvhDense(i).volumePoints, '--', 'LineWidth', 1.2);
    plot(dvhFinal(i).doseGrid, dvhFinal(i).volumePoints, '-', 'LineWidth', 1.6);
end

legendLabels = cell(2 * size(cst, 1), 1);
for i = 1:size(cst, 1)
    baseName = stringify_label(cst{i, 2}, sprintf('Structure %d', i));
    legendLabels{2 * i - 1} = sprintf('%s dense', baseName);
    legendLabels{2 * i} = sprintf('%s final', baseName);
end

legend(legendLabels, 'Interpreter', 'none', 'Location', 'bestoutside');
xlabel('Dose');
ylabel('Volume (%)');
title('Nominal DVH overlay for Task 28-30 proxy');
grid on;
box on;

plotPath = fullfile(outputDir, sprintf('%s_dvh_overlay.png', plotPrefix));
set(fig, 'PaperPositionMode', 'auto');
print(fig, plotPath, '-dpng', '-r150');
close(fig);
end

function plotPath = save_nominal_metric_plot(denseMetrics, finalMetrics, outputDir, plotPrefix)
metricNames = {'targetD95Min', 'targetD98Min', 'targetV95Min', 'oarDmaxMax'};
metricLabels = {'CTV D95', 'CTV D98', 'CTV V95', 'OAR Dmax'};
denseValues = cellfun(@(name) denseMetrics.summary.(name), metricNames);
finalValues = cellfun(@(name) finalMetrics.summary.(name), metricNames);

fig = figure('Visible', 'off', 'Color', 'w');
bar(categorical(metricLabels), [denseValues(:), finalValues(:)], 'grouped');
ylabel('Metric value');
title('Nominal summary metrics: Dense vs Final');
legend({'Dense', 'Final'}, 'Location', 'best');
grid on;

plotPath = fullfile(outputDir, sprintf('%s_metric_bars.png', plotPrefix));
set(fig, 'PaperPositionMode', 'auto');
print(fig, plotPath, '-dpng', '-r150');
close(fig);
end

function value = stringify_label(rawValue, fallbackValue)
if isstring(rawValue)
    value = char(rawValue);
elseif ischar(rawValue)
    value = rawValue;
else
    value = fallbackValue;
end

value = strtrim(value);
if isempty(value)
    value = fallbackValue;
end
end
