function plotFiles = plot_nominal_pipeline_from_artifact(artifactFile, patientFile, outputDir)
% plot_nominal_pipeline_from_artifact  Rebuild nominal pipeline plots from a saved Task 25-27 artifact.

narginchk(0, 3);

planningDir = fileparts(mfilename('fullpath'));
workspaceDir = fileparts(planningDir);
configDir = fullfile(workspaceDir, 'config');
deliveryDir = fullfile(workspaceDir, 'delivery');
evaluationDir = fullfile(workspaceDir, 'evaluation');
reductionDir = fullfile(workspaceDir, 'reduction');

addpath(configDir);
addpath(deliveryDir);
addpath(evaluationDir);
addpath(reductionDir);

resultsDir = fullfile(workspaceDir, 'results');

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
fprintf('Nominal Plot Rebuild\n');
fprintf('========================================\n');
fprintf('Artifact: %s\n', artifactFile);
fprintf('Patient:  %s\n', patientFile);
fprintf('Output:   %s\n\n', outputDir);

data = load(artifactFile, ...
    'cfg', 'denseResult', 'resultLayer', 'resultSpot', 'resultFinal', ...
    'map', 'groups', 'machine');
patientData = load(patientFile, 'ct', 'cst');

cfg = resolve_cfg(data);
machine = resolve_machine(data);

denseDose = resolve_dose_cube(data.denseResult, 'denseResult');
finalDose = resolve_dose_cube(data.resultFinal.resultGUI, 'resultFinal.resultGUI');

[~, artifactBaseName] = fileparts(artifactFile);
plotPrefix = sprintf('%s_nominal', artifactBaseName);

plotFiles = struct();
plotFiles.dvhOverlay = save_nominal_dvh_overlay( ...
    patientData.cst, denseDose, finalDose, outputDir, plotPrefix);
plotFiles.stageSummary = save_stage_summary_plot( ...
    data, cfg, machine, outputDir, plotPrefix);

fprintf('Saved DVH overlay: %s\n', plotFiles.dvhOverlay);
fprintf('Saved stage summary: %s\n\n', plotFiles.stageSummary);
end

function artifactFile = resolve_latest_artifact(resultsDir)
listing = dir(fullfile(resultsDir, 'task25_27_pipeline_*.mat'));
if isempty(listing)
    error('plot_nominal_pipeline_from_artifact:missingArtifact', ...
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
    error('plot_nominal_pipeline_from_artifact:missingPatient', ...
        'PROSTATE.mat was not found and matRad_rc is not on the MATLAB path.');
end

matRadRoot = fileparts(matRadRcPath);
candidate = fullfile(matRadRoot, 'data', 'PROSTATE.mat');
if exist(candidate, 'file') ~= 2
    error('plot_nominal_pipeline_from_artifact:missingPatient', ...
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

function machine = resolve_machine(data)
if isfield(data, 'machine') && isstruct(data.machine)
    machine = data.machine;
else
    machine = default_machine_model();
end
end

function doseCube = resolve_dose_cube(resultStruct, label)
if isfield(resultStruct, 'RBExDose') && ~isempty(resultStruct.RBExDose)
    doseCube = resultStruct.RBExDose;
elseif isfield(resultStruct, 'physicalDose') && ~isempty(resultStruct.physicalDose)
    doseCube = resultStruct.physicalDose;
else
    error('plot_nominal_pipeline_from_artifact:missingDose', ...
        '%s does not contain RBExDose or physicalDose.', label);
end
end

function plotPath = save_nominal_dvh_overlay(cst, denseDose, finalDose, outputDir, plotPrefix)
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
title('Nominal DVH overlay: Dense Arc vs Final');
grid on;
box on;

plotPath = fullfile(outputDir, sprintf('%s_dvh_overlay.png', plotPrefix));
set(fig, 'PaperPositionMode', 'auto');
print(fig, plotPath, '-dpng', '-r150');
close(fig);
end

function plotPath = save_stage_summary_plot(data, cfg, machine, outputDir, plotPrefix)
stageNames = {'Dense', 'After Task25', 'After Task26', 'Final'};
bdtTotals = nan(1, numel(stageNames));
activeLayers = nan(1, numel(stageNames));
activeSpots = nan(1, numel(stageNames));

if isfield(data, 'denseResult') && isfield(data.denseResult, 'w')
    [bdtTotals(1), activeLayers(1), activeSpots(1)] = compute_stage_delivery_metrics( ...
        data.denseResult.w, data.map, data.groups, cfg, machine);
end

if isfield(data, 'resultLayer')
    if isfield(data.resultLayer, 'bdtMetrics') && isstruct(data.resultLayer.bdtMetrics) && ...
            isfield(data.resultLayer.bdtMetrics, 'T_total')
        bdtTotals(2) = double(data.resultLayer.bdtMetrics.T_total);
    end
    if isfield(data.resultLayer, 'layerMask')
        activeLayers(2) = sum(logical(data.resultLayer.layerMask));
    end
    if isfield(data.resultLayer, 'bixelMask')
        activeSpots(2) = sum(logical(data.resultLayer.bixelMask));
    end
end

if isfield(data, 'resultSpot')
    if isfield(data.resultSpot, 'bdtMetrics') && isstruct(data.resultSpot.bdtMetrics) && ...
            isfield(data.resultSpot.bdtMetrics, 'T_total')
        bdtTotals(3) = double(data.resultSpot.bdtMetrics.T_total);
    end
    if isfield(data, 'groups') && isfield(data.resultSpot, 'w')
        activeLayers(3) = sum(build_layer_mask(data.resultSpot.w, data.groups, cfg.layerThreshold));
    end
    if isfield(data.resultSpot, 'spotMask')
        activeSpots(3) = sum(logical(data.resultSpot.spotMask));
    end
end

if isfield(data, 'resultFinal') && isfield(data.resultFinal, 'w')
    [bdtTotals(4), activeLayers(4), activeSpots(4)] = compute_stage_delivery_metrics( ...
        data.resultFinal.w, data.map, data.groups, cfg, machine);
end

fig = figure('Visible', 'off', 'Color', 'w');
tiledlayout(2, 1);

nexttile;
yyaxis left;
bar(categorical(stageNames), bdtTotals);
ylabel('Estimated BDT (s)');
yyaxis right;
plot(1:numel(stageNames), activeLayers, 'o-', 'LineWidth', 1.5, 'MarkerSize', 6);
ylabel('Active layers');
title('Nominal pipeline BDT and active layers');
grid on;

nexttile;
bar(categorical(stageNames), activeSpots);
ylabel('Active spots / retained bixels');
title('Nominal pipeline active spot summary');
grid on;

plotPath = fullfile(outputDir, sprintf('%s_stage_summary.png', plotPrefix));
set(fig, 'PaperPositionMode', 'auto');
print(fig, plotPath, '-dpng', '-r150');
close(fig);
end

function [bdtTotal, nActiveLayers, nActiveSpots] = compute_stage_delivery_metrics(w, map, groups, cfg, machine)
w = reshape(double(w), [], 1);
seq = extract_delivery_sequence(w, map, groups, cfg);
metrics = estimate_bdt(seq, machine);
bdtTotal = double(metrics.T_total);
nActiveLayers = sum(build_layer_mask(w, groups, cfg.layerThreshold));
nActiveSpots = sum(w > cfg.spotThreshold);
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
