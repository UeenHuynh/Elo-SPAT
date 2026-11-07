function outputs = run_task31_33_from_artifact(artifactFile, patientFile, outputDir)
% run_task31_33_from_artifact  Compute Tasks 31-33 outputs from a saved nominal artifact.
%
% This helper reuses an existing Task 25-27 nominal artifact and computes:
%   Task 31 - delivery metrics
%   Task 32 - solver metrics
%   Task 33 - comparison table
%
% It saves:
%   - a MAT file with methodEntries and methodTable
%   - a CSV export of the comparison table
%   - a PNG summary bar chart

narginchk(0, 3);

planningDir = fileparts(mfilename('fullpath'));
workspaceDir = fileparts(planningDir);
configDir = fullfile(workspaceDir, 'config');
deliveryDir = fullfile(workspaceDir, 'delivery');
evaluationDir = fullfile(workspaceDir, 'evaluation');
resultsDir = fullfile(workspaceDir, 'results');

addpath(configDir);
addpath(deliveryDir);
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
fprintf('Task 31-33 Evaluation From Artifact\n');
fprintf('========================================\n');
fprintf('Artifact: %s\n', artifactFile);
fprintf('Patient:  %s\n', patientFile);
fprintf('Output:   %s\n\n', outputDir);

data = load(artifactFile);
patientData = load(patientFile, 'ct', 'cst');

cfg = resolve_cfg(data);
machine = resolve_machine(data);
nominalDoseMatrix = get_nominal_dose_matrix(data.dij);

methodEntries = repmat(struct( ...
    'name', '', ...
    'doseMetrics', struct(), ...
    'deliveryMetrics', struct(), ...
    'solverMetrics', struct()), 0, 1);

methodEntries(end + 1, 1) = build_dense_entry(data, patientData, cfg, machine, nominalDoseMatrix); %#ok<AGROW>
methodEntries(end + 1, 1) = build_layer_entry(data, patientData, cfg, machine, nominalDoseMatrix); %#ok<AGROW>
methodEntries(end + 1, 1) = build_spot_entry(data, patientData, cfg, machine, nominalDoseMatrix); %#ok<AGROW>
methodEntries(end + 1, 1) = build_final_entry(data, patientData, cfg, machine, nominalDoseMatrix); %#ok<AGROW>

methodTable = compare_methods_table(methodEntries);

[~, artifactBaseName] = fileparts(artifactFile);
outputPrefix = sprintf('%s_task31_33', artifactBaseName);
matFile = fullfile(outputDir, sprintf('%s.mat', outputPrefix));
csvFile = fullfile(outputDir, sprintf('%s.csv', outputPrefix));
plotFile = fullfile(outputDir, sprintf('%s_summary.png', outputPrefix));

writetable(methodTable, csvFile);
save(matFile, 'methodEntries', 'methodTable', 'artifactFile', 'patientFile', '-v7.3');
save_summary_plot(methodTable, plotFile);

outputs = struct();
outputs.matFile = matFile;
outputs.csvFile = csvFile;
outputs.plotFile = plotFile;
outputs.methodTable = methodTable;

fprintf('Saved MAT:  %s\n', matFile);
fprintf('Saved CSV:  %s\n', csvFile);
fprintf('Saved plot: %s\n\n', plotFile);
end

function artifactFile = resolve_latest_artifact(resultsDir)
listing = dir(fullfile(resultsDir, 'task25_27_pipeline_*.mat'));
if isempty(listing)
    error('run_task31_33_from_artifact:missingArtifact', ...
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
    error('run_task31_33_from_artifact:missingPatient', ...
        'PROSTATE.mat was not found and matRad_rc is not on the MATLAB path.');
end

matRadRoot = fileparts(matRadRcPath);
candidate = fullfile(matRadRoot, 'data', 'PROSTATE.mat');
if exist(candidate, 'file') ~= 2
    error('run_task31_33_from_artifact:missingPatient', ...
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

function D = get_nominal_dose_matrix(dij)
if ~isstruct(dij) || ~isfield(dij, 'physicalDose')
    error('run_task31_33_from_artifact:missingDoseMatrix', ...
        'Artifact dij.physicalDose is required.');
end
if iscell(dij.physicalDose)
    D = dij.physicalDose{1};
else
    D = dij.physicalDose;
end
end

function entry = build_dense_entry(data, patientData, cfg, machine, D)
entry = struct();
entry.name = 'Dense Arc';
entry.doseMetrics = compute_stage_dose_metrics(data.denseResult.w, D, patientData, cfg);
entry.deliveryMetrics = compute_delivery_metrics(data.denseResult.w, data.map, data.groups, cfg, machine);
entry.solverMetrics = compute_solver_metrics(data.denseResult, data.denseMeta);
end

function entry = build_layer_entry(data, patientData, cfg, machine, D)
entry = struct();
entry.name = 'ELO-SPAT-like';
entry.doseMetrics = compute_stage_dose_metrics(data.resultLayer.w, D, patientData, cfg);
entry.deliveryMetrics = compute_delivery_metrics(data.resultLayer.w, data.map, data.groups, cfg, machine);
solverSeed = struct('statusMessage', 'Wrapper stage (solver summary not retained)');
entry.solverMetrics = compute_solver_metrics(struct(), solverSeed);
end

function entry = build_spot_entry(data, patientData, cfg, machine, D)
entry = struct();
entry.name = 'Spot Sparse';
entry.doseMetrics = compute_stage_dose_metrics(data.resultSpot.w, D, patientData, cfg);
entry.deliveryMetrics = compute_delivery_metrics(data.resultSpot.w, data.map, data.groups, cfg, machine);
solverSeed = struct('statusMessage', 'Wrapper stage (solver summary not retained)');
entry.solverMetrics = compute_solver_metrics(struct(), solverSeed);
end

function entry = build_final_entry(data, patientData, cfg, machine, D)
entry = struct();
entry.name = 'RMS-ELO-SPAT';
entry.doseMetrics = compute_stage_dose_metrics(data.resultFinal.w, D, patientData, cfg);
entry.deliveryMetrics = compute_delivery_metrics(data.resultFinal.w, data.map, data.groups, cfg, machine);
entry.solverMetrics = compute_solver_metrics(data.resultFinal.resultGUI, struct(), data.resultFinal.optimizer);
end

function doseMetrics = compute_stage_dose_metrics(w, D, patientData, cfg)
doseVector = full(D * reshape(double(w), [], 1));
doseMetrics = compute_dvh_metrics(doseVector, patientData.ct, patientData.cst, cfg);
end

function save_summary_plot(methodTable, plotFile)
methodLabels = string(methodTable.Method);

doseMetric = methodTable.dose_targetD95Min;
deliveryMetric = methodTable.delivery_T_total;
solverMetric = methodTable.solver_iterations;

fig = figure('Visible', 'off', 'Color', 'w');
tiledlayout(3, 1);

nexttile;
bar(categorical(methodLabels), doseMetric);
ylabel('CTV D95');
title('Task 31-33 dose summary');
grid on;

nexttile;
bar(categorical(methodLabels), deliveryMetric);
ylabel('Estimated BDT (s)');
title('Task 31-33 delivery summary');
grid on;

nexttile;
bar(categorical(methodLabels), solverMetric);
ylabel('Solver iterations');
xlabel('Method');
title('Task 31-33 solver summary');
grid on;

set(fig, 'PaperPositionMode', 'auto');
print(fig, plotFile, '-dpng', '-r150');
close(fig);
end
