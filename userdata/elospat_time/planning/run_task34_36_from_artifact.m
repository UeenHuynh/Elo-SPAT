function outputs = run_task34_36_from_artifact(artifactFile, patientFile, outputDir)
% run_task34_36_from_artifact  Run Tasks 34-36 from an existing nominal artifact.
%
% This helper:
%   - runs Task 36 post-filter baseline from the saved Dense Arc solution
%   - rebuilds a 5-method comparison table
%   - exports Task 34 LaTeX table
%   - exports Task 35 tradeoff CSV/MAT
%   - saves a simple tradeoff scatter PNG

narginchk(0, 3);

planningDir = fileparts(mfilename('fullpath'));
workspaceDir = fileparts(planningDir);
configDir = fullfile(workspaceDir, 'config');
deliveryDir = fullfile(workspaceDir, 'delivery');
evaluationDir = fullfile(workspaceDir, 'evaluation');
reductionDir = fullfile(workspaceDir, 'reduction');
resultsDir = fullfile(workspaceDir, 'results');

addpath(configDir);
addpath(deliveryDir);
addpath(evaluationDir);
addpath(reductionDir);

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
fprintf('Task 34-36 From Artifact\n');
fprintf('========================================\n');
fprintf('Artifact: %s\n', artifactFile);
fprintf('Patient:  %s\n', patientFile);
fprintf('Output:   %s\n\n', outputDir);

data = load(artifactFile);
patientData = load(patientFile, 'ct', 'cst');

cfg = resolve_cfg(data);
machine = resolve_machine(data);
D = get_nominal_dose_matrix(data.dij);

[postFilterResult, postFilterMask] = run_post_filter_baseline( ...
    data.dij, patientData.cst, data.pln, data.denseResult.w, ...
    data.map, data.groups, cfg, machine);

methodEntries = repmat(struct( ...
    'name', '', ...
    'doseMetrics', struct(), ...
    'deliveryMetrics', struct(), ...
    'solverMetrics', struct()), 0, 1);

methodEntries(end + 1, 1) = build_entry('Dense Arc', data.denseResult.w, data.denseResult, data.denseMeta, [], D, patientData, cfg, machine, data.map, data.groups); %#ok<AGROW>
methodEntries(end + 1, 1) = build_entry('ELO-SPAT-like', data.resultLayer.w, struct(), struct('statusMessage', 'Wrapper stage (solver summary not retained)'), [], D, patientData, cfg, machine, data.map, data.groups); %#ok<AGROW>
methodEntries(end + 1, 1) = build_entry('Spot Sparse', data.resultSpot.w, struct(), struct('statusMessage', 'Wrapper stage (solver summary not retained)'), [], D, patientData, cfg, machine, data.map, data.groups); %#ok<AGROW>
methodEntries(end + 1, 1) = build_entry('Post-filter', postFilterResult.w, postFilterResult.resultGUI, struct(), postFilterResult.optimizer, D, patientData, cfg, machine, data.map, data.groups); %#ok<AGROW>
methodEntries(end + 1, 1) = build_entry('RMS-ELO-SPAT', data.resultFinal.w, data.resultFinal.resultGUI, struct(), data.resultFinal.optimizer, D, patientData, cfg, machine, data.map, data.groups); %#ok<AGROW>

methodTable = compare_methods_table(methodEntries);

[~, artifactBaseName] = fileparts(artifactFile);
outputPrefix = sprintf('%s_task34_36', artifactBaseName);

matFile = fullfile(outputDir, sprintf('%s.mat', outputPrefix));
csvFile = fullfile(outputDir, sprintf('%s_methods.csv', outputPrefix));
texFile = fullfile(outputDir, sprintf('%s_table.tex', outputPrefix));
tradeoffPrefix = sprintf('%s_tradeoff', outputPrefix);
scatterFile = fullfile(outputDir, sprintf('%s_scatter.png', outputPrefix));

writetable(methodTable, csvFile);
[latexText, texFile] = export_latex_table(methodTable, texFile); %#ok<ASGLU>
tradeoffOutputs = export_tradeoff_data(methodTable, tradeoffPrefix, outputDir);
save_tradeoff_scatter(methodTable, scatterFile);

save(matFile, ...
    'artifactFile', 'patientFile', 'cfg', 'machine', ...
    'postFilterResult', 'postFilterMask', ...
    'methodEntries', 'methodTable', 'tradeoffOutputs', ...
    '-v7.3');

outputs = struct();
outputs.matFile = matFile;
outputs.csvFile = csvFile;
outputs.texFile = texFile;
outputs.tradeoffCsvFile = tradeoffOutputs.csvFile;
outputs.tradeoffMatFile = tradeoffOutputs.matFile;
outputs.scatterFile = scatterFile;

fprintf('Saved MAT:           %s\n', matFile);
fprintf('Saved method CSV:    %s\n', csvFile);
fprintf('Saved LaTeX table:   %s\n', texFile);
fprintf('Saved tradeoff CSV:  %s\n', tradeoffOutputs.csvFile);
fprintf('Saved tradeoff MAT:  %s\n', tradeoffOutputs.matFile);
fprintf('Saved tradeoff plot: %s\n\n', scatterFile);
end

function artifactFile = resolve_latest_artifact(resultsDir)
listing = dir(fullfile(resultsDir, 'task25_27_pipeline_*.mat'));
if isempty(listing)
    error('run_task34_36_from_artifact:missingArtifact', ...
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
    error('run_task34_36_from_artifact:missingPatient', ...
        'PROSTATE.mat was not found and matRad_rc is not on the MATLAB path.');
end

matRadRoot = fileparts(matRadRcPath);
candidate = fullfile(matRadRoot, 'data', 'PROSTATE.mat');
if exist(candidate, 'file') ~= 2
    error('run_task34_36_from_artifact:missingPatient', ...
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
    error('run_task34_36_from_artifact:missingDoseMatrix', ...
        'Artifact dij.physicalDose is required.');
end
if iscell(dij.physicalDose)
    D = dij.physicalDose{1};
else
    D = dij.physicalDose;
end
end

function entry = build_entry(name, w, resultGUI, meta, optimizer, D, patientData, cfg, machine, map, groups)
entry = struct();
entry.name = name;
entry.doseMetrics = compute_stage_dose_metrics(w, D, patientData, cfg);
entry.deliveryMetrics = compute_delivery_metrics(w, map, groups, cfg, machine);
entry.solverMetrics = compute_solver_metrics(resultGUI, meta, optimizer);
end

function doseMetrics = compute_stage_dose_metrics(w, D, patientData, cfg)
doseVector = full(D * reshape(double(w), [], 1));
doseMetrics = compute_dvh_metrics(doseVector, patientData.ct, patientData.cst, cfg);
end

function save_tradeoff_scatter(methodTable, plotFile)
fig = figure('Visible', 'off', 'Color', 'w');
scatter(methodTable.delivery_T_total, methodTable.dose_targetD95Min, 70, 'filled');
grid on;
xlabel('Estimated BDT (s)');
ylabel('CTV D95');
title('Task 35 tradeoff scatter');
text(methodTable.delivery_T_total + 0.5, methodTable.dose_targetD95Min, string(methodTable.Method), ...
    'FontSize', 9, 'Interpreter', 'none');
set(fig, 'PaperPositionMode', 'auto');
print(fig, plotFile, '-dpng', '-r150');
close(fig);
end
