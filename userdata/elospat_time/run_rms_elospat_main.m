function outputs = run_rms_elospat_main(cfgOverride)
% run_rms_elospat_main  End-to-end RMS-ELO-SPAT orchestrator.
%
% This research-only entrypoint stitches together:
%   - IMPT baseline
%   - Dense Arc baseline
%   - ELO-SPAT-like layer filtering
%   - spot-sparsity stage
%   - post-filter baseline
%   - final RMS-ELO-SPAT re-optimization
%   - metrics + exports
%
% INPUT
%   cfgOverride - optional struct with field overrides applied on top of
%                 default_config(). A nested cfgOverride.precomputed struct
%                 may be used to inject synthetic or cached objects for
%                 smoke tests and sweeps.
%
% OUTPUT
%   outputs - struct with stage results, comparison table, and saved files.

narginchk(0, 1);

workspaceDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(workspaceDir));
configDir = fullfile(workspaceDir, 'config');
planningDir = fullfile(workspaceDir, 'planning');
mappingDir = fullfile(workspaceDir, 'mapping');
objectiveDir = fullfile(workspaceDir, 'objective');
deliveryDir = fullfile(workspaceDir, 'delivery');
integratedDir = fullfile(workspaceDir, 'integrated');
reductionDir = fullfile(workspaceDir, 'reduction');
robustnessDir = fullfile(workspaceDir, 'robustness');
evaluationDir = fullfile(workspaceDir, 'evaluation');
resultsDir = fullfile(workspaceDir, 'results');

addpath(configDir);
addpath(planningDir);
addpath(mappingDir);
addpath(objectiveDir);
addpath(deliveryDir);
addpath(integratedDir);
addpath(reductionDir);
addpath(robustnessDir);
addpath(evaluationDir);

researchMeta = research_only_metadata();

fprintf('\n========================================\n');
fprintf('RMS-ELO-SPAT Main Orchestrator\n');
fprintf('========================================\n');
fprintf(['Research-only workflow. Outputs are exploratory computational ' ...
         'artifacts and not clinical treatment plans.\n\n']);

cfg = default_config();
if nargin >= 1 && ~isempty(cfgOverride)
    if ~isstruct(cfgOverride)
        error('run_rms_elospat_main:invalidInput', ...
            'cfgOverride must be a struct when provided.');
    end
    cfg = merge_structs(cfg, cfgOverride);
end

runtimeCleanup = apply_runtime_mode(cfg); %#ok<NASGU>

if is_sweep_enabled(cfg)
    outputs = run_parameter_sweep_mode(cfg, repoRoot, resultsDir, researchMeta);
    return;
end

machine = default_machine_model();
precomputed = resolve_precomputed(cfg);

if isfield(cfg, 'randomSeed') && isnumeric(cfg.randomSeed) && isscalar(cfg.randomSeed) && cfg.randomSeed > 0
    rng(double(cfg.randomSeed), 'twister');
end

if isempty(precomputed)
    validate_config(cfg, machine);
    ensure_matrad_available();
end

if isfield(cfg, 'outputDir') && ~isempty(cfg.outputDir)
    outputDir = cfg.outputDir;
else
    outputDir = resultsDir;
end
if ~is_absolute_path(outputDir)
    outputDir = fullfile(repoRoot, outputDir);
end
if exist(outputDir, 'dir') ~= 7
    mkdir(outputDir);
end

timestamp = datestr(now, 'yyyy-mm-dd_HHMMSS');
prefix = sprintf('rms_elospat_main_%s', timestamp);

[ct, cst, patientFile] = resolve_patient_data(cfg, precomputed);
fprintf('Patient source: %s\n', patientFile);
fprintf('Arc spacing: %.2f deg\n', cfg.arcSpacing);
fprintf('Robust enabled: %s\n\n', mat2str(logical(cfg.robust.enabled)));

% ---------------------------------------------------------------------
% Stage 1: IMPT baseline
% ---------------------------------------------------------------------
fprintf('--- Stage 1: IMPT baseline ---\n');
[impt, imptMapContext] = resolve_or_run_impt_stage(ct, cst, cfg, precomputed);

% ---------------------------------------------------------------------
% Stage 2: Discrete arc + Dense Arc baseline
% ---------------------------------------------------------------------
fprintf('--- Stage 2: Dense Arc baseline ---\n');
[arc, arcMapContext] = resolve_or_run_arc_stage(ct, cst, cfg, precomputed);

% ---------------------------------------------------------------------
% Stage 3-6: Wrapper / baselines
% ---------------------------------------------------------------------
fprintf('--- Stage 3: ELO-SPAT-like stage ---\n');
[resultLayer, layerMask] = run_elospat_like_stage( ...
    arc.dij, cst, arc.pln, arc.denseResult.w, arcMapContext.groups, arcMapContext.cpGroups, cfg, machine); %#ok<NASGU>

fprintf('--- Stage 4: Spot sparse stage ---\n');
[resultSpot, spotMask] = run_spot_sparse_stage( ...
    arc.dij, cst, arc.pln, resultLayer.w, arcMapContext.map, arcMapContext.groups, cfg, machine); %#ok<NASGU>

fprintf('--- Stage 5: Post-filter baseline ---\n');
[postFilterResult, postFilterMask] = run_post_filter_baseline( ...
    arc.dij, cst, arc.pln, arc.denseResult.w, ...
    arcMapContext.map, arcMapContext.groups, cfg, machine); %#ok<NASGU>

fprintf('--- Stage 6: Final re-optimization ---\n');
% Freeze the support selected by the sparse stages.  Final re-optimization
% is deliberately performed on D(:,S), then embedded back into full space.
supportMask = logical(spotMask) & (abs(reshape(double(resultLayer.w), [], 1)) > cfg.spotThreshold);
resultFinal = run_final_reoptimization(arc.dij, cst, arc.pln, resultSpot.w, ...
    cfg, supportMask, layerMask);

% ---------------------------------------------------------------------
% Optional new integrated path. The heuristic wrapper above is retained
% and evaluated independently as a comparator.
% ---------------------------------------------------------------------
integratedRun = [];
if strcmp(cfg.mode, 'integrated')
    fprintf('--- Integrated MATLAB RMS-ELO-SPAT optimizer ---\n');
    integratedRun = run_integrated_rms_elospat( ...
        arc.dij, cst, arc.pln, arc.denseResult.w, ...
        arcMapContext.map, arcMapContext.groups, arcMapContext.cpGroups, ...
        cfg, machine);
end

% ---------------------------------------------------------------------
% Stage 7: Metrics and comparison
% ---------------------------------------------------------------------
fprintf('--- Stage 7: Metrics and exports ---\n');
Darc = get_nominal_dose_matrix(arc.dij);
Dimpt = get_nominal_dose_matrix(impt.dij);

methodEntries = repmat(struct( ...
    'name', '', ...
    'doseMetrics', struct(), ...
    'deliveryMetrics', struct(), ...
    'solverMetrics', struct()), 0, 1);

methodEntries(end + 1, 1) = build_method_entry('IMPT', impt.resultGUI.w, impt.resultGUI, impt.meta, [], Dimpt, ct, cst, cfg, machine, imptMapContext.map, imptMapContext.groups); %#ok<AGROW>
methodEntries(end + 1, 1) = build_method_entry('Dense Arc', arc.denseResult.w, arc.denseResult, arc.denseMeta, [], Darc, ct, cst, cfg, machine, arcMapContext.map, arcMapContext.groups); %#ok<AGROW>
methodEntries(end + 1, 1) = build_method_entry('ELO-SPAT-like', resultLayer.w, struct(), struct('statusMessage', 'Wrapper stage (solver summary not retained)'), [], Darc, ct, cst, cfg, machine, arcMapContext.map, arcMapContext.groups); %#ok<AGROW>
methodEntries(end + 1, 1) = build_method_entry('Post-filter', postFilterResult.w, postFilterResult.resultGUI, struct('statusMessage', 'Post-filter baseline'), postFilterResult.optimizer, Darc, ct, cst, cfg, machine, arcMapContext.map, arcMapContext.groups); %#ok<AGROW>
methodEntries(end + 1, 1) = build_method_entry('RMS-ELO-SPAT', resultFinal.w, resultFinal.resultGUI, struct('statusMessage', 'Legacy heuristic wrapper final stage'), resultFinal.optimizer, Darc, ct, cst, cfg, machine, arcMapContext.map, arcMapContext.groups); %#ok<AGROW>
if ~isempty(integratedRun)
    integratedMeta = struct('statusMessage', integratedRun.optimizer.status, ...
        'iterations', integratedRun.optimizer.iterations, ...
        'objective', integratedRun.optimizer.objective);
    methodEntries(end + 1, 1) = build_method_entry( ...
        'Integrated RMS-ELO-SPAT', integratedRun.final.w, ...
        integratedRun.final.resultGUI, integratedMeta, ...
        integratedRun.final.optimizer, Darc, ct, cst, cfg, machine, ...
        arcMapContext.map, arcMapContext.groups); %#ok<AGROW>
end

methodTable = compare_methods_table(methodEntries);

summaryPlot = fullfile(outputDir, sprintf('%s_summary.png', prefix));
save_orchestrator_summary_plot(methodEntries, methodTable, summaryPlot);

methodCsv = fullfile(outputDir, sprintf('%s_methods.csv', prefix));
writetable(methodTable, methodCsv);

latexFile = fullfile(outputDir, sprintf('%s_table.tex', prefix));
[~, latexFile] = export_latex_table(methodTable, latexFile);

tradeoffOutputs = export_tradeoff_data(methodTable, sprintf('%s_tradeoff', prefix), outputDir);

% Lightweight nominal artifact for later robust resume.
nominalArtifactFile = fullfile(outputDir, sprintf('%s_nominal_artifact.mat', prefix));
pln = arc.pln; %#ok<NASGU>
stf = arc.stf; %#ok<NASGU>
denseResult = arc.denseResult; %#ok<NASGU>
save(nominalArtifactFile, ...
    'cfg', 'machine', 'pln', 'stf', 'denseResult', 'resultLayer', 'resultSpot', ...
    'resultFinal', 'integratedRun', 'researchMeta', ...
    '-v7.3');

robustOutputs = struct('status', 'not_requested');
if isfield(cfg, 'robust') && isstruct(cfg.robust) && isfield(cfg.robust, 'enabled') && logical(cfg.robust.enabled)
    try
        robustOutputs = run_task28_30_from_nominal(nominalArtifactFile, outputDir, false);
        robustOutputs.status = 'completed';
    catch ME
        robustOutputs = struct( ...
            'status', 'failed', ...
            'errorMessage', ME.message);
    end
end

matFile = fullfile(outputDir, sprintf('%s.mat', prefix));
stageSummaries = struct();
stageSummaries.resultLayer = strip_large_fields(resultLayer);
stageSummaries.resultSpot = strip_large_fields(resultSpot);
stageSummaries.postFilterResult = strip_large_fields(postFilterResult);
stageSummaries.resultFinal = strip_large_fields(resultFinal);
stageSummaries.integratedRun = strip_large_fields(integratedRun);
stageSummaries.imptMeta = impt.meta;
stageSummaries.denseMeta = arc.denseMeta;

save(matFile, ...
    'cfg', 'machine', 'patientFile', ...
    'methodEntries', 'methodTable', 'tradeoffOutputs', 'robustOutputs', ...
    'stageSummaries', 'summaryPlot', 'methodCsv', 'latexFile', 'nominalArtifactFile', 'researchMeta', ...
    '-v7.3');

outputs = struct();
outputs.cfg = cfg;
outputs.machine = machine;
outputs.researchMeta = researchMeta;
outputs.patientFile = patientFile;
outputs.methodEntries = methodEntries;
outputs.methodTable = methodTable;
outputs.files = struct( ...
    'summaryPlot', summaryPlot, ...
    'methodCsv', methodCsv, ...
    'latexFile', latexFile, ...
    'tradeoffCsv', tradeoffOutputs.csvFile, ...
    'tradeoffMat', tradeoffOutputs.matFile, ...
    'nominalArtifactFile', nominalArtifactFile, ...
    'matFile', matFile);
outputs.stages = struct( ...
    'impt', impt, ...
    'arc', arc, ...
    'resultLayer', resultLayer, ...
    'resultSpot', resultSpot, ...
    'postFilterResult', postFilterResult, ...
    'resultFinal', resultFinal, ...
    'integratedRun', integratedRun);
outputs.robust = robustOutputs;

fprintf('\n========================================\n');
fprintf('ORCHESTRATOR COMPLETE\n');
fprintf('========================================\n');
fprintf('Saved summary MAT: %s\n', matFile);
fprintf('Saved comparison CSV: %s\n', methodCsv);
fprintf('Saved LaTeX table: %s\n', latexFile);
fprintf('Saved tradeoff CSV: %s\n', tradeoffOutputs.csvFile);
fprintf('Saved summary plot: %s\n\n', summaryPlot);
end

function outputs = run_parameter_sweep_mode(cfg, repoRoot, resultsDir, researchMeta)
fprintf('\n========================================\n');
fprintf('RMS-ELO-SPAT Parameter Sweep Mode\n');
fprintf('========================================\n');

sweep = cfg.sweep;
[paramNames, valueSets] = resolve_sweep_grid(sweep);

if isfield(sweep, 'outputDir') && ~isempty(sweep.outputDir)
    outputDir = sweep.outputDir;
else
    outputDir = fullfile(resultsDir, sprintf('task40_sweep_%s', datestr(now, 'yyyy-mm-dd_HHMMSS')));
end
if ~is_absolute_path(outputDir)
    outputDir = fullfile(repoRoot, outputDir);
end
if exist(outputDir, 'dir') ~= 7
    mkdir(outputDir);
end

combinations = build_sweep_combinations(valueSets);
nRuns = size(combinations, 1);
fprintf('Sweep parameters: %s\n', strjoin(paramNames, ', '));
fprintf('Total combinations: %d\n', nRuns);
fprintf('Output root: %s\n\n', outputDir);

records = repmat(struct( ...
    'runIndex', 0, ...
    'status', "", ...
    'errorMessage', "", ...
    'outputDir', "", ...
    'summaryMatFile', "", ...
    'paramValues', struct(), ...
    'denseBDT', NaN, ...
    'rmsBDT', NaN, ...
    'deliveryReductionPercent', NaN, ...
    'denseCTVd95', NaN, ...
    'rmsCTVd95', NaN), nRuns, 1);
failureMessages = strings(0, 1);

for runIx = 1:nRuns
    cfgRun = cfg;
    cfgRun.sweep.enabled = false;

    paramStruct = struct();
    labelParts = strings(1, numel(paramNames));
    for j = 1:numel(paramNames)
        value = valueSets{j}(combinations(runIx, j));
        value = value{1};
        cfgRun = set_nested_field(cfgRun, paramNames{j}, value);
        paramStruct.(matlab.lang.makeValidName(paramNames{j})) = value;
        labelParts(j) = sprintf('%s-%s', sanitize_label_token(paramNames{j}), sanitize_value_token(value));
    end

    runDir = fullfile(outputDir, sprintf('run_%02d_%s', runIx, strjoin(cellstr(labelParts), '__')));
    if exist(runDir, 'dir') ~= 7
        mkdir(runDir);
    end
    cfgRun.outputDir = runDir;

    records(runIx).runIndex = runIx;
    records(runIx).paramValues = paramStruct;
    records(runIx).outputDir = string(runDir);

    fprintf('Sweep run %d/%d -> %s\n', runIx, nRuns, runDir);
    try
        runOutputs = run_rms_elospat_main(cfgRun);
        [denseRow, rmsRow] = extract_dense_and_rms_rows(runOutputs.methodTable);
        reductionPercent = 100 * (double(denseRow.delivery_T_total) - double(rmsRow.delivery_T_total)) / ...
            max(double(denseRow.delivery_T_total), eps);

        records(runIx).status = "completed";
        records(runIx).summaryMatFile = string(runOutputs.files.matFile);
        records(runIx).denseBDT = double(denseRow.delivery_T_total);
        records(runIx).rmsBDT = double(rmsRow.delivery_T_total);
        records(runIx).deliveryReductionPercent = reductionPercent;
        records(runIx).denseCTVd95 = double(denseRow.dose_targetD95Min);
        records(runIx).rmsCTVd95 = double(rmsRow.dose_targetD95Min);
    catch ME
        records(runIx).status = "failed";
        records(runIx).errorMessage = string(ME.message);
        failureMessages(end + 1, 1) = sprintf('Run %d failed: %s', runIx, ME.message); %#ok<AGROW>
        fprintf(2, '  FAILED: %s\n', ME.message);
    end
end

summaryTable = records_to_table(records, paramNames);
summaryMatFile = fullfile(outputDir, 'task40_sweep_summary.mat');
summaryCsvFile = fullfile(outputDir, 'task40_sweep_summary.csv');
failureLogFile = fullfile(outputDir, 'task40_sweep_failures.log');
heatmapFile = fullfile(outputDir, 'task40_sweep_heatmap.png');

writetable(summaryTable, summaryCsvFile);
save(summaryMatFile, 'records', 'summaryTable', 'paramNames', 'valueSets', 'researchMeta', '-v7.3');
write_failure_log(failureLogFile, failureMessages, researchMeta);
save_sweep_heatmap(summaryTable, paramNames, valueSets, heatmapFile);

outputs = struct();
outputs.mode = 'sweep';
outputs.researchMeta = researchMeta;
outputs.outputDir = outputDir;
outputs.summaryTable = summaryTable;
outputs.files = struct( ...
    'summaryMatFile', summaryMatFile, ...
    'summaryCsvFile', summaryCsvFile, ...
    'failureLogFile', failureLogFile, ...
    'heatmapFile', heatmapFile);

fprintf('\nSweep complete.\n');
fprintf('Saved summary MAT: %s\n', summaryMatFile);
fprintf('Saved summary CSV: %s\n', summaryCsvFile);
fprintf('Saved failure log: %s\n', failureLogFile);
fprintf('Saved heatmap: %s\n\n', heatmapFile);
end

function precomputed = resolve_precomputed(cfg)
precomputed = [];
if isfield(cfg, 'precomputed') && isstruct(cfg.precomputed)
    precomputed = cfg.precomputed;
end
end

function ensure_matrad_available()
requiredFns = {'matRad_generateStf', 'matRad_calcDoseInfluence', 'matRad_fluenceOptimization'};
for i = 1:numel(requiredFns)
    if exist(requiredFns{i}, 'file') ~= 2
        error('run_rms_elospat_main:missingMatRad', ...
            'Required matRad function %s is not on the MATLAB path.', requiredFns{i});
    end
end
end

function tf = is_sweep_enabled(cfg)
tf = isfield(cfg, 'sweep') && isstruct(cfg.sweep) && ...
    isfield(cfg.sweep, 'enabled') && logical(cfg.sweep.enabled);
end

function cleanupObj = apply_runtime_mode(cfg)
state = struct();
state.defaultFigureVisible = get(groot, 'DefaultFigureVisible');
state.hasMatRadConfig = false;
state.previousDisableGUI = [];
state.appliedSoftwareOpenGL = false;

runtimeCfg = struct();
if isfield(cfg, 'runtime') && isstruct(cfg.runtime)
    runtimeCfg = cfg.runtime;
end

headless = isfield(runtimeCfg, 'headless') && logical(runtimeCfg.headless);
forceSoftwareOpenGL = isfield(runtimeCfg, 'forceSoftwareOpenGL') && logical(runtimeCfg.forceSoftwareOpenGL);

if headless
    set(groot, 'DefaultFigureVisible', 'off');
    if exist('MatRad_Config', 'class') == 8
        try
            matRad_cfg = MatRad_Config.instance();
            state.hasMatRadConfig = true;
            state.previousDisableGUI = matRad_cfg.disableGUI;
            matRad_cfg.disableGUI = true;
        catch
            state.hasMatRadConfig = false;
        end
    end
end

if forceSoftwareOpenGL
    try
        opengl('software');
        state.appliedSoftwareOpenGL = true;
    catch
        state.appliedSoftwareOpenGL = false;
    end
end

cleanupObj = onCleanup(@() restore_runtime_mode(state));
end

function restore_runtime_mode(state)
try
    set(groot, 'DefaultFigureVisible', state.defaultFigureVisible);
catch
end

if isfield(state, 'hasMatRadConfig') && state.hasMatRadConfig
    try
        matRad_cfg = MatRad_Config.instance();
        matRad_cfg.disableGUI = state.previousDisableGUI;
    catch
    end
end
end

function [paramNames, valueSets] = resolve_sweep_grid(sweep)
requiredFields = {'parameters', 'values'};
for i = 1:numel(requiredFields)
    if ~isfield(sweep, requiredFields{i})
        error('run_rms_elospat_main:invalidSweep', ...
            'cfg.sweep.%s is required for sweep mode.', requiredFields{i});
    end
end

paramNames = cellstr(string(sweep.parameters(:)));
valueSets = sweep.values;
if ~iscell(valueSets) || numel(valueSets) ~= numel(paramNames)
    error('run_rms_elospat_main:invalidSweep', ...
        'cfg.sweep.values must be a cell array matching cfg.sweep.parameters.');
end

for i = 1:numel(valueSets)
    current = valueSets{i};
    if iscell(current)
        boxed = current(:);
    else
        boxed = num2cell(current(:));
    end
    if isempty(boxed)
        error('run_rms_elospat_main:invalidSweep', ...
            'cfg.sweep.values{%d} must contain at least one value.', i);
    end
    valueSets{i} = boxed;
end
end

function combinations = build_sweep_combinations(valueSets)
nParams = numel(valueSets);
indexRanges = cell(1, nParams);
for i = 1:nParams
    indexRanges{i} = 1:numel(valueSets{i});
end

gridCells = cell(1, nParams);
[gridCells{:}] = ndgrid(indexRanges{:});
combinations = zeros(numel(gridCells{1}), nParams);
for i = 1:nParams
    combinations(:, i) = gridCells{i}(:);
end
end

function [ct, cst, patientFile] = resolve_patient_data(cfg, precomputed)
if ~isempty(precomputed) && isfield(precomputed, 'ct') && isfield(precomputed, 'cst')
    ct = precomputed.ct;
    cst = precomputed.cst;
    patientFile = 'precomputed';
    return;
end

patientFile = cfg.patientFile;
loaded = load(patientFile, 'ct', 'cst');
ct = loaded.ct;
cst = loaded.cst;
end

function [impt, mapContext] = resolve_or_run_impt_stage(ct, cst, cfg, precomputed)
if ~isempty(precomputed) && isfield(precomputed, 'impt')
    impt = precomputed.impt;
else
    [impt.pln, impt.stf, impt.dij, impt.resultGUI, impt.meta] = build_impt_plan(ct, cst, cfg);
end

if isfield(impt, 'map') && isfield(impt, 'groups') && isfield(impt, 'cpGroups')
    mapContext.map = impt.map;
    mapContext.groups = impt.groups;
    mapContext.cpGroups = impt.cpGroups;
else
    mapContext.map = map_bixel_to_cp_energy_spot(impt.dij, impt.stf, impt.pln);
    [mapContext.groups, mapContext.cpGroups] = build_group_index(mapContext.map);
end
end

function [arc, mapContext] = resolve_or_run_arc_stage(ct, cst, cfg, precomputed)
if ~isempty(precomputed) && isfield(precomputed, 'arc')
    arc = precomputed.arc;
else
    [arc.pln, arc.stf, arc.dij] = build_discrete_arc_plan(ct, cst, cfg);
    [arc.denseResult, arc.denseMeta] = run_dense_arc_baseline(arc.dij, cst, arc.pln);
end

if ~isfield(arc, 'denseResult') || ~isfield(arc, 'denseMeta')
    [arc.denseResult, arc.denseMeta] = run_dense_arc_baseline(arc.dij, cst, arc.pln);
end

if isfield(arc, 'map') && isfield(arc, 'groups') && isfield(arc, 'cpGroups')
    mapContext.map = arc.map;
    mapContext.groups = arc.groups;
    mapContext.cpGroups = arc.cpGroups;
else
    mapContext.map = map_bixel_to_cp_energy_spot(arc.dij, arc.stf, arc.pln);
    [mapContext.groups, mapContext.cpGroups] = build_group_index(mapContext.map);
end

arc.map = mapContext.map;
arc.groups = mapContext.groups;
arc.cpGroups = mapContext.cpGroups;
end

function D = get_nominal_dose_matrix(dij)
if iscell(dij.physicalDose)
    D = dij.physicalDose{1};
else
    D = dij.physicalDose;
end
end

function entry = build_method_entry(name, w, resultLike, meta, optimizer, D, ct, cst, cfg, machine, map, groups)
entry = struct();
entry.name = name;
entry.doseMetrics = compute_dvh_metrics(full(D * reshape(double(w), [], 1)), ct, cst, cfg);
entry.deliveryMetrics = compute_delivery_metrics(w, map, groups, cfg, machine);
entry.solverMetrics = compute_solver_metrics(resultLike, meta, optimizer);
end

function save_orchestrator_summary_plot(methodEntries, methodTable, plotFile)
targetIx = find_first_target_index(methodEntries(1).doseMetrics);
fig = figure('Visible', 'off', 'Color', 'w');
tiledlayout(2, 1);

nexttile;
hold on;
for i = 1:numel(methodEntries)
    dm = methodEntries(i).doseMetrics;
    if isempty(targetIx) || targetIx > numel(dm.structures)
        continue;
    end
    plot(dm.doseGrid, dm.structures(targetIx).dvhCurve, 'LineWidth', 1.5);
end
legend(string(methodTable.Method), 'Interpreter', 'none', 'Location', 'bestoutside');
xlabel('Dose');
ylabel('Volume (%)');
title('Target DVH summary across methods');
grid on;

nexttile;
yyaxis left;
bar(categorical(string(methodTable.Method)), methodTable.dose_targetD95Min);
ylabel('CTV D95');
yyaxis right;
plot(1:height(methodTable), methodTable.delivery_T_total, 'o-', 'LineWidth', 1.5);
ylabel('Estimated BDT (s)');
title('Dose and delivery summary');
grid on;

set(fig, 'PaperPositionMode', 'auto');
print(fig, plotFile, '-dpng', '-r150');
close(fig);
end

function idx = find_first_target_index(doseMetrics)
idx = [];
if ~isstruct(doseMetrics) || ~isfield(doseMetrics, 'structures')
    return;
end
for i = 1:numel(doseMetrics.structures)
    if isfield(doseMetrics.structures(i), 'isTarget') && logical(doseMetrics.structures(i).isTarget)
        idx = i;
        return;
    end
end
end

function stripped = strip_large_fields(resultStruct)
stripped = resultStruct;
if ~isstruct(stripped)
    return;
end
if isfield(stripped, 'resultGUI') && isstruct(stripped.resultGUI)
    stripped.resultGUI = keep_solver_and_weight_fields(stripped.resultGUI);
end
end

function s = keep_solver_and_weight_fields(s)
keepFields = {'w', 'fval', 'info', 'usedOptimizer', 'converged', 'exitStatus', 'statusMessage'};
fields = fieldnames(s);
for i = 1:numel(fields)
    if ~ismember(fields{i}, keepFields)
        s = rmfield(s, fields{i});
    end
end
end

function merged = merge_structs(base, override)
merged = base;
overrideFields = fieldnames(override);
for i = 1:numel(overrideFields)
    key = overrideFields{i};
    if isfield(base, key) && isstruct(base.(key)) && isstruct(override.(key))
        merged.(key) = merge_structs(base.(key), override.(key));
    else
        merged.(key) = override.(key);
    end
end
end

function s = set_nested_field(s, dottedName, value)
parts = strsplit(char(dottedName), '.');
if numel(parts) == 1
    s.(parts{1}) = value;
    return;
end

fieldName = parts{1};
if ~isfield(s, fieldName) || ~isstruct(s.(fieldName))
    s.(fieldName) = struct();
end
s.(fieldName) = set_nested_field(s.(fieldName), strjoin(parts(2:end), '.'), value);
end

function token = sanitize_label_token(value)
token = regexprep(char(string(value)), '[^A-Za-z0-9]+', '_');
token = regexprep(token, '^_+|_+$', '');
if isempty(token)
    token = 'param';
end
end

function token = sanitize_value_token(value)
if islogical(value)
    token = lower(mat2str(value));
    return;
end
if isnumeric(value) && isscalar(value)
    token = regexprep(sprintf('%.6g', double(value)), '[^A-Za-z0-9]+', '_');
else
    token = regexprep(char(string(value)), '[^A-Za-z0-9]+', '_');
end
token = regexprep(token, '^_+|_+$', '');
if isempty(token)
    token = 'value';
end
end

function [denseRow, rmsRow] = extract_dense_and_rms_rows(methodTable)
denseMask = strcmp(string(methodTable.Method), "Dense Arc");
rmsMask = strcmp(string(methodTable.Method), "RMS-ELO-SPAT");
if ~any(denseMask) || ~any(rmsMask)
    error('run_rms_elospat_main:missingMethods', ...
        'Method table must contain Dense Arc and RMS-ELO-SPAT rows.');
end
denseRow = methodTable(find(denseMask, 1, 'first'), :);
rmsRow = methodTable(find(rmsMask, 1, 'first'), :);
end

function summaryTable = records_to_table(records, paramNames)
nRuns = numel(records);
summaryTable = table((1:nRuns)', strings(nRuns, 1), strings(nRuns, 1), strings(nRuns, 1), ...
    NaN(nRuns, 1), NaN(nRuns, 1), NaN(nRuns, 1), NaN(nRuns, 1), NaN(nRuns, 1), ...
    'VariableNames', {'runIndex', 'status', 'errorMessage', 'outputDir', ...
    'denseBDT', 'rmsBDT', 'deliveryReductionPercent', 'denseCTVd95', 'rmsCTVd95'});

for i = 1:nRuns
    summaryTable.status(i) = string(records(i).status);
    summaryTable.errorMessage(i) = string(records(i).errorMessage);
    summaryTable.outputDir(i) = string(records(i).outputDir);
    summaryTable.denseBDT(i) = records(i).denseBDT;
    summaryTable.rmsBDT(i) = records(i).rmsBDT;
    summaryTable.deliveryReductionPercent(i) = records(i).deliveryReductionPercent;
    summaryTable.denseCTVd95(i) = records(i).denseCTVd95;
    summaryTable.rmsCTVd95(i) = records(i).rmsCTVd95;
end

for j = 1:numel(paramNames)
    columnName = matlab.lang.makeValidName(paramNames{j});
    column = strings(nRuns, 1);
    numericColumn = NaN(nRuns, 1);
    allNumeric = true;
    for i = 1:nRuns
        value = records(i).paramValues.(columnName);
        if isnumeric(value) && isscalar(value)
            numericColumn(i) = double(value);
            column(i) = string(value);
        elseif islogical(value) && isscalar(value)
            numericColumn(i) = double(value);
            column(i) = string(value);
        else
            allNumeric = false;
            column(i) = string(value);
        end
    end
    if allNumeric
        summaryTable.(columnName) = numericColumn;
    else
        summaryTable.(columnName) = column;
    end
end
end

function write_failure_log(failureLogFile, failureMessages, researchMeta)
fid = fopen(failureLogFile, 'w');
if fid == -1
    error('run_rms_elospat_main:fileOpenFailed', ...
        'Could not open %s for writing.', failureLogFile);
end
cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid, '%s\n', researchMeta.disclaimer);
fprintf(fid, 'Machine model note: %s\n', researchMeta.machineModelNote);
fprintf(fid, 'Citation placeholder: %s\n\n', researchMeta.citationPlaceholder);
if isempty(failureMessages)
    fprintf(fid, 'No sweep failures were recorded.\n');
else
    for i = 1:numel(failureMessages)
        fprintf(fid, '%s\n', failureMessages(i));
    end
end
end

function save_sweep_heatmap(summaryTable, paramNames, valueSets, plotFile)
fig = figure('Visible', 'off', 'Color', 'w');

if numel(paramNames) ~= 2 || ~are_numeric_value_sets(valueSets)
    text(0.1, 0.5, 'Heatmap available only for 2-parameter sweeps.', 'FontSize', 11);
    axis off;
else
    xName = matlab.lang.makeValidName(paramNames{1});
    yName = matlab.lang.makeValidName(paramNames{2});
    xValues = unique(cellfun(@double, valueSets{1}));
    yValues = unique(cellfun(@double, valueSets{2}));
    heatmapValues = NaN(numel(yValues), numel(xValues));

    for i = 1:height(summaryTable)
        if summaryTable.status(i) ~= "completed"
            continue;
        end
        [~, xIx] = min(abs(xValues - double(summaryTable.(xName)(i))));
        [~, yIx] = min(abs(yValues - double(summaryTable.(yName)(i))));
        heatmapValues(yIx, xIx) = summaryTable.deliveryReductionPercent(i);
    end

    imagesc(xValues, yValues, heatmapValues);
    set(gca, 'YDir', 'normal');
    colorbar;
    xlabel(paramNames{1}, 'Interpreter', 'none');
    ylabel(paramNames{2}, 'Interpreter', 'none');
    title('Task 40 parameter-sweep heatmap (BDT reduction %)');
end

set(fig, 'PaperPositionMode', 'auto');
print(fig, plotFile, '-dpng', '-r150');
close(fig);
end

function tf = are_numeric_value_sets(valueSets)
tf = true;
for i = 1:numel(valueSets)
    current = valueSets{i};
    for j = 1:numel(current)
        value = current{j};
        if ~(isnumeric(value) && isscalar(value)) && ~(islogical(value) && isscalar(value))
            tf = false;
            return;
        end
    end
end
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
