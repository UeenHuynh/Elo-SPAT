function outputs = run_task28_30_from_nominal(nominalArtifactFile, outputDir, saveFullDij)
% run_task28_30_from_nominal  Resume Tasks 28-30 from a saved nominal 25-27 artifact.
%
% INPUTS
%   nominalArtifactFile - optional path to a saved task25_27_pipeline_*.mat artifact.
%                         If omitted, the newest matching artifact in ../results is used.
%   outputDir           - optional directory for resumed robust outputs.
%   saveFullDij         - optional logical flag. Default false to avoid saving
%                         multi-GB robust dij artifacts unless explicitly requested.
%
% OUTPUT
%   outputs - struct with saved file paths and key summary fields.

narginchk(0, 3);

fprintf('\n========================================\n');
fprintf('Task 28-30: Resume Robust From Nominal\n');
fprintf('========================================\n\n');

planningDir = fileparts(mfilename('fullpath'));
workspaceDir = fileparts(planningDir);
configDir = fullfile(workspaceDir, 'config');
mappingDir = fullfile(workspaceDir, 'mapping');
deliveryDir = fullfile(workspaceDir, 'delivery');
reductionDir = fullfile(workspaceDir, 'reduction');
objectiveDir = fullfile(workspaceDir, 'objective');
evaluationDir = fullfile(workspaceDir, 'evaluation');
robustnessDir = fullfile(workspaceDir, 'robustness');
resultsDir = fullfile(workspaceDir, 'results');

addpath(configDir);
addpath(mappingDir);
addpath(deliveryDir);
addpath(reductionDir);
addpath(objectiveDir);
addpath(evaluationDir);
addpath(robustnessDir);

if nargin < 1 || isempty(nominalArtifactFile)
    nominalArtifactFile = resolve_latest_nominal_artifact(resultsDir);
end

if nargin < 2 || isempty(outputDir)
    outputDir = resultsDir;
end

if nargin < 3 || isempty(saveFullDij)
    saveFullDij = false;
end

if exist(outputDir, 'dir') ~= 7
    mkdir(outputDir);
end

validateattributes(saveFullDij, {'logical', 'numeric'}, {'scalar'}, ...
    mfilename, 'saveFullDij');
saveFullDij = logical(saveFullDij);

fprintf('Nominal artifact: %s\n', nominalArtifactFile);
fprintf('Output directory: %s\n', outputDir);
fprintf('Save full robust dij: %s\n\n', mat2str(saveFullDij));

patientFile = resolve_patient_file();
fprintf('Patient file: %s\n\n', patientFile);

artifactData = load(nominalArtifactFile, ...
    'cfg', 'pln', 'stf', 'resultSpot', 'resultFinal', 'resultLayer', 'denseResult', 'machine');
load(patientFile, 'ct', 'cst');

cfg = resolve_cfg(artifactData);
cfg.captureConsoleOutput = false;
cfg.robust.enabled = true;

plnNom = resolve_pln(artifactData);
machine = resolve_machine(artifactData);
wEvalSource = '';
[wEval, wEvalSource] = resolve_warm_start(artifactData);

fprintf('Evaluation weight length: %d\n', numel(wEval));
fprintf('Evaluation weight source: %s\n', wEvalSource);
fprintf('Arc spacing: %.2f degrees\n', cfg.arcSpacing);
fprintf('Robust mode: %s\n\n', cfg.robust.mode);

startTime = tic;

plnRob = configure_robust_scenarios(plnNom, ct, cfg);
cstRob = configure_robust_objectives(cst, cfg);

fprintf('Resolving nominal-compatible STF...\n');
[stfRob, stfSource] = resolve_evaluation_stf(artifactData, ct, cst, plnNom, wEval);
fprintf('STF source: %s\n', stfSource);
fprintf('STF bixel count: %d\n\n', count_stf_bixels(stfRob));

timestamp = datestr(now, 'yyyy-mm-dd_HHMMSS');
[~, nominalBaseName] = fileparts(nominalArtifactFile);
prefix = sprintf('%s_robust_resume_%s', nominalBaseName, timestamp);

fprintf('Evaluating robust scenarios sequentially (memory-safe mode)...\n');
[robustEval, nominalMetrics, nominalDoseCube, worstCaseDoseCube] = ...
    evaluate_robust_scenarios_sequential(ct, cstRob, stfRob, plnRob, wEval, cfg, saveFullDij, outputDir, prefix);

meta = struct();
meta.methodName = 'Task 28-30 robust resume';
meta.sourceNominalArtifact = nominalArtifactFile;
meta.timestamp = datetime('now');
meta.runtimeSeconds = toc(startTime);
meta.outputPaths = struct('runDir', outputDir);
meta.robustMode = cfg.robust.mode;
meta.nScenarios = robustEval.nScenarios;
meta.saveFullDij = saveFullDij;
meta.executionMode = 'sequential_evaluation_only';
meta.runFinalReoptimization = false;
meta.warmStartSource = wEvalSource;
meta.stfSource = stfSource;

dvhOverlayPlot = save_robust_dvh_overlay(cstRob, nominalDoseCube, worstCaseDoseCube, outputDir, prefix);
summaryPlot = save_robust_summary_plot(robustEval, outputDir, prefix);

matFile = fullfile(outputDir, sprintf('%s.mat', prefix));
runInfo = struct();
runInfo.nominalArtifactFile = nominalArtifactFile;
runInfo.patientFile = patientFile;
runInfo.timestamp = timestamp;
runInfo.runtimeSeconds = meta.runtimeSeconds;
runInfo.saveFullDij = saveFullDij;
runInfo.warmStartLength = numel(wEval);
runInfo.warmStartSource = wEvalSource;
runInfo.robustMode = cfg.robust.mode;
runInfo.executionMode = meta.executionMode;
runInfo.stfSource = stfSource;

if saveFullDij
    save(matFile, ...
        'cfg', 'machine', 'plnRob', 'cstRob', 'stfRob', ...
        'wEval', 'robustEval', 'nominalMetrics', ...
        'nominalDoseCube', 'worstCaseDoseCube', ...
        'meta', 'runInfo', ...
        '-v7.3');
else
    save(matFile, ...
        'cfg', 'machine', 'plnRob', 'cstRob', 'stfRob', ...
        'wEval', 'robustEval', 'nominalMetrics', 'meta', 'runInfo', ...
        '-v7.3');
end

outputs = struct();
outputs.nominalArtifactFile = nominalArtifactFile;
outputs.matFile = matFile;
outputs.dvhPlot = dvhOverlayPlot;
outputs.summaryPlot = summaryPlot;
outputs.runtimeSeconds = meta.runtimeSeconds;
outputs.nScenarios = robustEval.nScenarios;
outputs.passRate = robustEval.passRate;
outputs.worstCase = robustEval.worstCase;
outputs.warmStartSource = wEvalSource;
outputs.stfSource = stfSource;

fprintf('\n========================================\n');
fprintf('RESUME COMPLETE\n');
fprintf('========================================\n');
fprintf('Saved MAT: %s\n', matFile);
fprintf('Saved DVH: %s\n', dvhOverlayPlot);
fprintf('Saved summary plot: %s\n', summaryPlot);
fprintf('Worst-case target D95: %.6f\n', robustEval.worstCase.targetD95Min);
fprintf('Worst-case OAR Dmax: %.6f\n', robustEval.worstCase.oarDmaxMax);
fprintf('Runtime: %.2f minutes\n\n', meta.runtimeSeconds / 60);
end

function nominalArtifactFile = resolve_latest_nominal_artifact(resultsDir)
listing = dir(fullfile(resultsDir, 'task25_27_pipeline_*.mat'));
if isempty(listing)
    error('run_task28_30_from_nominal:missingArtifact', ...
        'No task25_27_pipeline_*.mat artifact found in %s.', resultsDir);
end

[~, ix] = max([listing.datenum]);
nominalArtifactFile = fullfile(listing(ix).folder, listing(ix).name);
end

function patientFile = resolve_patient_file()
patientFile = which('PROSTATE.mat');
if ~isempty(patientFile)
    return;
end

matRadRcPath = which('matRad_rc');
if isempty(matRadRcPath)
    error('run_task28_30_from_nominal:missingPatient', ...
        'PROSTATE.mat was not found and matRad_rc is not on the MATLAB path.');
end

matRadRoot = fileparts(matRadRcPath);
candidate = fullfile(matRadRoot, 'data', 'PROSTATE.mat');
if exist(candidate, 'file') ~= 2
    error('run_task28_30_from_nominal:missingPatient', ...
        'Could not locate PROSTATE.mat at %s.', candidate);
end

patientFile = candidate;
end

function cfg = resolve_cfg(artifactData)
if isfield(artifactData, 'cfg') && isstruct(artifactData.cfg)
    cfg = artifactData.cfg;
else
    cfg = default_config();
end
end

function pln = resolve_pln(artifactData)
if ~isfield(artifactData, 'pln') || ~isstruct(artifactData.pln)
    error('run_task28_30_from_nominal:missingPlan', ...
        'The nominal artifact must contain a valid pln struct.');
end

pln = normalize_plan_machine_fields(artifactData.pln);
end

function machine = resolve_machine(artifactData)
if isfield(artifactData, 'machine') && isstruct(artifactData.machine)
    machine = artifactData.machine;
else
    machine = default_machine_model();
end
end

function [wInit, sourceName] = resolve_warm_start(artifactData)
candidates = {
    'resultFinal', ...
    'resultSpot', ...
    'resultLayer', ...
    'denseResult'
    };

for i = 1:numel(candidates)
    name = candidates{i};
    if ~isfield(artifactData, name) || ~isstruct(artifactData.(name))
        continue;
    end

    entry = artifactData.(name);
    if isfield(entry, 'w') && isnumeric(entry.w) && isvector(entry.w) && ~isempty(entry.w)
        wInit = reshape(double(entry.w), [], 1);
        sourceName = sprintf('%s.w', name);
        return;
    end
end

error('run_task28_30_from_nominal:missingWarmStart', ...
    ['The nominal artifact does not contain any usable warm start. ', ...
     'Expected one of resultFinal.w, resultSpot.w, resultLayer.w, or denseResult.w.']);
end

function [stf, sourceName] = resolve_evaluation_stf(artifactData, ct, cst, plnNom, wEval)
if isfield(artifactData, 'stf') && isstruct(artifactData.stf) && ~isempty(artifactData.stf)
    stf = artifactData.stf;
    sourceName = 'nominalArtifact.stf';
else
    fprintf('Nominal artifact has no saved STF; regenerating STF from nominal plan.\n');
    stf = matRad_generateStf(ct, cst, plnNom);
    sourceName = 'regenerated nominal plan STF';
end

nBixels = count_stf_bixels(stf);
if ~isfinite(nBixels) || nBixels <= 0
    error('run_task28_30_from_nominal:invalidStf', ...
        'Could not infer a positive bixel count from the nominal-compatible STF.');
end

if nBixels ~= numel(wEval)
    error('run_task28_30_from_nominal:stfWeightMismatch', ...
        ['Nominal-compatible STF has %d bixels but evaluation weight length is %d. ', ...
         'Use the exact nominal STF saved with the artifact or regenerate the nominal artifact.'], ...
        nBixels, numel(wEval));
end
end

function nBixels = count_stf_bixels(stf)
nBixels = 0;
if ~isstruct(stf) || isempty(stf)
    nBixels = NaN;
    return;
end

for b = 1:numel(stf)
    if ~isfield(stf(b), 'ray') || ~isstruct(stf(b).ray)
        continue;
    end
    for r = 1:numel(stf(b).ray)
        if isfield(stf(b).ray(r), 'energy') && ~isempty(stf(b).ray(r).energy)
            nBixels = nBixels + numel(stf(b).ray(r).energy);
        elseif isfield(stf(b).ray(r), 'rayPos_bev') && ~isempty(stf(b).ray(r).rayPos_bev)
            nBixels = nBixels + size(stf(b).ray(r).rayPos_bev, 1);
        end
    end
end
end

function plotFile = save_robust_summary_plot(robustEval, outputDir, prefix)
targetD95 = reshape([robustEval.scenarioMetrics.targetD95Min], [], 1);
oarDmax = reshape([robustEval.scenarioMetrics.oarDmaxMax], [], 1);

fig = figure('Visible', 'off', 'Color', 'w');
tiledlayout(2, 1);

nexttile;
bar(targetD95);
xlabel('Scenario index');
ylabel('Target D95');
title('Robust target D95 by scenario');
grid on;

nexttile;
plot(oarDmax, 'o-', 'LineWidth', 1.5);
xlabel('Scenario index');
ylabel('OAR Dmax');
title('Robust OAR Dmax by scenario');
grid on;

plotFile = fullfile(outputDir, sprintf('%s_summary.png', prefix));
set(fig, 'PaperPositionMode', 'auto');
print(fig, plotFile, '-dpng', '-r150');
close(fig);
end

function [robustEval, nominalMetrics, nominalDoseCube, worstCaseDoseCube] = evaluate_robust_scenarios_sequential(ct, cstRob, stfRob, plnRob, wEval, cfg, saveFullDij, outputDir, prefix)
scenarioCount = resolve_scenario_count(plnRob);
scenarioMetrics = repmat(struct( ...
    'scenarioIndex', 0, ...
    'doseMetrics', struct(), ...
    'targetD95Min', NaN, ...
    'targetD98Min', NaN, ...
    'targetV95Min', NaN, ...
    'oarDmaxMax', NaN, ...
    'scenarioInfo', struct()), scenarioCount, 1);

doseDims = [];
sumDose = [];
sumSqDose = [];
nominalDoseCube = [];
worstCaseDoseCube = [];
worstScenarioIx = NaN;
perScenarioFiles = cell(scenarioCount, 1);

for s = 1:scenarioCount
    fprintf('  Scenario %d of %d...\n', s, scenarioCount);
    plnSingle = plnRob;
    plnSingle.multScen = extract_single_scenario(plnRob.multScen, s);

    dijSingle = matRad_calcDoseInfluence(ct, cstRob, stfRob, plnSingle);
    if numel(wEval) ~= double(dijSingle.totalNumOfBixels)
        error('run_task28_30_from_nominal:dimensionMismatch', ...
            ['Evaluation weight length %d does not match scenario %d dij.totalNumOfBixels %d. ', ...
             'The nominal artifact may not correspond to the robust geometry.'], ...
            numel(wEval), s, double(dijSingle.totalNumOfBixels));
    end

    doseVector = compute_scenario_dose_vector(dijSingle, wEval);
    if isempty(doseDims)
        doseDims = infer_dose_dims_from_single(ct, dijSingle, doseVector);
        sumDose = zeros(doseDims);
        sumSqDose = zeros(doseDims);
    end

    doseCube = reshape(doseVector, doseDims);
    doseMetrics = compute_dvh_metrics(doseCube, ct, cstRob, cfg);

    scenarioMetrics(s).scenarioIndex = s;
    scenarioMetrics(s).doseMetrics = doseMetrics;
    scenarioMetrics(s).targetD95Min = doseMetrics.summary.targetD95Min;
    scenarioMetrics(s).targetD98Min = doseMetrics.summary.targetD98Min;
    scenarioMetrics(s).targetV95Min = doseMetrics.summary.targetV95Min;
    scenarioMetrics(s).oarDmaxMax = doseMetrics.summary.oarDmaxMax;
    scenarioMetrics(s).scenarioInfo = struct('scenarioIndex', s);

    if s == 1
        nominalDoseCube = doseCube;
        nominalMetrics = doseMetrics;
    end

    sumDose = sumDose + doseCube;
    sumSqDose = sumSqDose + doseCube .^ 2;

    if saveFullDij
        perScenarioFile = fullfile(outputDir, sprintf('%s_scen%02d.mat', prefix, s));
        save(perScenarioFile, 'doseCube', 'doseMetrics', 'plnSingle', '-v7.3');
        perScenarioFiles{s} = perScenarioFile;
    end
end

worstCase = compute_worst_case_metrics(scenarioMetrics);
if isfinite(worstCase.worstTargetD95Scenario)
    worstScenarioIx = worstCase.worstTargetD95Scenario;
else
    worstScenarioIx = 1;
end
worstCaseDoseCube = reconstruct_saved_or_rerun_scenario(ct, cstRob, stfRob, plnRob, wEval, worstScenarioIx, saveFullDij, perScenarioFiles);

scenarioPass = compute_scenario_pass_mask(scenarioMetrics, cfg);
meanDoseCube = sumDose / scenarioCount;
stdDoseCube = sqrt(max(sumSqDose / scenarioCount - meanDoseCube .^ 2, 0));

robustEval = struct();
robustEval.scenarioMetrics = scenarioMetrics;
robustEval.worstCase = worstCase;
robustEval.passRate = mean(double(scenarioPass));
robustEval.scenarioPass = scenarioPass;
robustEval.nominalScenarioIndex = 1;
robustEval.stdDoseCube = stdDoseCube;
robustEval.meanDoseCube = meanDoseCube;
robustEval.nScenarios = scenarioCount;
robustEval.mode = 'sequential_evaluation_only';
robustEval.perScenarioFiles = perScenarioFiles;
robustEval.worstCaseScenarioIndex = worstScenarioIx;
end

function scenarioCount = resolve_scenario_count(plnRob)
scenarioCount = 1;
if isstruct(plnRob) && isfield(plnRob, 'multScen')
    ms = plnRob.multScen;
    if isobject(ms) && isprop(ms, 'totNumScen') && ~isempty(ms.totNumScen)
        scenarioCount = double(ms.totNumScen);
    elseif isstruct(ms)
        % Plans loaded from MAT files can contain the scenario model as a
        % struct.  Use the explicit scenario table in that case; otherwise
        % silently evaluating only scenario 1 would invalidate the robust
        % check.
        if isfield(ms, 'scenForProb') && ~isempty(ms.scenForProb)
            scenarioCount = size(ms.scenForProb, 1);
        elseif isfield(ms, 'totNumScen') && isscalar(ms.totNumScen)
            scenarioCount = double(ms.totNumScen);
        end
    end
end
end

function singleScenario = extract_single_scenario(multScen, scenarioIndex)
if isobject(multScen) && ismethod(multScen, 'extractSingleScenario')
    singleScenario = multScen.extractSingleScenario(scenarioIndex);
else
    singleScenario = multScen;
end
end

function doseVector = compute_scenario_dose_vector(dijSingle, wEval)
physicalDose = dijSingle.physicalDose;
if iscell(physicalDose)
    physicalDose = physicalDose{1};
end

if isvector(physicalDose)
    doseVector = reshape(double(physicalDose), [], 1);
else
    doseVector = reshape(full(physicalDose * reshape(double(wEval), [], 1)), [], 1);
end
end

function doseDims = infer_dose_dims_from_single(ct, dijSingle, doseVector)
if isfield(dijSingle, 'doseGrid') && isfield(dijSingle.doseGrid, 'dimensions') && numel(dijSingle.doseGrid.dimensions) == 3
    doseDims = reshape(double(dijSingle.doseGrid.dimensions), 1, []);
elseif isstruct(ct) && isfield(ct, 'cubeDim') && numel(ct.cubeDim) == 3
    doseDims = reshape(double(ct.cubeDim), 1, []);
else
    doseDims = [numel(doseVector), 1, 1];
end
end

function scenarioPass = compute_scenario_pass_mask(scenarioMetrics, cfg)
scenarioPass = true(numel(scenarioMetrics), 1);
nominalTargetD95 = scenarioMetrics(1).targetD95Min;
toleranceFraction = get_robust_tolerance_fraction(cfg);

if ~(isfinite(nominalTargetD95) && nominalTargetD95 > 0)
    scenarioPass(:) = false;
    return;
end

for s = 1:numel(scenarioMetrics)
    currentD95 = scenarioMetrics(s).targetD95Min;
    if ~isfinite(currentD95)
        scenarioPass(s) = false;
        continue;
    end

    degradation = (nominalTargetD95 - currentD95) / nominalTargetD95;
    scenarioPass(s) = degradation <= toleranceFraction;
end
end

function toleranceFraction = get_robust_tolerance_fraction(cfg)
tolerancePercent = 2.0;
if isstruct(cfg) && isfield(cfg, 'evaluation') && isstruct(cfg.evaluation) && ...
        isfield(cfg.evaluation, 'robustD95TolerancePercent') && ...
        isnumeric(cfg.evaluation.robustD95TolerancePercent) && ...
        isscalar(cfg.evaluation.robustD95TolerancePercent) && ...
        isfinite(cfg.evaluation.robustD95TolerancePercent)
    tolerancePercent = double(cfg.evaluation.robustD95TolerancePercent);
end
toleranceFraction = tolerancePercent / 100;
end

function doseCube = reconstruct_saved_or_rerun_scenario(ct, cstRob, stfRob, plnRob, wEval, scenarioIndex, saveFullDij, perScenarioFiles)
if saveFullDij && scenarioIndex >= 1 && scenarioIndex <= numel(perScenarioFiles) && ~isempty(perScenarioFiles{scenarioIndex})
    data = load(perScenarioFiles{scenarioIndex}, 'doseCube');
    doseCube = data.doseCube;
    return;
end

plnSingle = plnRob;
plnSingle.multScen = extract_single_scenario(plnRob.multScen, scenarioIndex);
dijSingle = matRad_calcDoseInfluence(ct, cstRob, stfRob, plnSingle);
doseVector = compute_scenario_dose_vector(dijSingle, wEval);
doseDims = infer_dose_dims_from_single(ct, dijSingle, doseVector);
doseCube = reshape(doseVector, doseDims);
end

function plotFile = save_robust_dvh_overlay(cstRob, nominalDoseCube, worstCaseDoseCube, outputDir, prefix)
nominalDvh = matRad_calcDVH(cstRob, nominalDoseCube, 'cum');
worstDvh = matRad_calcDVH(cstRob, worstCaseDoseCube, 'cum');

fig = figure('Visible', 'off', 'Color', 'w');
hold on;

for i = 1:numel(nominalDvh)
    plot(nominalDvh(i).doseGrid, nominalDvh(i).volumePoints, '-', 'LineWidth', 1.5);
    plot(worstDvh(i).doseGrid, worstDvh(i).volumePoints, '--', 'LineWidth', 1.2);
end

legendLabels = cell(2 * size(cstRob, 1), 1);
for i = 1:size(cstRob, 1)
    baseName = stringify_label(cstRob{i, 2}, sprintf('Structure %d', i));
    legendLabels{2 * i - 1} = sprintf('%s nominal', baseName);
    legendLabels{2 * i} = sprintf('%s worst-case', baseName);
end

legend(legendLabels, 'Interpreter', 'none', 'Location', 'bestoutside');
xlabel('Dose');
ylabel('Volume (%)');
title('Worst-case versus nominal DVH overlay');
grid on;
box on;

plotFile = fullfile(outputDir, sprintf('%s_dvh_overlay.png', prefix));
set(fig, 'PaperPositionMode', 'auto');
print(fig, plotFile, '-dpng', '-r150');
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
