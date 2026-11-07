function summaryTable = run_support_quality_sweep_20deg(outputDir, arcSpacing, sourceDir)
% Revisit an arc-spacing spot history and compare locked-support final plans.
% Research-only: BDT is a proxy and quality tolerances are exploratory.

narginchk(0, 3);
workspaceDir = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(workspaceDir, 'config'));
addpath(fullfile(workspaceDir, 'planning'));
addpath(fullfile(workspaceDir, 'mapping'));
addpath(fullfile(workspaceDir, 'objective'));
addpath(fullfile(workspaceDir, 'delivery'));
addpath(fullfile(workspaceDir, 'reduction'));
addpath(fullfile(workspaceDir, 'robustness'));
addpath(fullfile(workspaceDir, 'evaluation'));

if nargin < 1 || isempty(outputDir)
    outputDir = fullfile(workspaceDir, 'results', 'support_quality_sweep_20deg');
end
if nargin < 2 || isempty(arcSpacing)
    arcSpacing = 20;
end
if nargin < 3 || isempty(sourceDir)
    if arcSpacing == 20
        sourceDir = fullfile(workspaceDir, 'results', 'support_fix_20deg', 'arc_20deg');
    else
        sourceDir = fullfile(workspaceDir, ...
            'results', 'full_nominal_then_robust_2026-06-02_003113', ...
            'nominal', 'task41', sprintf('arc_%gdeg', arcSpacing));
    end
end
if exist(outputDir, 'dir') ~= 7
    mkdir(outputDir);
end

artifactFiles = dir(fullfile(sourceDir, '*_nominal_artifact.mat'));
methodFiles = dir(fullfile(sourceDir, '*_methods.csv'));
if numel(artifactFiles) ~= 1 || numel(methodFiles) ~= 1
    error('run_support_quality_sweep_20deg:invalidSource', ...
        'Expected exactly one nominal artifact and one methods CSV in %s.', sourceDir);
end
source = load(fullfile(sourceDir, artifactFiles(1).name), ...
    'cfg', 'machine', 'denseResult', 'resultLayer', 'resultSpot');
reference = readtable(fullfile(sourceDir, methodFiles(1).name));
denseRow = reference(strcmp(string(reference.Method), "Dense Arc"), :);
if height(denseRow) ~= 1 || source.cfg.arcSpacing ~= arcSpacing
    error('run_support_quality_sweep_20deg:invalidSource', ...
        'The reference must contain one Dense Arc row from the requested spacing.');
end

cfg = source.cfg;
cfg.outputDir = outputDir;
cfg.captureConsoleOutput = false;
cfg.robust.enabled = false;
loadedPatient = load(cfg.patientFile, 'ct', 'cst');
ct = loadedPatient.ct;
cst = loadedPatient.cst;

fprintf('Rebuilding the %g-degree dose matrix...\n', arcSpacing);
[pln, stf, dij] = build_discrete_arc_plan(ct, cst, cfg);
if double(dij.totalNumOfBixels) ~= numel(source.denseResult.w) || ...
        numel(source.resultLayer.w) ~= double(dij.totalNumOfBixels)
    error('run_support_quality_sweep_20deg:geometryMismatch', ...
        'Rebuilt matrix and saved reference weights have different sizes.');
end

if iscell(dij.physicalDose)
    D = dij.physicalDose{1};
else
    D = dij.physicalDose;
end
denseDose = compute_dvh_metrics(full(D * double(source.denseResult.w(:))), ct, cst, cfg);
check_metric('Dense D95', denseDose.summary.targetD95Min, ...
    double(denseRow.dose_targetD95Min), 0.05);
check_metric('Dense D98', denseDose.summary.targetD98Min, ...
    double(denseRow.dose_targetD98Min), 0.05);
check_metric('Dense OAR Dmax', denseDose.summary.oarDmaxMax, ...
    double(denseRow.dose_oarDmaxMax), 0.05);
spotDose = compute_dvh_metrics(full(D * double(source.resultSpot.w(:))), ct, cst, cfg);
check_metric('Saved sparse D95', spotDose.summary.targetD95Min, ...
    double(source.resultSpot.history(end).ctvD95), 0.05);
fprintf('Dose matrix matched the September reference metrics.\n');

map = map_bixel_to_cp_energy_spot(dij, stf, pln);
[groups, ~] = build_group_index(map);
targetFractions = [0.20, 0.10, 0.05, 0.03, 0.015];
targets = max(1, round(targetFractions * double(dij.totalNumOfBixels)));
existingFiles = dir(fullfile(outputDir, 'support_*.mat'));
if numel(existingFiles) == numel(targets)
    fprintf('Resuming from %d saved support snapshots.\n', numel(existingFiles));
    snapshots = repmat(struct('iteration', 0, 'nActiveSpots', 0, ...
        'w', [], 'supportMask', []), numel(existingFiles), 1);
    for i = 1:numel(existingFiles)
        saved = load(fullfile(outputDir, existingFiles(i).name), 'snap');
        snapshots(i) = saved.snap;
        if numel(snapshots(i).w) ~= double(dij.totalNumOfBixels) || ...
                numel(snapshots(i).supportMask) ~= double(dij.totalNumOfBixels) || ...
                nnz(snapshots(i).supportMask) ~= snapshots(i).nActiveSpots || ...
                any(snapshots(i).w(~snapshots(i).supportMask) ~= 0)
            error('run_support_quality_sweep_20deg:invalidSnapshot', ...
                'Saved support snapshot %s is inconsistent.', existingFiles(i).name);
        end
    end
    counts = [snapshots.nActiveSpots];
    selected = zeros(numel(targets), 1);
    for i = 1:numel(targets)
        [~, selected(i)] = min(abs(counts - targets(i)));
    end
else
    fprintf('Replaying spot stage while capturing accepted support states...\n');
    [spotResult, ~] = run_spot_sparse_stage(dij, cst, pln, ...
        source.resultLayer.w, map, groups, cfg, source.machine, true);
    snapshots = spotResult.supportSnapshots;
    if isempty(snapshots)
        error('run_support_quality_sweep_20deg:noSnapshots', ...
            'The spot stage produced no accepted support states.');
    end
    selected = zeros(numel(targets), 1);
    counts = [snapshots.nActiveSpots];
    for i = 1:numel(targets)
        [~, selected(i)] = min(abs(counts - targets(i)));
    end
    historyTable = struct2table(spotResult.history);
    writetable(historyTable, fullfile(outputDir, 'spot_history.csv'));
end
if numel(unique(selected)) ~= numel(targets)
    error('run_support_quality_sweep_20deg:duplicateSupport', ...
        'Nearest history points are not distinct; inspect the support files.');
end

records = repmat(struct( ...
    'targetSpots', NaN, 'historyIteration', NaN, 'supportSpots', NaN, ...
    'finalActiveSpots', NaN, 'finalActiveLayers', NaN, ...
    'bdtSeconds', NaN, 'bdtReductionPercent', NaN, ...
    'D95', NaN, 'D98', NaN, 'oarDmaxMax', NaN, ...
    'D95ChangePercent', NaN, 'D98ChangePercent', NaN, ...
    'oarDmaxChangePercent', NaN, ...
    'matchedWithin1Percent', false, 'matchedWithin2Percent', false, ...
    'solverSeconds', NaN, 'status', "pending", 'errorMessage', ""), ...
    numel(targets), 1);

for i = 1:numel(targets)
    snap = snapshots(selected(i));
    records(i).targetSpots = targets(i);
    records(i).historyIteration = snap.iteration;
    records(i).supportSpots = snap.nActiveSpots;
    save(fullfile(outputDir, sprintf('support_%05d.mat', snap.nActiveSpots)), ...
        'snap', '-v7.3');
    fprintf('Final locked-support re-opt: target %d, actual %d spots...\n', ...
        targets(i), snap.nActiveSpots);
    try
        startTime = tic;
        final = run_final_reoptimization(dij, cst, pln, snap.w, ...
            cfg, snap.supportMask);
        records(i).solverSeconds = toc(startTime);
        if any(abs(final.w(~snap.supportMask)) > 1e-10)
            error('run_support_quality_sweep_20deg:supportLeakage', ...
                'A final plan activated weights outside its locked support.');
        end
        dose = compute_dvh_metrics(full(D * double(final.w(:))), ct, cst, cfg);
        delivery = compute_delivery_metrics(final.w, map, groups, cfg, source.machine);

        records(i).finalActiveSpots = delivery.nActiveSpots;
        records(i).finalActiveLayers = delivery.nActiveLayers;
        records(i).bdtSeconds = delivery.T_total;
        records(i).bdtReductionPercent = 100 * ...
            (double(denseRow.delivery_T_total) - delivery.T_total) / ...
            double(denseRow.delivery_T_total);
        records(i).D95 = dose.summary.targetD95Min;
        records(i).D98 = dose.summary.targetD98Min;
        records(i).oarDmaxMax = dose.summary.oarDmaxMax;
        records(i).D95ChangePercent = percent_change( ...
            records(i).D95, denseDose.summary.targetD95Min);
        records(i).D98ChangePercent = percent_change( ...
            records(i).D98, denseDose.summary.targetD98Min);
        records(i).oarDmaxChangePercent = percent_change( ...
            records(i).oarDmaxMax, denseDose.summary.oarDmaxMax);
        records(i).matchedWithin1Percent = quality_pass(records(i), 1);
        records(i).matchedWithin2Percent = quality_pass(records(i), 2);
        records(i).status = "completed";

        wFinal = final.w; %#ok<NASGU>
        supportMask = snap.supportMask; %#ok<NASGU>
        doseSummary = dose.summary; %#ok<NASGU>
        deliveryMetrics = delivery; %#ok<NASGU>
        save(fullfile(outputDir, sprintf('final_%05d.mat', snap.nActiveSpots)), ...
            'wFinal', 'supportMask', 'doseSummary', 'deliveryMetrics', '-v7.3');
    catch ME
        records(i).status = "failed";
        records(i).errorMessage = string(ME.message);
        fprintf(2, 'Support %d failed: %s\n', snap.nActiveSpots, ME.message);
    end
    summaryTable = struct2table(records);
    writetable(summaryTable, fullfile(outputDir, 'support_quality_summary.csv'));
end

summaryTable = struct2table(records);
if any(strcmp(summaryTable.status, "failed"))
    error('run_support_quality_sweep_20deg:incomplete', ...
        'One or more final support optimizations failed; see support_quality_summary.csv.');
end
eligible = strcmp(summaryTable.status, "completed") & ...
    summaryTable.matchedWithin2Percent;
if any(eligible)
    indices = find(eligible);
    [~, localBest] = max(summaryTable.bdtReductionPercent(indices));
    winner = summaryTable(indices(localBest), :);
    writetable(winner, fullfile(outputDir, 'selected_matched_quality.csv'));
    fprintf('Selected support: %d spots, BDT reduction %.2f%%.\n', ...
        winner.supportSpots, winner.bdtReductionPercent);
else
    fprintf('No support met the exploratory 2%% dose-quality rule.\n');
end
end

function check_metric(name, actual, expected, tolerance)
if ~isfinite(actual) || abs(actual - expected) > tolerance
    error('run_support_quality_sweep_20deg:doseMismatch', ...
        '%s: rebuilt %.6f differs from reference %.6f.', name, actual, expected);
end
end

function change = percent_change(value, baseline)
change = 100 * (double(value) - double(baseline)) / double(baseline);
end

function pass = quality_pass(record, tolerancePercent)
% Exploratory rule: target D95/D98 may fall by at most the tolerance;
% maximum OAR Dmax may rise by at most the same percentage.
pass = isfinite(record.D95ChangePercent) && ...
    isfinite(record.D98ChangePercent) && ...
    isfinite(record.oarDmaxChangePercent) && ...
    record.D95ChangePercent >= -tolerancePercent && ...
    record.D98ChangePercent >= -tolerancePercent && ...
    record.oarDmaxChangePercent <= tolerancePercent;
end
