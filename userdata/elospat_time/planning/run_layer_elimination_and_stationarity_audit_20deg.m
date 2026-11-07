function outputs = run_layer_elimination_and_stationarity_audit_20deg(outputDir, sweepDir)
% Audit eliminated layers and force controlled multi-iteration validation.
%
% The group-half proximal implementation is not modified. This runner only
% disables the two relative stopping tolerances for controlled experiments.

narginchk(0, 2);
root = fileparts(fileparts(mfilename('fullpath')));
add_module_paths(root);
if nargin < 1 || isempty(outputDir)
    outputDir = fullfile(root, 'results', ...
        ['layer_stationarity_audit_20deg_' datestr(now, 'yyyymmdd_HHMMSS')]);
end
if nargin < 2 || isempty(sweepDir)
    sweepDir = fullfile(root, 'results', ...
        'representative_layer_sweep_20deg_job21878');
end
create_empty_output_directory(outputDir);

artifactFile = fullfile(root, 'results', 'task25_27_pipeline_20deg.mat');
baselineFile = fullfile(root, 'results', ...
    'support_representativeness_20deg_job21861', ...
    'support_86400_checkpoint.mat');
selectedFile = fullfile(sweepDir, 'lambda_alpha_0p01_checkpoint.mat');
if exist(baselineFile, 'file') ~= 2 || exist(selectedFile, 'file') ~= 2
    error('run_layer_elimination_and_stationarity_audit_20deg:missingInput', ...
        'Required zero-lambda or alpha=0.01 checkpoint is missing.');
end

fprintf('Loading existing dij and saved checkpoints; no dose calculation rebuild.\n');
source = load(artifactFile, 'cfg', 'machine', 'denseResult', ...
    'dij', 'pln', 'map', 'groups');
base = load(baselineFile, 'supportMask', 'final', 'finalDose', ...
    'finalDelivery', 'record');
selected = load(selectedFile, 'optimizer', 'layerStats', 'cfgTrial', ...
    'record', 'supportFull', 'finalDose', 'finalDelivery');
patient = load(source.cfg.patientFile, 'ct', 'cst');

cfg = merge_structs_local(default_config(), source.cfg);
cfg.mode = 'integrated'; cfg.robust.enabled = false;
cfg.integrated.componentStage = 'layer';
cfg.lambdaB = 0; cfg.lambdaS = 0; cfg.lambdaES = 0; cfg.lambdaT = 0;
cfg.integrated.constraintPenalty = 1e3;
cfg.integrated.runGradientCheck = false;
cfg.integrated.verbose = true;
cfg.captureConsoleOutput = false;
cfg.enforceSelectedMinimumWeight = cfg.integrated.enforceFinalMinimumWeight;

dij = source.dij; pln = source.pln; map = source.map;
groupsFull = source.groups; machine = source.machine;
denseW = reshape(double(source.denseResult.w), [], 1);
nFull = double(dij.totalNumOfBixels);
validationMask = reshape(logical(base.supportMask), [], 1);
[dijReduced, reducedToFull, ~, mapReduced] = ...
    reduce_dij_to_support(dij, validationMask, map);
for i = 1:numel(mapReduced)
    mapReduced(i).originalGlobalIx = mapReduced(i).globalIx;
    mapReduced(i).globalIx = i;
end
[groupsReduced, cpGroupsReduced] = build_group_index(mapReduced);
doseModel = build_matrad_dose_model( ...
    dijReduced, patient.cst, pln, cfg.integrated.constraintPenalty);
context = build_integrated_context( ...
    groupsReduced, cpGroupsReduced, mapReduced, machine, cfg);
wZeroFull = reshape(double(base.final.w), [], 1);
w0 = wZeroFull(reducedToFull);
D = nominal_matrix(dij);

%% 1. Four-layer audit.
zeroGroups = find(selected.layerStats.exactZeroMask);
if numel(zeroGroups) ~= 4
    error('run_layer_elimination_and_stationarity_audit_20deg:unexpectedZeros', ...
        'Expected four alpha=0.01 exact-zero groups, found %d.', numel(zeroGroups));
end
denseDelivery = compute_delivery_metrics(denseW, map, groupsFull, cfg, machine);
zeroDelivery = compute_delivery_metrics(wZeroFull, map, groupsFull, cfg, machine);
layerRows = repmat(empty_layer_row(), numel(zeroGroups), 1);
structureRows = repmat(empty_structure_row(), 0, 1);
for k = 1:numel(zeroGroups)
    gReduced = zeroGroups(k);
    reducedGroup = groupsReduced(gReduced);
    gFull = match_full_group(reducedGroup, groupsFull);
    fullIx = reshape(double(groupsFull(gFull).ix), [], 1);
    reducedIx = reshape(double(reducedGroup.ix), [], 1);
    retainedFullIx = reducedToFull(reducedIx);
    if numel(retainedFullIx) ~= 1
        warning('Layer %d has %d retained spots, not one.', ...
            gReduced, numel(retainedFullIx));
    end

    denseLayerW = denseW(fullIx);
    reducedInitialW = w0(reducedIx);
    denseContribution = D(:, fullIx) * denseLayerW;
    reducedContribution = D(:, retainedFullIx) * reducedInitialW;
    denseMarginal = marginal_bdt(denseW, fullIx, map, groupsFull, cfg, machine);
    reducedMarginal = marginal_bdt( ...
        wZeroFull, retainedFullIx, map, groupsFull, cfg, machine);

    row = empty_layer_row();
    row.reducedGroupIndex = gReduced;
    row.fullGroupIndex = gFull;
    row.cpIndex = double(reducedGroup.cpIndex);
    row.energyIndex = double(reducedGroup.energyIndex);
    row.energyMeV = double(reducedGroup.energy);
    row.reducedSpotCount = numel(reducedIx);
    row.denseSpotCount = numel(fullIx);
    row.denseActiveSpotCount = nnz(denseLayerW > cfg.spotThreshold);
    row.retainedFullSpotIndex = retainedFullIx(1);
    row.reducedInitialWeight = reducedInitialW(1);
    row.reducedInitialWeightSum = sum(reducedInitialW);
    row.reducedInitialLayerNorm = norm(reducedInitialW, 2);
    row.denseWeightSum = sum(denseLayerW);
    row.denseLayerNorm = norm(denseLayerW, 2);
    row.denseWeightMin = min(denseLayerW);
    row.denseWeightMax = max(denseLayerW);
    row.alpha001FinalLayerNorm = selected.layerStats.groupNorm(gReduced);
    row.reducedDoseL2 = norm(reducedContribution, 2);
    row.reducedDoseMax = max(reducedContribution);
    row.denseDoseL2 = norm(denseContribution, 2);
    row.denseDoseMax = max(denseContribution);
    row.reducedMarginalBDT = reducedMarginal.T_total;
    row.reducedMarginalEnergyTime = reducedMarginal.T_energy;
    row.reducedMarginalSpotTime = reducedMarginal.T_spot;
    row.reducedMarginalIrradiationTime = reducedMarginal.T_irradiation;
    row.denseMarginalBDT = denseMarginal.T_total;
    row.denseMarginalEnergyTime = denseMarginal.T_energy;
    row.denseMarginalSpotTime = denseMarginal.T_spot;
    row.denseMarginalIrradiationTime = denseMarginal.T_irradiation;
    layerRows(k) = row;

    structureRows = [structureRows; contribution_by_structure( ...
        k, row, denseContribution, reducedContribution, patient.cst)]; %#ok<AGROW>
end
layerTable = struct2table(layerRows, 'AsArray', true);
structureTable = struct2table(structureRows, 'AsArray', true);
writetable(layerTable, fullfile(outputDir, 'eliminated_layer_audit.csv'));
writetable(structureTable, fullfile(outputDir, ...
    'eliminated_layer_structure_dose_contributions.csv'));

%% 2. Audit original one-update stopping and stationarity.
originalCfg = selected.cfgTrial;
[originalInitialObjective, ~] = evaluate_integrated_objective( ...
    w0, doseModel, context, originalCfg);
originalStationarity = stationarity_report( ...
    selected.optimizer.w, doseModel, context, originalCfg, ...
    selected.optimizer.history.stepSize(end));
originalStopping = struct( ...
    'alpha', selected.record.alpha, ...
    'lambdaEL', selected.record.lambdaEL, ...
    'iterations', selected.optimizer.iterations, ...
    'status', selected.optimizer.status, ...
    'stepSize', selected.optimizer.history.stepSize(end), ...
    'relativeStep', selected.optimizer.history.relativeStep(end), ...
    'relativeStepTolerance', originalCfg.integrated.relativeStepTolerance, ...
    'relativeObjective', selected.optimizer.history.relativeObjective(end), ...
    'relativeObjectiveTolerance', originalCfg.integrated.relativeObjectiveTolerance, ...
    'objectiveDecrease', originalInitialObjective - selected.optimizer.objective, ...
    'bothStoppingTestsPass', ...
        selected.optimizer.history.relativeStep(end) <= ...
            originalCfg.integrated.relativeStepTolerance && ...
        selected.optimizer.history.relativeObjective(end) <= ...
            originalCfg.integrated.relativeObjectiveTolerance, ...
    'proxDisplacementNorm', originalStationarity.proxDisplacementNorm, ...
    'relativeProxDisplacement', originalStationarity.relativeProxDisplacement, ...
    'proxGradientResidualNorm', originalStationarity.proxGradientResidualNorm, ...
    'relativeProxGradientResidual', ...
        originalStationarity.relativeProxGradientResidual);
writetable(struct2table(originalStopping), ...
    fullfile(outputDir, 'original_one_iteration_stationarity.csv'));

%% Controlled runs: same mathematics, stopping tolerances forced to zero.
controlledAlpha = [0.01; 1];
controlledRecords = repmat(empty_controlled_record(), numel(controlledAlpha), 1);
scale = selected.record.lambdaEL / selected.record.alpha;
for iTrial = 1:numel(controlledAlpha)
    alpha = controlledAlpha(iTrial);
    cfgTrial = cfg;
    cfgTrial.lambdaEL = alpha * scale;
    cfgTrial.integrated.maxIterations = 40;
    cfgTrial.integrated.relativeStepTolerance = 0;
    cfgTrial.integrated.relativeObjectiveTolerance = 0;
    token = alpha_token(alpha);
    fprintf('\n=== Forced 40-iteration audit alpha=%g ===\n', alpha);
    try
        [initialObjective, initialParts] = evaluate_integrated_objective( ...
            w0, doseModel, context, cfgTrial);
        optimizer = run_integrated_optimizer( ...
            doseModel, context, w0, configure_integrated_stage(cfgTrial));
        stats = exact_layer_stats_local(optimizer.w, groupsReduced);
        extracted = extract_integrated_support(optimizer.w, w0, context, cfgTrial);
        supportFull = false(nFull, 1);
        supportFull(reducedToFull(extracted.mask)) = true;
        wSupportFull = zeros(nFull, 1);
        wSupportFull(reducedToFull) = extracted.w;
        final = run_final_reoptimization(dij, patient.cst, pln, ...
            wSupportFull, cfgTrial, supportFull, extracted.layerMask);
        finalDose = dose_metrics(D, final.w, patient, cfgTrial);
        finalDelivery = compute_delivery_metrics( ...
            final.w, map, groupsFull, cfgTrial, machine);
        finalSolver = compute_solver_metrics( ...
            final.resultGUI, struct(), final.optimizer);
        supportPreserved = ~any(abs(final.w(~supportFull)) > 1e-10);
        stationarity = stationarity_report(optimizer.w, doseModel, context, ...
            cfgTrial, optimizer.history.stepSize(end));
        trajectory = controlled_trajectory(initialObjective, initialParts, optimizer);
        [dosePass, doseTable, worstAdverse] = dose_quality_detail( ...
            finalDose, dose_metrics(D, denseW, patient, cfgTrial), ...
            cfg.evaluation.ctvD95TolerancePercent);
        controlledRecords(iTrial) = make_controlled_record(alpha, cfgTrial, ...
            optimizer, stats, extracted, finalDelivery, zeroDelivery, ...
            finalDose, finalSolver, supportPreserved, stationarity, ...
            initialObjective, ...
            dosePass, worstAdverse);
        record = controlledRecords(iTrial); %#ok<NASGU>
        writetable(trajectory, fullfile(outputDir, ...
            sprintf('controlled_alpha_%s_trajectory.csv', token)));
        writetable(layer_group_table(stats), fullfile(outputDir, ...
            sprintf('controlled_alpha_%s_layer_groups.csv', token)));
        writetable(doseTable, fullfile(outputDir, ...
            sprintf('controlled_alpha_%s_final_ipopt_dose.csv', token)));
        save(fullfile(outputDir, sprintf( ...
            'controlled_alpha_%s_checkpoint.mat', token)), ...
            'record', 'optimizer', 'stats', 'extracted', 'supportFull', ...
            'final', 'finalDose', 'finalDelivery', 'finalSolver', ...
            'stationarity', 'trajectory', 'doseTable', 'cfgTrial', '-v7.3');
    catch trialError
        controlledRecords(iTrial) = failed_controlled_record( ...
            alpha, alpha * scale, trialError);
        record = controlledRecords(iTrial); %#ok<NASGU>
        errorInfo = struct('identifier', trialError.identifier, ...
            'message', trialError.message); %#ok<NASGU>
        save(fullfile(outputDir, sprintf( ...
            'controlled_alpha_%s_checkpoint.mat', token)), ...
            'record', 'errorInfo', 'cfgTrial', '-v7.3');
        warning('Controlled alpha=%g failed: %s. Continuing.', ...
            alpha, trialError.message);
    end
    writetable(struct2table(controlledRecords(1:iTrial), 'AsArray', true), ...
        fullfile(outputDir, 'controlled_multi_iteration_summary.csv'));
end

save(fullfile(outputDir, 'layer_stationarity_audit_complete.mat'), ...
    'layerTable', 'structureTable', 'originalStopping', ...
    'originalStationarity', 'controlledRecords', 'denseDelivery', ...
    'zeroDelivery', 'artifactFile', 'baselineFile', 'selectedFile', '-v7.3');
outputs = struct('outputDir', outputDir, 'layerTable', layerTable, ...
    'originalStopping', originalStopping, ...
    'controlledRecords', controlledRecords);
fprintf('\nLAYER ELIMINATION AND STATIONARITY AUDIT COMPLETE.\n');
end

function add_module_paths(root)
folders = {'config','planning','mapping','objective','delivery','reduction', ...
    'robustness','evaluation','integrated'};
for i = 1:numel(folders), addpath(fullfile(root, folders{i})); end
end

function create_empty_output_directory(outputDir)
if exist(outputDir, 'dir') == 7
    entries = dir(outputDir);
    entries = entries(~ismember({entries.name}, {'.', '..'}));
    if ~isempty(entries), error('Output directory must be empty: %s', outputDir); end
else
    mkdir(outputDir);
end
end

function D = nominal_matrix(dij)
if iscell(dij.physicalDose), D = dij.physicalDose{1}; else, D = dij.physicalDose; end
end

function g = match_full_group(reducedGroup, fullGroups)
cp = double(reducedGroup.cpIndex); ei = double(reducedGroup.energyIndex);
mask = reshape(double([fullGroups.cpIndex]), [], 1) == cp & ...
    reshape(double([fullGroups.energyIndex]), [], 1) == ei;
g = find(mask);
if numel(g) ~= 1
    error('Could not uniquely match CP %d energy index %d.', cp, ei);
end
end

function delta = marginal_bdt(w, ix, map, groups, cfg, machine)
withLayer = compute_delivery_metrics(w, map, groups, cfg, machine);
withoutW = w; withoutW(ix) = 0;
withoutLayer = compute_delivery_metrics(withoutW, map, groups, cfg, machine);
delta = struct('T_total', withLayer.T_total - withoutLayer.T_total, ...
    'T_energy', withLayer.T_energy - withoutLayer.T_energy, ...
    'T_spot', withLayer.T_spot - withoutLayer.T_spot, ...
    'T_irradiation', withLayer.T_irradiation - withoutLayer.T_irradiation, ...
    'T_gantry', withLayer.T_gantry - withoutLayer.T_gantry, ...
    'T_buffer', withLayer.T_buffer - withoutLayer.T_buffer);
end

function rows = contribution_by_structure(layerNumber, layerRow, denseDose, ...
    reducedDose, cst)
n = size(cst, 1);
rows = repmat(empty_structure_row(), n, 1);
for i = 1:n
    voxels = cst{i,4};
    if iscell(voxels), voxels = voxels{1}; end
    voxels = unique(double(voxels(:)));
    voxels = voxels(voxels >= 1 & voxels <= numel(denseDose));
    row = empty_structure_row();
    row.layerNumber = layerNumber;
    row.reducedGroupIndex = layerRow.reducedGroupIndex;
    row.cpIndex = layerRow.cpIndex;
    row.energyMeV = layerRow.energyMeV;
    row.structure = cst{i,2};
    if ~isempty(voxels)
        row.reducedMeanDoseContribution = mean(reducedDose(voxels));
        row.reducedMaxDoseContribution = max(reducedDose(voxels));
        row.denseMeanDoseContribution = mean(denseDose(voxels));
        row.denseMaxDoseContribution = max(denseDose(voxels));
    end
    rows(i) = row;
end
end

function s = stationarity_report(w, doseModel, context, cfg, step)
[~, gradient] = evaluate_integrated_smooth(w, doseModel, context, cfg, true);
proxPoint = prox_rms_elospat(w - step * gradient, step, ...
    context.groups, cfg.lambdaEL, zeros(size(w)));
displacement = w - proxPoint;
s = struct('stepSize', step, 'weightNorm', norm(w), ...
    'gradientNorm', norm(gradient), ...
    'proxDisplacementNorm', norm(displacement), ...
    'relativeProxDisplacement', norm(displacement) / max(1, norm(w)), ...
    'proxGradientResidualNorm', norm(displacement) / step, ...
    'relativeProxGradientResidual', ...
        norm(displacement) / (step * max(1, norm(w))));
end

function m = dose_metrics(D, w, patient, cfg)
m = compute_dvh_metrics(full(D * reshape(double(w), [], 1)), ...
    patient.ct, patient.cst, cfg);
end

function stats = exact_layer_stats_local(w, groups)
n = numel(groups); norms = zeros(n, 1); zeroMask = false(n, 1);
for g = 1:n
    values = w(groups(g).ix); norms(g) = norm(values, 2);
    zeroMask(g) = all(values == 0);
end
stats = struct('nGroups', n, 'nExactZeroLayers', nnz(zeroMask), ...
    'nExactActiveLayers', n - nnz(zeroMask), 'exactZeroMask', zeroMask, ...
    'groupNorm', norms, 'groupIndex', (1:n)', ...
    'cpIndex', reshape(double([groups.cpIndex]), [], 1), ...
    'energyIndex', reshape(double([groups.energyIndex]), [], 1), ...
    'energy', reshape(double([groups.energy]), [], 1));
end

function tbl = controlled_trajectory(initialObjective, initialParts, optimizer)
n = optimizer.iterations + 1; iteration = (0:n-1)';
objective = [initialObjective; optimizer.history.objective(:)];
doseObjective = [initialParts.dose; optimizer.history.doseObjective(:)];
layerPenalty = [initialParts.layerL21half; optimizer.history.layerPenalty(:)];
maxConstraintViolation = [initialParts.maxConstraintViolation; ...
    optimizer.history.maxConstraintViolation(:)];
relativeStep = [NaN; optimizer.history.relativeStep(:)];
relativeObjective = [NaN; optimizer.history.relativeObjective(:)];
stepSize = [NaN; optimizer.history.stepSize(:)];
tbl = table(iteration, objective, doseObjective, layerPenalty, ...
    maxConstraintViolation, relativeStep, relativeObjective, stepSize);
end

function tbl = layer_group_table(stats)
tbl = table(stats.groupIndex, stats.cpIndex, stats.energyIndex, stats.energy, ...
    stats.groupNorm, stats.exactZeroMask, 'VariableNames', ...
    {'groupIndex','cpIndex','energyIndex','energy','groupNorm','exactZero'});
end

function [pass, tbl, worstAdverse] = dose_quality_detail(candidate, dense, tolerance)
n = numel(dense.structures);
rows = repmat(struct('structure', '', 'type', '', 'isTarget', false, ...
    'denseD98', NaN, 'finalD98', NaN, 'D98ChangePercent', NaN, ...
    'denseD95', NaN, 'finalD95', NaN, 'D95ChangePercent', NaN, ...
    'denseDmean', NaN, 'finalDmean', NaN, 'DmeanChangePercent', NaN, ...
    'denseDmax', NaN, 'finalDmax', NaN, 'DmaxChangePercent', NaN, ...
    'adverseChangePercent', NaN, 'structurePass', true), n, 1);
worstAdverse = 0;
for i = 1:n
    d = dense.structures(i); f = candidate.structures(i);
    rows(i).structure = d.name; rows(i).type = d.type; rows(i).isTarget = d.isTarget;
    rows(i).denseD98 = d.D98; rows(i).finalD98 = f.D98;
    rows(i).denseD95 = d.D95; rows(i).finalD95 = f.D95;
    rows(i).denseDmean = d.Dmean; rows(i).finalDmean = f.Dmean;
    rows(i).denseDmax = d.Dmax; rows(i).finalDmax = f.Dmax;
    rows(i).D98ChangePercent = pct(f.D98, d.D98);
    rows(i).D95ChangePercent = pct(f.D95, d.D95);
    rows(i).DmeanChangePercent = pct(f.Dmean, d.Dmean);
    rows(i).DmaxChangePercent = pct(f.Dmax, d.Dmax);
    if d.isTarget
        adverse = max([0, -rows(i).D95ChangePercent, -rows(i).D98ChangePercent]);
    elseif strcmpi(d.type, 'OAR')
        adverse = max([0, rows(i).DmeanChangePercent, rows(i).DmaxChangePercent]);
    else
        adverse = 0;
    end
    rows(i).adverseChangePercent = adverse;
    rows(i).structurePass = adverse <= tolerance;
    worstAdverse = max(worstAdverse, adverse);
end
tbl = struct2table(rows); pass = all([rows.structurePass]);
end

function r = empty_layer_row()
r = struct('reducedGroupIndex', NaN, 'fullGroupIndex', NaN, ...
    'cpIndex', NaN, 'energyIndex', NaN, 'energyMeV', NaN, ...
    'reducedSpotCount', NaN, 'denseSpotCount', NaN, ...
    'denseActiveSpotCount', NaN, 'retainedFullSpotIndex', NaN, ...
    'reducedInitialWeight', NaN, 'reducedInitialWeightSum', NaN, ...
    'reducedInitialLayerNorm', NaN, 'denseWeightSum', NaN, ...
    'denseLayerNorm', NaN, 'denseWeightMin', NaN, 'denseWeightMax', NaN, ...
    'alpha001FinalLayerNorm', NaN, 'reducedDoseL2', NaN, ...
    'reducedDoseMax', NaN, 'denseDoseL2', NaN, 'denseDoseMax', NaN, ...
    'reducedMarginalBDT', NaN, 'reducedMarginalEnergyTime', NaN, ...
    'reducedMarginalSpotTime', NaN, 'reducedMarginalIrradiationTime', NaN, ...
    'denseMarginalBDT', NaN, 'denseMarginalEnergyTime', NaN, ...
    'denseMarginalSpotTime', NaN, 'denseMarginalIrradiationTime', NaN);
end

function r = empty_structure_row()
r = struct('layerNumber', NaN, 'reducedGroupIndex', NaN, 'cpIndex', NaN, ...
    'energyMeV', NaN, 'structure', '', ...
    'reducedMeanDoseContribution', NaN, 'reducedMaxDoseContribution', NaN, ...
    'denseMeanDoseContribution', NaN, 'denseMaxDoseContribution', NaN);
end

function r = empty_controlled_record()
r = struct('alpha', NaN, 'lambdaEL', NaN, 'iterations', NaN, ...
    'integratedStatus', '', 'objectiveInitial', NaN, 'objectiveFinal', NaN, ...
    'objectiveDecrease', NaN, 'finalRelativeStep', NaN, ...
    'finalRelativeObjective', NaN, 'finalStepSize', NaN, ...
    'proxGradientResidualNorm', NaN, 'relativeProxGradientResidual', NaN, ...
    'exactZeroLayers', NaN, 'finalActiveLayers', NaN, ...
    'finalActiveSpots', NaN, 'estimatedBDT', NaN, ...
    'bdtChangeFromZeroPercent', NaN, 'targetD95', NaN, 'targetD98', NaN, ...
    'worstAdverseDoseChangePercent', NaN, 'doseQualityPass', false, ...
    'maxConstraintViolation', NaN, 'supportPreserved', false, ...
    'finalIpoptStatus', '', 'runStatus', 'pending', ...
    'errorIdentifier', '', 'errorMessage', '');
end

function r = make_controlled_record(alpha, cfg, optimizer, stats, extracted, ...
    delivery, zeroDelivery, dose, finalSolver, supportPreserved, stationarity, ...
    initialObjective, dosePass, worstAdverse)
r = empty_controlled_record();
r.alpha = alpha; r.lambdaEL = cfg.lambdaEL;
r.iterations = optimizer.iterations; r.integratedStatus = optimizer.status;
r.objectiveInitial = initialObjective;
r.objectiveFinal = optimizer.objective;
r.objectiveDecrease = initialObjective - optimizer.objective;
r.finalRelativeStep = optimizer.history.relativeStep(end);
r.finalRelativeObjective = optimizer.history.relativeObjective(end);
r.finalStepSize = optimizer.history.stepSize(end);
r.proxGradientResidualNorm = stationarity.proxGradientResidualNorm;
r.relativeProxGradientResidual = stationarity.relativeProxGradientResidual;
r.exactZeroLayers = stats.nExactZeroLayers;
r.finalActiveLayers = delivery.nActiveLayers;
r.finalActiveSpots = delivery.nActiveSpots;
r.estimatedBDT = delivery.T_total;
r.bdtChangeFromZeroPercent = pct(delivery.T_total, zeroDelivery.T_total);
r.targetD95 = dose.summary.targetD95Min; r.targetD98 = dose.summary.targetD98Min;
r.worstAdverseDoseChangePercent = worstAdverse; r.doseQualityPass = dosePass;
r.maxConstraintViolation = optimizer.parts.maxConstraintViolation;
r.supportPreserved = supportPreserved && extracted.nLayers == delivery.nActiveLayers;
r.finalIpoptStatus = finalSolver.statusMessage; r.runStatus = 'completed';
end

function r = failed_controlled_record(alpha, lambda, err)
r = empty_controlled_record(); r.alpha = alpha; r.lambdaEL = lambda;
r.runStatus = 'failed'; r.errorIdentifier = err.identifier; r.errorMessage = err.message;
end

function token = alpha_token(alpha)
token = strrep(sprintf('%.4g', alpha), '.', 'p'); token = strrep(token, '-', 'm');
end

function value = pct(x, ref)
value = 100 * (double(x) - double(ref)) / max(abs(double(ref)), eps);
end

function merged = merge_structs_local(base, override)
merged = base; names = fieldnames(override);
for i = 1:numel(names)
    key = names{i};
    if isfield(merged, key) && isstruct(merged.(key)) && isstruct(override.(key))
        merged.(key) = merge_structs_local(merged.(key), override.(key));
    else
        merged.(key) = override.(key);
    end
end
end
