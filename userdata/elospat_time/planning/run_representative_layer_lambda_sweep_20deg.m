function outputs = run_representative_layer_lambda_sweep_20deg(outputDir, baselineDir)
% Nonzero-lambda layer validation on the validated 4800-spots/CP support.
%
% Reuses the preserved Dense Arc dij and the zero-lambda final-IPOPT result
% from the support-representativeness audit. The dose matrix is not rebuilt.
% Every lambda trial runs the unchanged integrated group-half optimizer,
% extracts a fixed support, and performs support-constrained final IPOPT.

narginchk(0, 2);
root = fileparts(fileparts(mfilename('fullpath')));
add_module_paths(root);
if nargin < 1 || isempty(outputDir)
    outputDir = fullfile(root, 'results', ...
        ['representative_layer_sweep_20deg_' datestr(now, 'yyyymmdd_HHMMSS')]);
end
if nargin < 2 || isempty(baselineDir)
    baselineDir = fullfile(root, 'results', ...
        'support_representativeness_20deg_job21861');
end
create_empty_output_directory(outputDir);

artifactFile = fullfile(root, 'results', 'task25_27_pipeline_20deg.mat');
baselineFile = fullfile(baselineDir, 'support_86400_checkpoint.mat');
if exist(baselineFile, 'file') ~= 2
    error('run_representative_layer_lambda_sweep_20deg:missingBaseline', ...
        'Missing validated 4800/CP checkpoint: %s', baselineFile);
end

fprintf('Loading preserved Dense Arc artifact (no dose influence rebuild)...\n');
source = load(artifactFile, 'cfg', 'machine', 'denseResult', ...
    'dij', 'pln', 'map', 'groups');
base = load(baselineFile, 'supportMask', 'selectionDiagnostics', ...
    'final', 'finalDose', 'finalDelivery', 'solverMetrics', 'record');
patient = load(source.cfg.patientFile, 'ct', 'cst');
cfg = merge_structs_local(default_config(), source.cfg);
cfg.mode = 'integrated';
cfg.robust.enabled = false;
cfg.integrated.componentStage = 'layer';
cfg.lambdaB = 0; cfg.lambdaS = 0; cfg.lambdaES = 0; cfg.lambdaT = 0;
cfg.integrated.constraintPenalty = 1e3;
cfg.integrated.maxIterations = 80;
cfg.integrated.relativeStepTolerance = 1e-5;
cfg.integrated.relativeObjectiveTolerance = 1e-6;
cfg.integrated.runGradientCheck = false;
cfg.integrated.verbose = true;
cfg.captureConsoleOutput = false;
cfg.enforceSelectedMinimumWeight = cfg.integrated.enforceFinalMinimumWeight;

dij = source.dij;
pln = source.pln;
map = source.map;
groupsFull = source.groups;
machine = source.machine;
denseW = reshape(double(source.denseResult.w), [], 1);
nFull = double(dij.totalNumOfBixels);
supportMask = reshape(logical(base.supportMask), [], 1);
expectedSupportSpots = double(base.record.supportSpots);
if numel(supportMask) ~= nFull || nnz(supportMask) ~= expectedSupportSpots
    error('run_representative_layer_lambda_sweep_20deg:invalidBaselineSupport', ...
        'Expected the saved validated support (%d spots), found %d spots.', ...
        expectedSupportSpots, nnz(supportMask));
end
if ~base.record.representativenessPass || ~base.record.supportPreserved
    error('run_representative_layer_lambda_sweep_20deg:invalidBaseline', ...
        'Saved zero-lambda support did not pass representativeness/preservation.');
end
if any(abs(base.final.w(~supportMask)) > 1e-10)
    error('run_representative_layer_lambda_sweep_20deg:baselineLeakage', ...
        'Saved zero-lambda solution leaks outside its support.');
end

% Build the validated reduced problem once. No full sparse dose matrix is
% converted to dense, and all trials use exactly the same column mapping.
[dijReduced, reducedToFull, fullToReduced, mapReduced] = ...
    reduce_dij_to_support(dij, supportMask, map);
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
denseDose = dose_metrics(D, denseW, patient, cfg);
denseDelivery = compute_delivery_metrics( ...
    denseW, map, groupsFull, cfg, machine);
zeroDose = base.finalDose;
zeroDelivery = base.finalDelivery;

cfgScale = cfg;
cfgScale.lambdaEL = 0;
[initialObjective, initialParts] = evaluate_integrated_objective( ...
    w0, doseModel, context, cfgScale);
scale = initialParts.dose / max(initialParts.layerL21half, eps);
alphaGrid = [1e-4; 3e-4; 1e-3; 3e-3; 1e-2; 3e-2; 1e-1; 3e-1; 1];
records = repmat(empty_record(), numel(alphaGrid), 1);
fprintf('Validated support: %d spots, %d energy layers.\n', ...
    nnz(supportMask), numel(groupsReduced));
fprintf('Zero-lambda scale Fdose/R_EL = %.10g.\n', scale);

baselineSummary = struct('supportSpots', nnz(supportMask), ...
    'supportLayers', numel(groupsReduced), 'initialObjective', initialObjective, ...
    'initialDoseObjective', initialParts.dose, ...
    'initialLayerPenalty', initialParts.layerL21half, ...
    'lambdaScale', scale, 'denseBDT', denseDelivery.T_total, ...
    'zeroBDT', zeroDelivery.T_total, ...
    'zeroIpoptStatus', base.solverMetrics.statusMessage, ...
    'zeroRepresentativenessPass', base.record.representativenessPass);
save(fullfile(outputDir, 'sweep_setup_checkpoint.mat'), ...
    'baselineSummary', 'alphaGrid', 'supportMask', 'reducedToFull', ...
    'fullToReduced', 'artifactFile', 'baselineFile', 'cfg', '-v7.3');
writetable(struct2table(baselineSummary), ...
    fullfile(outputDir, 'zero_lambda_reused_baseline.csv'));

for iTrial = 1:numel(alphaGrid)
    alpha = alphaGrid(iTrial);
    cfgTrial = cfg;
    cfgTrial.lambdaEL = alpha * scale;
    token = alpha_token(alpha);
    checkpointFile = fullfile(outputDir, ...
        sprintf('lambda_alpha_%s_checkpoint.mat', token));
    fprintf('\n=== Representative layer trial alpha=%g, lambdaEL=%.10g ===\n', ...
        alpha, cfgTrial.lambdaEL);
    try
        [trialInitialObjective, trialInitialParts] = ...
            evaluate_integrated_objective(w0, doseModel, context, cfgTrial);
        optimizer = run_integrated_optimizer( ...
            doseModel, context, w0, configure_integrated_stage(cfgTrial));
        trajectory = trajectory_table( ...
            trialInitialObjective, trialInitialParts, optimizer);
        descentTolerance = 1e-10 * max(1, max(abs(trajectory.objective)));
        objectiveDescentPass = all(diff(trajectory.objective) <= descentTolerance);
        layerStats = exact_layer_stats_local(optimizer.w, groupsReduced);
        extracted = extract_integrated_support( ...
            optimizer.w, w0, context, cfgTrial);

        supportFull = false(nFull, 1);
        supportFull(reducedToFull(extracted.mask)) = true;
        wSupportFull = zeros(nFull, 1);
        wSupportFull(reducedToFull) = extracted.w;
        if any(wSupportFull(~supportFull) ~= 0)
            error('run_representative_layer_lambda_sweep_20deg:supportEmbedding', ...
                'Reduced-to-full support embedding leaked.');
        end

        startIpopt = tic;
        final = run_final_reoptimization(dij, patient.cst, pln, ...
            wSupportFull, cfgTrial, supportFull, extracted.layerMask);
        ipoptSeconds = toc(startIpopt);
        supportPreserved = ~any(abs(final.w(~supportFull)) > 1e-10);
        finalDose = dose_metrics(D, final.w, patient, cfgTrial);
        finalDelivery = compute_delivery_metrics( ...
            final.w, map, groupsFull, cfgTrial, machine);
        finalSolver = compute_solver_metrics( ...
            final.resultGUI, struct(), final.optimizer);
        [doseQualityPass, structureTable, worstAdverse] = ...
            dose_quality_detail(finalDose, denseDose, ...
            cfg.evaluation.ctvD95TolerancePercent);
        layerTable = layer_group_table(layerStats);
        writetable(trajectory, fullfile(outputDir, ...
            sprintf('lambda_alpha_%s_objective_trajectory.csv', token)));
        writetable(layerTable, fullfile(outputDir, ...
            sprintf('lambda_alpha_%s_layer_groups.csv', token)));
        writetable(structureTable, fullfile(outputDir, ...
            sprintf('lambda_alpha_%s_final_ipopt_dose.csv', token)));

        records(iTrial) = completed_record(alpha, cfgTrial.lambdaEL, ...
            optimizer, layerStats, extracted, finalDelivery, denseDelivery, ...
            zeroDelivery, finalDose, trialInitialObjective, ...
            objectiveDescentPass, supportPreserved, finalSolver, ...
            ipoptSeconds, doseQualityPass, worstAdverse);
        record = records(iTrial); %#ok<NASGU>
        save(checkpointFile, 'record', 'optimizer', 'layerStats', ...
            'extracted', 'supportFull', 'wSupportFull', 'final', ...
            'finalDose', 'finalDelivery', 'finalSolver', 'trajectory', ...
            'layerTable', 'structureTable', 'cfgTrial', '-v7.3');
    catch trialError
        records(iTrial) = failed_record(alpha, cfgTrial.lambdaEL, trialError);
        errorInfo = struct('identifier', trialError.identifier, ...
            'message', trialError.message); %#ok<NASGU>
        record = records(iTrial); %#ok<NASGU>
        save(checkpointFile, 'record', 'errorInfo', 'cfgTrial', '-v7.3');
        warning('run_representative_layer_lambda_sweep_20deg:trialFailed', ...
            'alpha=%g failed: %s. Continuing safely.', alpha, trialError.message);
    end
    writetable(struct2table(records(1:iTrial), 'AsArray', true), ...
        fullfile(outputDir, 'lambda_sweep.csv'));
end

recordTable = struct2table(records, 'AsArray', true);
eligible = [records.doseQualityPass] & [records.objectiveDescentPass] & ...
    [records.constraintPass] & [records.supportPreserved] & ...
    [records.finalIpoptAcceptable] & ...
    [records.exactZeroLayers] > 0 & ...
    ~strcmp({records.integratedStatus}, 'line_search_failed') & ...
    ~strcmp({records.runStatus}, 'failed');
recordTable.operatingPointEligible = eligible(:);
writetable(recordTable, fullfile(outputDir, 'lambda_sweep.csv'));
save(fullfile(outputDir, 'representative_layer_sweep_complete.mat'), ...
    'recordTable', 'records', 'eligible', 'alphaGrid', 'scale', ...
    'baselineSummary', 'denseDose', 'denseDelivery', 'zeroDose', ...
    'zeroDelivery', 'artifactFile', 'baselineFile', 'cfg', '-v7.3');
outputs = struct('outputDir', outputDir, 'recordTable', recordTable, ...
    'eligibleIndices', find(eligible), ...
    'operatingPointSelected', false);
fprintf('\nSWEEP COMPLETE. Operating point selection is deferred to the recorded\n');
fprintf('dose-feasible layer/BDT trade-off; maximum sparsity is not automatic.\n');
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
    if ~isempty(entries)
        error('run_representative_layer_lambda_sweep_20deg:nonemptyOutput', ...
            'Use an empty job-specific output directory: %s', outputDir);
    end
else
    mkdir(outputDir);
end
end

function D = nominal_matrix(dij)
if iscell(dij.physicalDose), D = dij.physicalDose{1}; else, D = dij.physicalDose; end
end

function m = dose_metrics(D, w, patient, cfg)
m = compute_dvh_metrics(full(D * reshape(double(w), [], 1)), ...
    patient.ct, patient.cst, cfg);
end

function stats = exact_layer_stats_local(w, groups)
n = numel(groups);
norms = zeros(n, 1); zerosMask = false(n, 1);
for g = 1:n
    values = w(groups(g).ix);
    norms(g) = norm(values, 2);
    zerosMask(g) = all(values == 0);
end
stats = struct('nGroups', n, 'nExactZeroLayers', nnz(zerosMask), ...
    'nExactActiveLayers', n - nnz(zerosMask), 'exactZeroMask', zerosMask, ...
    'groupNorm', norms, 'groupIndex', (1:n)', ...
    'cpIndex', reshape(double([groups.cpIndex]), [], 1), ...
    'energyIndex', reshape(double([groups.energyIndex]), [], 1), ...
    'energy', reshape(double([groups.energy]), [], 1));
end

function tbl = trajectory_table(initialObjective, initialParts, optimizer)
n = optimizer.iterations + 1;
iteration = (0:n-1)';
objective = [initialObjective; optimizer.history.objective(:)];
doseObjective = [initialParts.dose; optimizer.history.doseObjective(:)];
layerPenalty = [initialParts.layerL21half; optimizer.history.layerPenalty(:)];
maxConstraintViolation = [initialParts.maxConstraintViolation; ...
    optimizer.history.maxConstraintViolation(:)];
stepSize = [NaN; optimizer.history.stepSize(:)];
tbl = table(iteration, objective, doseObjective, layerPenalty, ...
    maxConstraintViolation, stepSize);
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
        rows(i).adverseChangePercent = max([0, -rows(i).D95ChangePercent, ...
            -rows(i).D98ChangePercent]);
    elseif strcmpi(d.type, 'OAR')
        rows(i).adverseChangePercent = max([0, rows(i).DmeanChangePercent, ...
            rows(i).DmaxChangePercent]);
    else
        rows(i).adverseChangePercent = 0;
    end
    rows(i).structurePass = rows(i).adverseChangePercent <= tolerance;
    worstAdverse = max(worstAdverse, rows(i).adverseChangePercent);
end
tbl = struct2table(rows);
pass = all([rows.structurePass]);
end

function r = empty_record()
r = struct('alpha', NaN, 'lambdaEL', NaN, 'objectiveInitial', NaN, ...
    'objectiveFinal', NaN, 'doseObjectiveFinal', NaN, ...
    'exactActiveLayers', NaN, 'exactZeroLayers', NaN, ...
    'finalActiveLayers', NaN, 'finalActiveSpots', NaN, ...
    'layerReductionPercent', NaN, 'estimatedBDT', NaN, ...
    'bdtChangeFromZeroPercent', NaN, 'bdtChangeFromDensePercent', NaN, ...
    'targetD95', NaN, 'targetD98', NaN, 'oarDmaxMax', NaN, ...
    'worstAdverseDoseChangePercent', NaN, 'doseQualityPass', false, ...
    'maxConstraintViolation', NaN, 'constraintPass', false, ...
    'objectiveDescentPass', false, 'supportPreserved', false, ...
    'integratedIterations', NaN, 'integratedStatus', '', ...
    'finalIpoptStatus', '', 'finalIpoptConverged', false, ...
    'finalIpoptAcceptable', false, 'finalIpoptSeconds', NaN, ...
    'runStatus', 'pending', 'errorIdentifier', '', 'errorMessage', '');
end

function r = completed_record(alpha, lambda, optimizer, stats, extracted, ...
    delivery, denseDelivery, zeroDelivery, dose, initialObjective, ...
    descentPass, supportPreserved, finalSolver, ipoptSeconds, dosePass, ...
    worstAdverse)
r = empty_record();
r.alpha = alpha; r.lambdaEL = lambda;
r.objectiveInitial = initialObjective;
r.objectiveFinal = optimizer.objective;
r.doseObjectiveFinal = optimizer.parts.dose;
r.exactActiveLayers = stats.nExactActiveLayers;
r.exactZeroLayers = stats.nExactZeroLayers;
r.finalActiveLayers = delivery.nActiveLayers;
r.finalActiveSpots = delivery.nActiveSpots;
r.layerReductionPercent = 100 * stats.nExactZeroLayers / stats.nGroups;
r.estimatedBDT = delivery.T_total;
r.bdtChangeFromZeroPercent = pct(delivery.T_total, zeroDelivery.T_total);
r.bdtChangeFromDensePercent = pct(delivery.T_total, denseDelivery.T_total);
r.targetD95 = dose.summary.targetD95Min;
r.targetD98 = dose.summary.targetD98Min;
r.oarDmaxMax = dose.summary.oarDmaxMax;
r.worstAdverseDoseChangePercent = worstAdverse;
r.doseQualityPass = dosePass;
r.maxConstraintViolation = optimizer.parts.maxConstraintViolation;
r.constraintPass = r.maxConstraintViolation <= 1e-4;
r.objectiveDescentPass = descentPass;
r.supportPreserved = supportPreserved && extracted.nLayers == delivery.nActiveLayers;
r.integratedIterations = optimizer.iterations;
r.integratedStatus = optimizer.status;
r.finalIpoptStatus = finalSolver.statusMessage;
r.finalIpoptConverged = finalSolver.converged;
r.finalIpoptAcceptable = finalSolver.converged || ...
    finalSolver.acceptableButNotConverged || ...
    contains(lower(finalSolver.statusMessage), 'acceptable');
r.finalIpoptSeconds = ipoptSeconds;
r.runStatus = 'completed';
end

function r = failed_record(alpha, lambda, err)
r = empty_record();
r.alpha = alpha; r.lambdaEL = lambda; r.runStatus = 'failed';
r.errorIdentifier = err.identifier; r.errorMessage = err.message;
end

function token = alpha_token(alpha)
token = strrep(sprintf('%.4g', alpha), '.', 'p');
token = strrep(token, '-', 'm');
end

function value = pct(x, ref)
value = 100 * (double(x) - double(ref)) / max(abs(double(ref)), eps);
end

function merged = merge_structs_local(base, override)
merged = base;
names = fieldnames(override);
for i = 1:numel(names)
    key = names{i};
    if isfield(merged, key) && isstruct(merged.(key)) && isstruct(override.(key))
        merged.(key) = merge_structs_local(merged.(key), override.(key));
    else
        merged.(key) = override.(key);
    end
end
end
