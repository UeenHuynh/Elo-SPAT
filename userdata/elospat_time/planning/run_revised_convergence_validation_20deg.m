function outputs = run_revised_convergence_validation_20deg(outputDir)
% Validate revised stationarity-aware stopping on representative support.
%
% Historical jobs 21878/21891 are read-only inputs. The objective and exact
% group-half proximal operator are unchanged. Only layer sparsity is active.

narginchk(0, 1);
root = fileparts(fileparts(mfilename('fullpath')));
add_module_paths(root);
if nargin < 1 || isempty(outputDir)
    outputDir = fullfile(root, 'results', ...
        ['revised_convergence_validation_20deg_' datestr(now, 'yyyymmdd_HHMMSS')]);
end
create_empty_output_directory(outputDir);

artifactFile = fullfile(root, 'results', 'task25_27_pipeline_20deg.mat');
baselineFile = fullfile(root, 'results', ...
    'support_representativeness_20deg_job21861', ...
    'support_86400_checkpoint.mat');
historicalFile = fullfile(root, 'results', ...
    'representative_layer_sweep_20deg_job21878', ...
    'lambda_alpha_0p01_checkpoint.mat');
source = load(artifactFile, 'cfg', 'machine', 'denseResult', ...
    'dij', 'pln', 'map', 'groups');
base = load(baselineFile, 'supportMask', 'final', 'finalDelivery');
historical = load(historicalFile, 'record');
patient = load(source.cfg.patientFile, 'ct', 'cst');

cfg = merge_structs_local(default_config(), source.cfg);
cfg.mode = 'integrated'; cfg.robust.enabled = false;
cfg.integrated.componentStage = 'layer';
cfg.lambdaB = 0; cfg.lambdaS = 0; cfg.lambdaES = 0; cfg.lambdaT = 0;
cfg.integrated.constraintPenalty = 1e3;
cfg.integrated.maxIterations = 200;
cfg.integrated.minimumIterations = 10;
cfg.integrated.stationarityTolerance = 1e-7;
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
denseDose = dose_metrics(D, denseW, patient, cfg);
denseDelivery = compute_delivery_metrics(denseW, map, groupsFull, cfg, machine);
zeroDelivery = base.finalDelivery;

scale = historical.record.lambdaEL / historical.record.alpha;
alphaGrid = [0.01; 1];
records = repmat(empty_record(), numel(alphaGrid), 1);
save(fullfile(outputDir, 'revised_convergence_setup.mat'), ...
    'alphaGrid', 'scale', 'cfg', 'validationMask', 'reducedToFull', ...
    'artifactFile', 'baselineFile', 'historicalFile', '-v7.3');

for iTrial = 1:numel(alphaGrid)
    alpha = alphaGrid(iTrial);
    cfgTrial = cfg; cfgTrial.lambdaEL = alpha * scale;
    token = alpha_token(alpha);
    fprintf('\n=== Revised convergence validation alpha=%g ===\n', alpha);
    try
        [initialObjective, initialParts] = evaluate_integrated_objective( ...
            w0, doseModel, context, cfgTrial);
        optimizer = run_integrated_optimizer( ...
            doseModel, context, w0, configure_integrated_stage(cfgTrial));
        trajectory = trajectory_table(initialObjective, initialParts, optimizer);
        objectiveDescentPass = all(diff(trajectory.objective) <= ...
            1e-10 * max(1, max(abs(trajectory.objective))));
        [sparsityStable, tailObjectiveChange] = tail_stability( ...
            optimizer.history, 10);
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
        [dosePass, doseTable, worstAdverse] = dose_quality_detail( ...
            finalDose, denseDose, cfg.evaluation.ctvD95TolerancePercent);
        record = make_record(alpha, cfgTrial, optimizer, stats, extracted, ...
            finalDelivery, zeroDelivery, denseDelivery, finalDose, ...
            finalSolver, supportPreserved, objectiveDescentPass, ...
            sparsityStable, tailObjectiveChange, dosePass, worstAdverse, ...
            initialObjective);
        records(iTrial) = record;
        writetable(trajectory, fullfile(outputDir, ...
            sprintf('alpha_%s_trajectory.csv', token)));
        writetable(layer_group_table(stats), fullfile(outputDir, ...
            sprintf('alpha_%s_layer_groups.csv', token)));
        writetable(doseTable, fullfile(outputDir, ...
            sprintf('alpha_%s_final_ipopt_dose.csv', token)));
        save(fullfile(outputDir, sprintf('alpha_%s_checkpoint.mat', token)), ...
            'record', 'optimizer', 'stats', 'extracted', 'supportFull', ...
            'final', 'finalDose', 'finalDelivery', 'finalSolver', ...
            'trajectory', 'doseTable', 'cfgTrial', '-v7.3');
    catch trialError
        records(iTrial) = failed_record(alpha, cfgTrial.lambdaEL, trialError);
        record = records(iTrial); %#ok<NASGU>
        errorInfo = struct('identifier', trialError.identifier, ...
            'message', trialError.message); %#ok<NASGU>
        save(fullfile(outputDir, sprintf('alpha_%s_checkpoint.mat', token)), ...
            'record', 'errorInfo', 'cfgTrial', '-v7.3');
        warning('Revised alpha=%g failed: %s. Continuing.', ...
            alpha, trialError.message);
    end
    writetable(struct2table(records(1:iTrial), 'AsArray', true), ...
        fullfile(outputDir, 'revised_convergence_summary.csv'));
end

validationPassed = all([records.validationPass]);
save(fullfile(outputDir, 'revised_convergence_validation_complete.mat'), ...
    'validationPassed', 'records', 'alphaGrid', 'scale', ...
    'denseDose', 'denseDelivery', 'zeroDelivery', 'cfg', ...
    'artifactFile', 'baselineFile', 'historicalFile', '-v7.3');
outputs = struct('validationPassed', validationPassed, ...
    'records', records, 'outputDir', outputDir);
fprintf('\nREVISED CONVERGENCE VALIDATION: %s\n', pass_fail(validationPassed));
fprintf('No operating point is selected by this validation runner.\n');
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

function m = dose_metrics(D, w, patient, cfg)
m = compute_dvh_metrics(full(D * reshape(double(w), [], 1)), ...
    patient.ct, patient.cst, cfg);
end

function [stable, relativeTailChange] = tail_stability(history, window)
n = numel(history.objective);
if n < window
    stable = false; relativeTailChange = Inf; return;
end
ix = (n - window + 1):n;
stable = numel(unique(history.exactZeroLayers(ix))) == 1;
relativeTailChange = abs(history.objective(ix(1)) - history.objective(ix(end))) / ...
    max(1, abs(history.objective(ix(1))));
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

function tbl = trajectory_table(initialObjective, initialParts, optimizer)
h = optimizer.history; n = optimizer.iterations + 1; iteration = (0:n-1)';
objective = [initialObjective; h.objective(:)];
doseObjective = [initialParts.dose; h.doseObjective(:)];
layerPenalty = [initialParts.layerL21half; h.layerPenalty(:)];
objectiveDecrease = [NaN; h.objectiveDecrease(:)];
stepSize = [NaN; h.stepSize(:)];
relativeStep = [NaN; h.relativeStep(:)];
relativeObjective = [NaN; h.relativeObjective(:)];
proxGradientResidualNorm = [NaN; h.proxGradientResidualNorm(:)];
normalizedStationarity = [NaN; h.normalizedStationarity(:)];
exactZeroLayers = [NaN; h.exactZeroLayers(:)];
activeSpots = [NaN; h.activeSpots(:)];
maxConstraintViolation = [initialParts.maxConstraintViolation; ...
    h.maxConstraintViolation(:)];
tbl = table(iteration, objective, doseObjective, layerPenalty, ...
    objectiveDecrease, stepSize, relativeStep, relativeObjective, ...
    proxGradientResidualNorm, normalizedStationarity, exactZeroLayers, ...
    activeSpots, maxConstraintViolation);
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

function r = empty_record()
r = struct('alpha', NaN, 'lambdaEL', NaN, 'iterations', NaN, ...
    'integratedStatus', '', 'stoppingReason', '', ...
    'stationarityTolerance', NaN, 'finalNormalizedStationarity', NaN, ...
    'finalStepSize', NaN, 'objectiveInitial', NaN, 'objectiveFinal', NaN, ...
    'objectiveDecrease', NaN, 'objectiveDescentPass', false, ...
    'tailObjectiveRelativeChange', NaN, 'sparsityStableLast10', false, ...
    'exactZeroLayers', NaN, 'finalActiveLayers', NaN, ...
    'finalActiveSpots', NaN, 'estimatedBDT', NaN, ...
    'bdtChangeFromZeroPercent', NaN, 'bdtChangeFromDensePercent', NaN, ...
    'targetD95', NaN, 'targetD98', NaN, ...
    'worstAdverseDoseChangePercent', NaN, 'doseQualityPass', false, ...
    'maxConstraintViolation', NaN, 'constraintPass', false, ...
    'supportPreserved', false, 'finalIpoptStatus', '', ...
    'finalIpoptAcceptable', false, 'validationPass', false, ...
    'runStatus', 'pending', 'errorIdentifier', '', 'errorMessage', '');
end

function r = make_record(alpha, cfg, optimizer, stats, extracted, delivery, ...
    zeroDelivery, denseDelivery, dose, finalSolver, supportPreserved, ...
    descentPass, sparsityStable, tailObjectiveChange, dosePass, worstAdverse, ...
    initialObjective)
r = empty_record();
r.alpha = alpha; r.lambdaEL = cfg.lambdaEL;
r.iterations = optimizer.iterations; r.integratedStatus = optimizer.status;
r.stoppingReason = optimizer.stoppingReason;
r.stationarityTolerance = optimizer.stationarityTolerance;
r.finalNormalizedStationarity = optimizer.finalNormalizedStationarity;
r.finalStepSize = optimizer.history.stepSize(end);
r.objectiveInitial = initialObjective;
r.objectiveFinal = optimizer.objective;
r.objectiveDecrease = r.objectiveInitial - r.objectiveFinal;
r.objectiveDescentPass = descentPass;
r.tailObjectiveRelativeChange = tailObjectiveChange;
r.sparsityStableLast10 = sparsityStable;
r.exactZeroLayers = stats.nExactZeroLayers;
r.finalActiveLayers = delivery.nActiveLayers;
r.finalActiveSpots = delivery.nActiveSpots;
r.estimatedBDT = delivery.T_total;
r.bdtChangeFromZeroPercent = pct(delivery.T_total, zeroDelivery.T_total);
r.bdtChangeFromDensePercent = pct(delivery.T_total, denseDelivery.T_total);
r.targetD95 = dose.summary.targetD95Min; r.targetD98 = dose.summary.targetD98Min;
r.worstAdverseDoseChangePercent = worstAdverse; r.doseQualityPass = dosePass;
r.maxConstraintViolation = optimizer.parts.maxConstraintViolation;
r.constraintPass = r.maxConstraintViolation <= 1e-4;
r.supportPreserved = supportPreserved && extracted.nLayers == delivery.nActiveLayers;
r.finalIpoptStatus = finalSolver.statusMessage;
r.finalIpoptAcceptable = finalSolver.converged || ...
    finalSolver.acceptableButNotConverged || ...
    contains(lower(finalSolver.statusMessage), 'acceptable');
r.validationPass = strcmp(optimizer.status, 'converged') && ...
    strcmp(optimizer.stoppingReason, 'normalized_proximal_gradient_tolerance') && ...
    optimizer.finalNormalizedStationarity <= optimizer.stationarityTolerance && ...
    descentPass && sparsityStable && r.constraintPass && r.supportPreserved && ...
    dosePass && r.finalIpoptAcceptable;
r.runStatus = 'completed';
end

function r = failed_record(alpha, lambda, err)
r = empty_record(); r.alpha = alpha; r.lambdaEL = lambda;
r.runStatus = 'failed'; r.errorIdentifier = err.identifier; r.errorMessage = err.message;
end

function token = alpha_token(alpha)
token = strrep(sprintf('%.4g', alpha), '.', 'p'); token = strrep(token, '-', 'm');
end

function value = pct(x, ref)
value = 100 * (double(x) - double(ref)) / max(abs(double(ref)), eps);
end

function label = pass_fail(value)
if value, label = 'PASS'; else, label = 'FAIL'; end
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
