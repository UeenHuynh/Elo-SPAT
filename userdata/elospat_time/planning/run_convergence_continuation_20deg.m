function outputs = run_convergence_continuation_20deg(outputDir)
% Continue job 21894 exactly from iteration 200 in checkpointed blocks.
%
% Layer-only validation: spot, sequencing, and time penalties remain off.
% maxIterations is interpreted as an additional budget by the resume API.

narginchk(0, 1);
root = fileparts(fileparts(mfilename('fullpath')));
add_module_paths(root);
if nargin < 1 || isempty(outputDir)
    outputDir = fullfile(root, 'results', ...
        ['convergence_continuation_20deg_' datestr(now, 'yyyymmdd_HHMMSS')]);
end
create_empty_output_directory(outputDir);

artifactFile = fullfile(root, 'results', 'task25_27_pipeline_20deg.mat');
baselineFile = fullfile(root, 'results', ...
    'support_representativeness_20deg_job21861', ...
    'support_86400_checkpoint.mat');
sourceDir = fullfile(root, 'results', ...
    'revised_convergence_validation_20deg_job21894');
source = load(artifactFile, 'cfg', 'machine', 'denseResult', ...
    'dij', 'pln', 'map', 'groups');
base = load(baselineFile, 'supportMask', 'final', 'finalDelivery');
patient = load(source.cfg.patientFile, 'ct', 'cst');

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
D = nominal_matrix(dij);
denseDose = dose_metrics(D, denseW, patient, source.cfg);
denseDelivery = compute_delivery_metrics( ...
    denseW, map, groupsFull, source.cfg, machine);
zeroDelivery = base.finalDelivery;
wReferenceFull = reshape(double(base.final.w), [], 1);
wReferenceReduced = wReferenceFull(reducedToFull);

alphaGrid = [0.01; 1];
chunkSize = 100;
additionalBudget = 400;
stationarityTolerance = 1e-7;
records = repmat(empty_record(), numel(alphaGrid), 1);
save(fullfile(outputDir, 'continuation_setup.mat'), ...
    'alphaGrid', 'chunkSize', 'additionalBudget', ...
    'stationarityTolerance', 'validationMask', 'reducedToFull', ...
    'artifactFile', 'baselineFile', 'sourceDir', '-v7.3');

for iTrial = 1:numel(alphaGrid)
    alpha = alphaGrid(iTrial); token = alpha_token(alpha);
    seedFile = fullfile(sourceDir, ...
        sprintf('alpha_%s_checkpoint.mat', token));
    seed = load(seedFile, 'optimizer', 'record', 'cfgTrial');
    cfgTrial = seed.cfgTrial;
    assert_layer_only_configuration(cfgTrial, stationarityTolerance);
    cfgTrial.integrated.maxIterations = chunkSize;
    cfgTrial.integrated.verbose = true;
    cfgTrial.captureConsoleOutput = false;
    cfgTrial = configure_integrated_stage(cfgTrial);
    doseModel = build_matrad_dose_model( ...
        dijReduced, patient.cst, pln, cfgTrial.integrated.constraintPenalty);
    context = build_integrated_context( ...
        groupsReduced, cpGroupsReduced, mapReduced, machine, cfgTrial);

    state = seed.optimizer;
    verify_seed_state(state, seed.record, cfgTrial, doseModel, context);
    startIteration = state.iterations;
    fprintf('\n=== Exact continuation alpha=%g from iteration %d ===\n', ...
        alpha, startIteration);
    try
        while state.iterations < startIteration + additionalBudget
            state = run_integrated_optimizer( ...
                doseModel, context, state.w, cfgTrial, state);
            trajectory = trajectory_table(seed.record.objectiveInitial, state);
            iter = state.iterations;
            writetable(trajectory, fullfile(outputDir, ...
                sprintf('alpha_%s_trajectory.csv', token)));
            optimizer = state; %#ok<NASGU>
            checkpointInfo = checkpoint_info(optimizer); %#ok<NASGU>
            save(fullfile(outputDir, sprintf( ...
                'alpha_%s_iter_%04d_checkpoint.mat', token, iter)), ...
                'optimizer', 'checkpointInfo', 'cfgTrial', 'seedFile', '-v7.3');
            fprintf(['CHECKPOINT alpha=%g iter=%d J=%.12g residual=%.6e ' ...
                'step=%.3e zeroLayers=%d status=%s reason=%s\n'], ...
                alpha, iter, state.objective, ...
                state.finalNormalizedStationarity, state.acceptedStep, ...
                state.history.exactZeroLayers(end), state.status, ...
                state.stoppingReason);
            if ~strcmp(state.status, 'maximum_iterations')
                break;
            end
        end

        optimizer = state;
        stats = exact_layer_stats_local(optimizer.w, groupsReduced);
        extracted = extract_integrated_support( ...
            optimizer.w, wReferenceReduced, context, cfgTrial);
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
            finalDose, denseDose, cfgTrial.evaluation.ctvD95TolerancePercent);
        record = make_record(alpha, cfgTrial, seed, optimizer, stats, ...
            extracted, finalDelivery, zeroDelivery, denseDelivery, ...
            finalDose, finalSolver, supportPreserved, dosePass, worstAdverse);
        records(iTrial) = record;
        trajectory = trajectory_table(seed.record.objectiveInitial, optimizer);
        writetable(trajectory, fullfile(outputDir, ...
            sprintf('alpha_%s_trajectory.csv', token)));
        writetable(layer_group_table(stats), fullfile(outputDir, ...
            sprintf('alpha_%s_layer_groups.csv', token)));
        writetable(doseTable, fullfile(outputDir, ...
            sprintf('alpha_%s_final_ipopt_dose.csv', token)));
        save(fullfile(outputDir, sprintf('alpha_%s_final.mat', token)), ...
            'record', 'optimizer', 'stats', 'extracted', 'supportFull', ...
            'final', 'finalDose', 'finalDelivery', 'finalSolver', ...
            'trajectory', 'doseTable', 'cfgTrial', 'seedFile', '-v7.3');
    catch trialError
        records(iTrial) = failed_record(alpha, seed.cfgTrial.lambdaEL, trialError);
        record = records(iTrial); %#ok<NASGU>
        errorInfo = struct('identifier', trialError.identifier, ...
            'message', trialError.message, 'stack', trialError.stack); %#ok<NASGU>
        save(fullfile(outputDir, sprintf('alpha_%s_failure.mat', token)), ...
            'record', 'errorInfo', 'cfgTrial', 'seedFile', '-v7.3');
        warning('Continuation alpha=%g failed: %s. Continuing.', ...
            alpha, trialError.message);
    end
    writetable(struct2table(records(1:iTrial), 'AsArray', true), ...
        fullfile(outputDir, 'continuation_summary.csv'));
end

save(fullfile(outputDir, 'continuation_complete.mat'), ...
    'records', 'alphaGrid', 'chunkSize', 'additionalBudget', ...
    'stationarityTolerance', 'denseDose', 'denseDelivery', ...
    'zeroDelivery', 'artifactFile', 'baselineFile', 'sourceDir', '-v7.3');
outputs = struct('records', records, 'outputDir', outputDir);
fprintf('\nCONVERGENCE CONTINUATION COMPLETE. No operating point selected.\n');
end

function verify_seed_state(state, record, cfg, doseModel, context)
required = {'w','objective','parts','history','iterations'};
for i = 1:numel(required), assert(isfield(state, required{i})); end
assert(state.iterations == 200, 'Expected job 21894 iteration-200 seed.');
assert(numel(state.w) == context.nWeights);
assert(numel(state.history.iteration) == state.iterations);
assert(state.history.iteration(end) == state.iterations);
assert(abs(state.history.stepSize(end) - 1) <= 1e-14);
assert(abs(cfg.integrated.stationarityTolerance - 1e-7) <= eps);
[objective, ~] = evaluate_integrated_objective( ...
    state.w, doseModel, context, cfg);
relativeError = abs(objective - state.objective) / max(1, abs(state.objective));
assert(relativeError <= 1e-11, ...
    'Seed objective mismatch (relative error %.3e).', relativeError);
assert(abs(record.objectiveFinal - state.objective) <= ...
    1e-10 * max(1, abs(state.objective)));
end

function assert_layer_only_configuration(cfg, tolerance)
assert(cfg.lambdaEL > 0, 'lambdaEL must remain nonzero.');
assert(cfg.lambdaB == 0 && cfg.lambdaS == 0 && ...
    cfg.lambdaES == 0 && cfg.lambdaT == 0, ...
    'Only energy-layer regularization is allowed in this validation.');
assert(strcmpi(cfg.integrated.componentStage, 'layer'));
assert(abs(cfg.integrated.stationarityTolerance - tolerance) <= eps);
end

function info = checkpoint_info(optimizer)
info = struct('iteration', optimizer.iterations, ...
    'objective', optimizer.objective, ...
    'normalizedStationarity', optimizer.finalNormalizedStationarity, ...
    'absoluteProxGradientResidual', ...
        optimizer.history.proxGradientResidualNorm(end), ...
    'acceptedStep', optimizer.acceptedStep, ...
    'exactZeroLayers', optimizer.history.exactZeroLayers(end), ...
    'activeLayers', optimizer.history.exactActiveLayers(end), ...
    'activeSpots', optimizer.history.activeSpots(end), ...
    'status', optimizer.status, 'stoppingReason', optimizer.stoppingReason, ...
    'resumeObjectiveRelativeError', optimizer.resumeObjectiveRelativeError);
end

function tbl = trajectory_table(initialObjective, optimizer)
h = optimizer.history;
iteration = [0; h.iteration(:)];
objective = [initialObjective; h.objective(:)];
doseObjective = [NaN; h.doseObjective(:)];
layerPenalty = [NaN; h.layerPenalty(:)];
objectiveDecrease = [NaN; h.objectiveDecrease(:)];
stepSize = [NaN; h.stepSize(:)];
relativeStep = [NaN; h.relativeStep(:)];
relativeObjective = [NaN; h.relativeObjective(:)];
proxGradientResidualNorm = [NaN; h.proxGradientResidualNorm(:)];
normalizedStationarity = [NaN; h.normalizedStationarity(:)];
exactZeroLayers = [NaN; h.exactZeroLayers(:)];
exactActiveLayers = [NaN; h.exactActiveLayers(:)];
activeSpots = [NaN; h.activeSpots(:)];
maxConstraintViolation = [NaN; h.maxConstraintViolation(:)];
tbl = table(iteration, objective, doseObjective, layerPenalty, ...
    objectiveDecrease, stepSize, relativeStep, relativeObjective, ...
    proxGradientResidualNorm, normalizedStationarity, exactZeroLayers, ...
    exactActiveLayers, activeSpots, maxConstraintViolation);
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

function tbl = layer_group_table(stats)
tbl = table(stats.groupIndex, stats.cpIndex, stats.energyIndex, stats.energy, ...
    stats.groupNorm, stats.exactZeroMask, 'VariableNames', ...
    {'groupIndex','cpIndex','energyIndex','energy','groupNorm','exactZero'});
end

function r = empty_record()
r = struct('alpha', NaN, 'lambdaEL', NaN, 'startIteration', NaN, ...
    'iterations', NaN, 'segmentIterations', NaN, 'integratedStatus', '', ...
    'stoppingReason', '', 'stationarityTolerance', NaN, ...
    'finalNormalizedStationarity', NaN, ...
    'absoluteProxGradientResidual', NaN, 'finalStepSize', NaN, ...
    'objectiveAtIteration200', NaN, 'objectiveFinal', NaN, ...
    'continuationObjectiveDecrease', NaN, 'objectiveDescentPass', false, ...
    'exactZeroLayers', NaN, 'integratedActiveLayers', NaN, ...
    'finalActiveLayers', NaN, 'finalActiveSpots', NaN, ...
    'estimatedBDT', NaN, 'bdtChangeFromZeroPercent', NaN, ...
    'bdtChangeFromDensePercent', NaN, 'targetD95', NaN, 'targetD98', NaN, ...
    'worstAdverseDoseChangePercent', NaN, 'doseQualityPass', false, ...
    'maxConstraintViolation', NaN, 'constraintPass', false, ...
    'supportPreserved', false, 'finalIpoptStatus', '', ...
    'finalIpoptAcceptable', false, 'runStatus', 'pending', ...
    'errorIdentifier', '', 'errorMessage', '');
end

function r = make_record(alpha, cfg, seed, optimizer, stats, extracted, ...
    delivery, zeroDelivery, denseDelivery, dose, finalSolver, ...
    supportPreserved, dosePass, worstAdverse)
r = empty_record();
r.alpha = alpha; r.lambdaEL = cfg.lambdaEL;
r.startIteration = seed.optimizer.iterations;
r.iterations = optimizer.iterations;
r.segmentIterations = optimizer.iterations - seed.optimizer.iterations;
r.integratedStatus = optimizer.status; r.stoppingReason = optimizer.stoppingReason;
r.stationarityTolerance = optimizer.stationarityTolerance;
r.finalNormalizedStationarity = optimizer.finalNormalizedStationarity;
r.absoluteProxGradientResidual = optimizer.history.proxGradientResidualNorm(end);
r.finalStepSize = optimizer.history.stepSize(end);
r.objectiveAtIteration200 = seed.optimizer.objective;
r.objectiveFinal = optimizer.objective;
r.continuationObjectiveDecrease = seed.optimizer.objective - optimizer.objective;
r.objectiveDescentPass = all(diff(optimizer.history.objective) <= ...
    1e-10 * max(1, max(abs(optimizer.history.objective))));
r.exactZeroLayers = stats.nExactZeroLayers;
r.integratedActiveLayers = stats.nExactActiveLayers;
r.finalActiveLayers = delivery.nActiveLayers;
r.finalActiveSpots = delivery.nActiveSpots;
r.estimatedBDT = delivery.T_total;
r.bdtChangeFromZeroPercent = pct(delivery.T_total, zeroDelivery.T_total);
r.bdtChangeFromDensePercent = pct(delivery.T_total, denseDelivery.T_total);
r.targetD95 = dose.summary.targetD95Min; r.targetD98 = dose.summary.targetD98Min;
r.worstAdverseDoseChangePercent = worstAdverse; r.doseQualityPass = dosePass;
r.maxConstraintViolation = optimizer.parts.maxConstraintViolation;
r.constraintPass = r.maxConstraintViolation <= 1e-4;
r.supportPreserved = supportPreserved;
r.finalIpoptStatus = finalSolver.statusMessage;
r.finalIpoptAcceptable = finalSolver.converged || ...
    finalSolver.acceptableButNotConverged || ...
    contains(lower(finalSolver.statusMessage), 'acceptable');
r.runStatus = 'completed';
end

function r = failed_record(alpha, lambda, err)
r = empty_record(); r.alpha = alpha; r.lambdaEL = lambda;
r.runStatus = 'failed'; r.errorIdentifier = err.identifier; r.errorMessage = err.message;
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

function add_module_paths(root)
folders = {'config','planning','mapping','objective','delivery','reduction', ...
    'robustness','evaluation','integrated'};
for i = 1:numel(folders), addpath(fullfile(root, folders{i})); end
end

function create_empty_output_directory(outputDir)
if exist(outputDir, 'dir') == 7
    entries = dir(outputDir); entries = entries(~ismember({entries.name}, {'.','..'}));
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

function token = alpha_token(alpha)
token = strrep(sprintf('%.4g', alpha), '.', 'p'); token = strrep(token, '-', 'm');
end

function value = pct(x, ref)
value = 100 * (double(x) - double(ref)) / max(abs(double(ref)), eps);
end
