function outputs = run_reduced_support_representativeness_20deg(outputDir, priorValidationDir)
% Diagnose whether reduced-support construction causes dose mismatch.
%
% This runner does not invoke the integrated optimizer. It increases a
% nested, layer-stratified Dense-reference support and performs only locked-
% support matRad/IPOPT reoptimization. Nonzero sparsity, sequencing and time
% penalties are disabled.

narginchk(0, 2);
workspaceDir = fileparts(fileparts(mfilename('fullpath')));
add_module_paths(workspaceDir);
if nargin < 1 || isempty(outputDir)
    outputDir = fullfile(workspaceDir, 'results', ...
        ['support_representativeness_20deg_' datestr(now, 'yyyymmdd_HHMMSS')]);
end
if nargin < 2 || isempty(priorValidationDir)
    priorValidationDir = fullfile(workspaceDir, 'results', ...
        'integrated_layer_validation_20deg_job21835');
end
create_empty_output_directory(outputDir);

artifactFile = fullfile(workspaceDir, 'results', 'task25_27_pipeline_20deg.mat');
source = load(artifactFile, 'cfg', 'machine', 'denseResult', ...
    'dij', 'pln', 'stf', 'map', 'groups');
patient = load(source.cfg.patientFile, 'ct', 'cst');
cfg = merge_structs_local(default_config(), source.cfg);
cfg.robust.enabled = false;
cfg.lambdaEL = 0; cfg.lambdaB = 0; cfg.lambdaS = 0;
cfg.lambdaES = 0; cfg.lambdaT = 0;
cfg.enforceSelectedMinimumWeight = cfg.integrated.enforceFinalMinimumWeight;

dij = source.dij;
pln = source.pln;
map = source.map;
groups = source.groups;
machine = source.machine;
denseW = reshape(double(source.denseResult.w), [], 1);
D = nominal_matrix(dij);
denseDose = dose_metrics(D, denseW, patient, cfg);
denseDelivery = compute_delivery_metrics(denseW, map, groups, cfg, machine);
tolerance = cfg.evaluation.ctvD95TolerancePercent;

% Attribute the existing 1200/CP failure using all three zero-lambda stages.
priorCheckpoint = fullfile(priorValidationDir, 'zero_lambda_checkpoint.mat');
if exist(priorCheckpoint, 'file') ~= 2
    error('run_reduced_support_representativeness_20deg:missingPriorCheckpoint', ...
        'Missing prior zero-lambda checkpoint: %s', priorCheckpoint);
end
prior = load(priorCheckpoint, 'zero');
zero = prior.zero;
[mask1200, diagnostics1200] = build_layer_stratified_support(denseW, map, 1200);
if ~isequal(mask1200, zero.validationMask)
    error('run_reduced_support_representativeness_20deg:selectorMismatch', ...
        'Refactored 1200/CP support does not reproduce job 21835.');
end
stageDose = struct();
stageDose.dense = denseDose;
stageDose.restrictedInitial = dose_metrics(D, zero.wInitialFull, patient, cfg);
stageDose.lambdaZeroOptimizer = dose_metrics(D, zero.wOptimizerFull, patient, cfg);
stageDose.finalIpopt = dose_metrics(D, zero.final.w, patient, cfg);
stageTable = build_stage_comparison(stageDose);
writetable(stageTable, fullfile(outputDir, 'support_01200_stage_comparison.csv'));
writetable(struct2table(diagnostics1200), fullfile(outputDir, ...
    'support_01200_selection_by_cp.csv'));

caps = [1800; 2400; 3600; 4800];
records = repmat(empty_record(), numel(caps), 1);
previousMask = mask1200;
representativeIndex = NaN;
for i = 1:numel(caps)
    cap = caps(i);
    fprintf('\n=== Zero-lambda support representativeness: %d spots/CP ===\n', cap);
    [supportMask, selectionDiagnostics] = ...
        build_layer_stratified_support(denseW, map, cap);
    if any(previousMask & ~supportMask)
        error('run_reduced_support_representativeness_20deg:notNested', ...
            'Support at %d spots/CP is not a superset of the preceding support.', cap);
    end
    previousMask = supportMask;
    wInit = zeros(size(denseW));
    wInit(supportMask) = denseW(supportMask);
    initialDose = dose_metrics(D, wInit, patient, cfg);
    selectionTable = struct2table(selectionDiagnostics);
    writetable(selectionTable, fullfile(outputDir, ...
        sprintf('support_%05d_selection_by_cp.csv', nnz(supportMask))));

    try
        layerMask = support_layer_mask(supportMask, groups);
        startTime = tic;
        final = run_final_reoptimization(dij, patient.cst, pln, wInit, ...
            cfg, supportMask, layerMask);
        solverSeconds = toc(startTime);
        supportPreserved = ~any(abs(final.w(~supportMask)) > 1e-10);
        finalDose = dose_metrics(D, final.w, patient, cfg);
        finalDelivery = compute_delivery_metrics( ...
            final.w, map, groups, cfg, machine);
        solverMetrics = compute_solver_metrics( ...
            final.resultGUI, struct(), final.optimizer);
        representativenessPass = quality_pass_detail( ...
            finalDose, denseDose, tolerance);
        structureTable = build_structure_comparison( ...
            denseDose, initialDose, finalDose, tolerance);
        writetable(structureTable, fullfile(outputDir, ...
            sprintf('support_%05d_structure_comparison.csv', nnz(supportMask))));

        records(i) = make_record(cap, supportMask, selectionDiagnostics, ...
            denseW, initialDose, finalDose, denseDose, finalDelivery, denseDelivery, ...
            representativenessPass, supportPreserved, solverMetrics, solverSeconds);
        record = records(i); %#ok<NASGU>
        save(fullfile(outputDir, sprintf('support_%05d_checkpoint.mat', ...
            nnz(supportMask))), 'record', 'supportMask', 'selectionDiagnostics', ...
            'wInit', 'final', 'initialDose', 'finalDose', 'finalDelivery', ...
            'solverMetrics', 'structureTable', '-v7.3');
    catch trialError
        records(i) = failed_record(cap, supportMask, trialError);
        errorInfo = struct('identifier', trialError.identifier, ...
            'message', trialError.message); %#ok<NASGU>
        record = records(i); %#ok<NASGU>
        save(fullfile(outputDir, sprintf('support_%05d_checkpoint.mat', ...
            nnz(supportMask))), 'record', 'supportMask', ...
            'selectionDiagnostics', 'wInit', 'errorInfo', '-v7.3');
        warning('run_reduced_support_representativeness_20deg:trialFailed', ...
            '%d spots/CP failed: %s. Continuing safely.', cap, trialError.message);
    end
    writetable(struct2table(records(1:i), 'AsArray', true), ...
        fullfile(outputDir, 'support_representativeness_summary.csv'));
    if records(i).representativenessPass && records(i).supportPreserved
        representativeIndex = i;
        fprintf('Representativeness gate passed at %d spots/CP. Stopping sweep.\n', cap);
        break;
    end
end

testedRecords = records(1:i);
passed = isfinite(representativeIndex);
save(fullfile(outputDir, 'support_representativeness_audit.mat'), ...
    'passed', 'representativeIndex', 'testedRecords', 'caps', ...
    'stageDose', 'stageTable', 'denseDose', 'denseDelivery', ...
    'artifactFile', 'priorCheckpoint', 'cfg', '-v7.3');
outputs = struct('passed', passed, 'representativeIndex', representativeIndex, ...
    'records', testedRecords, 'outputDir', outputDir);
if passed
    fprintf('\nREDUCED SUPPORT REPRESENTATIVENESS PASSED.\n');
else
    fprintf('\nREDUCED SUPPORT REPRESENTATIVENESS NOT YET ACHIEVED.\n');
end
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
        error('run_reduced_support_representativeness_20deg:nonemptyOutput', ...
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

function mask = support_layer_mask(supportMask, groups)
mask = false(numel(groups), 1);
for g = 1:numel(groups), mask(g) = any(supportMask(groups(g).ix)); end
end

function r = empty_record()
r = struct('spotsPerCP', NaN, 'supportSpots', NaN, ...
    'selectedWeightFraction', NaN, 'minimumLayersPerCP', NaN, ...
    'initialD95', NaN, 'initialD98', NaN, 'initialLtFemoralDmax', NaN, ...
    'finalD95', NaN, 'finalD98', NaN, 'finalLtFemoralDmean', NaN, ...
    'finalLtFemoralDmax', NaN, 'ltFemoralDmaxChangePercent', NaN, ...
    'activeLayers', NaN, 'activeSpots', NaN, 'estimatedBDT', NaN, ...
    'bdtChangePercent', NaN, ...
    'representativenessPass', false, 'supportPreserved', false, ...
    'solverStatus', '', 'solverSeconds', NaN, 'status', 'pending', ...
    'errorIdentifier', '', 'errorMessage', '');
end

function r = make_record(cap, supportMask, selectionDiagnostics, denseW, ...
    initialDose, finalDose, denseDose, finalDelivery, denseDelivery, ...
    representativenessPass, supportPreserved, solverMetrics, solverSeconds)
r = empty_record();
r.spotsPerCP = cap;
r.supportSpots = nnz(supportMask);
r.selectedWeightFraction = sum(max(denseW(supportMask), 0)) / ...
    max(sum(max(denseW, 0)), eps);
r.minimumLayersPerCP = min([selectionDiagnostics.selectedLayers]);
r.initialD95 = initialDose.summary.targetD95Min;
r.initialD98 = initialDose.summary.targetD98Min;
r.initialLtFemoralDmax = structure_metric(initialDose, 'Lt femoral head', 'Dmax');
r.finalD95 = finalDose.summary.targetD95Min;
r.finalD98 = finalDose.summary.targetD98Min;
r.finalLtFemoralDmean = structure_metric(finalDose, 'Lt femoral head', 'Dmean');
r.finalLtFemoralDmax = structure_metric(finalDose, 'Lt femoral head', 'Dmax');
denseLtDmax = structure_metric(denseDose, 'Lt femoral head', 'Dmax');
r.ltFemoralDmaxChangePercent = percent_change(r.finalLtFemoralDmax, denseLtDmax);
r.activeLayers = finalDelivery.nActiveLayers;
r.activeSpots = finalDelivery.nActiveSpots;
r.estimatedBDT = finalDelivery.T_total;
r.bdtChangePercent = percent_change(finalDelivery.T_total, denseDelivery.T_total);
r.representativenessPass = representativenessPass;
r.supportPreserved = supportPreserved;
r.solverStatus = solverMetrics.statusMessage;
r.solverSeconds = solverSeconds;
r.status = 'completed';
end

function r = failed_record(cap, supportMask, trialError)
r = empty_record();
r.spotsPerCP = cap;
r.supportSpots = nnz(supportMask);
r.status = 'failed';
r.errorIdentifier = trialError.identifier;
r.errorMessage = trialError.message;
end

function value = structure_metric(dose, name, field)
names = {dose.structures.name};
ix = find(strcmpi(names, name), 1);
if isempty(ix), error('Missing structure %s.', name); end
value = dose.structures(ix).(field);
end

function pass = quality_pass_detail(candidate, baseline, tolerance)
pass = true;
for i = 1:numel(candidate.structures)
    current = candidate.structures(i);
    reference = baseline.structures(i);
    if current.isTarget
        pass = pass && percent_change(current.D95, reference.D95) >= -tolerance && ...
            percent_change(current.D98, reference.D98) >= -tolerance;
    elseif strcmpi(current.type, 'OAR')
        pass = pass && percent_change(current.Dmean, reference.Dmean) <= tolerance && ...
            percent_change(current.Dmax, reference.Dmax) <= tolerance;
    end
end
end

function tbl = build_structure_comparison(dense, initial, final, tolerance)
n = numel(dense.structures);
rows = repmat(struct('structure', '', 'type', '', 'isTarget', false, ...
    'denseD98', NaN, 'initialD98', NaN, 'finalD98', NaN, ...
    'denseD95', NaN, 'initialD95', NaN, 'finalD95', NaN, ...
    'denseDmean', NaN, 'initialDmean', NaN, 'finalDmean', NaN, ...
    'denseDmax', NaN, 'initialDmax', NaN, 'finalDmax', NaN, ...
    'finalD95ChangePercent', NaN, 'finalD98ChangePercent', NaN, ...
    'finalDmeanChangePercent', NaN, 'finalDmaxChangePercent', NaN, ...
    'structurePass', true), n, 1);
for i = 1:n
    d = dense.structures(i); a = initial.structures(i); f = final.structures(i);
    rows(i).structure = d.name; rows(i).type = d.type; rows(i).isTarget = d.isTarget;
    rows(i).denseD98 = d.D98; rows(i).initialD98 = a.D98; rows(i).finalD98 = f.D98;
    rows(i).denseD95 = d.D95; rows(i).initialD95 = a.D95; rows(i).finalD95 = f.D95;
    rows(i).denseDmean = d.Dmean; rows(i).initialDmean = a.Dmean; rows(i).finalDmean = f.Dmean;
    rows(i).denseDmax = d.Dmax; rows(i).initialDmax = a.Dmax; rows(i).finalDmax = f.Dmax;
    rows(i).finalD95ChangePercent = percent_change(f.D95, d.D95);
    rows(i).finalD98ChangePercent = percent_change(f.D98, d.D98);
    rows(i).finalDmeanChangePercent = percent_change(f.Dmean, d.Dmean);
    rows(i).finalDmaxChangePercent = percent_change(f.Dmax, d.Dmax);
    if d.isTarget
        rows(i).structurePass = rows(i).finalD95ChangePercent >= -tolerance && ...
            rows(i).finalD98ChangePercent >= -tolerance;
    elseif strcmpi(d.type, 'OAR')
        rows(i).structurePass = rows(i).finalDmeanChangePercent <= tolerance && ...
            rows(i).finalDmaxChangePercent <= tolerance;
    end
end
tbl = struct2table(rows);
end

function tbl = build_stage_comparison(stageDose)
dense = stageDose.dense;
stages = {'dense','restrictedInitial','lambdaZeroOptimizer','finalIpopt'};
n = numel(dense.structures);
rows = repmat(struct('structure', '', 'type', '', 'metric', '', ...
    'dense', NaN, 'restrictedInitial', NaN, ...
    'lambdaZeroOptimizer', NaN, 'finalIpopt', NaN), 4*n, 1);
fields = {'D98','D95','Dmean','Dmax'};
k = 0;
for i = 1:n
    for j = 1:numel(fields)
        k = k + 1;
        rows(k).structure = dense.structures(i).name;
        rows(k).type = dense.structures(i).type;
        rows(k).metric = fields{j};
        for s = 1:numel(stages)
            rows(k).(stages{s}) = stageDose.(stages{s}).structures(i).(fields{j});
        end
    end
end
tbl = struct2table(rows);
end

function p = percent_change(value, baseline)
p = 100 * (double(value) - double(baseline)) / max(abs(double(baseline)), eps);
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
