function outputs = run_stage2_admm_variant_h_reduced(outputDir, auditFile)
%RUN_STAGE2_ADMM_VARIANT_H_REDUCED Reduced-only solver benchmark.
%
% No full 136,409-variable ADMM solve and no Variant I code path exists here.

narginchk(2, 2);
root = fileparts(fileparts(mfilename('fullpath')));
add_module_paths(root);
create_empty_output_directory(outputDir);

testResult = runtests(fullfile(root, 'tests', 'test_prox_l0_mmu.m'));
assert(all([testResult.Passed]), 'Exact L0+MMU operator unit tests failed.');

artifactFile = fullfile(root, 'results', 'task25_27_pipeline_20deg.mat');
frozenFile = fullfile(root, 'results', 'full_integrated_layer_20deg_job22180', ...
    'alpha_0p3', 'alpha_0p3_final.mat');
expectedHash = '9c77e38f26f43515fc1515a29523b94ec2c52265df74e9686f0196efc7097f18';
verify_sha256(frozenFile, expectedHash);
source = load(artifactFile, 'cfg', 'machine', 'dij', 'pln', 'map');
layer = load(frozenFile, 'record', 'final');
units = load(auditFile, 'audit', 'planningMUPerSpot', 'nppm', 'machineMinMU');
assert(layer.record.planAccepted && abs(layer.record.alpha - 0.3) <= 1e-12);
assert(units.audit.nSpots == numel(source.map));

patient = load(source.cfg.patientFile, 'cst');
cfg = merge_structs_local(default_config(), source.cfg);
cfg.lambdaB = 0; cfg.lambdaES = 0; cfg.lambdaT = 0; cfg.lambdaS = 0;
wStage1 = reshape(double(layer.final.w), [], 1);
allowedMask = wStage1 > 0;
assert(nnz(allowedMask) == 136409);

% Small layer-stratified subset of the allowed Stage-1 set. The 20 excluded
% Stage-1 layers can never enter the reduced coordinates.
capPerCP = 100;
[candidateMask, selectionDiagnostics] = ...
    build_layer_stratified_support(wStage1, source.map, capPerCP);
reducedMask = candidateMask & allowedMask;
assert(~any(reducedMask & ~allowedMask));
[dijReduced, reducedToFull, fullToReduced, mapReduced] = ...
    reduce_dij_to_support(source.dij, reducedMask, source.map);
for i = 1:numel(mapReduced)
    mapReduced(i).originalGlobalIx = mapReduced(i).globalIx;
    mapReduced(i).globalIx = i;
end
[groupsReduced, ~] = build_group_index(mapReduced);
[groupsFull, ~] = build_group_index(source.map);
w0 = wStage1(reducedToFull);

% The project planning threshold is cfg.mmuThreshold in weight units. Its
% physical-MU image has now been explicitly computed by the audit. It is
% distinct from the machine lower bound unless the audit proves equality.
planningMU = units.planningMUPerSpot(reducedToFull);
nppm = units.nppm(reducedToFull);
minimumWeight = nppm .* planningMU ./ 1e6;
roundTripError = max(abs(planningMU - minimumWeight .* 1e6 ./ nppm));
assert(roundTripError <= 1e-14 * max(1, max(abs(planningMU))));
assert(all(isfinite(minimumWeight) & minimumWeight >= 0));

doseModel = build_matrad_dose_model( ...
    dijReduced, patient.cst, source.pln, cfg.integrated.constraintPenalty);
positive = sort(w0(w0 > 0));
tau = positive(max(1, ceil(0.10 * numel(positive))));
rho = 1;
lambda0 = 0.5 * rho * tau^2;

problem = struct('doseModel', doseModel, 'w0', w0, ...
    'minimumWeight', minimumWeight, 'groups', groupsReduced);
problem.bdtFunction = @(z) reduced_bdt(z, reducedToFull, ...
    numel(wStage1), source.map, groupsFull, cfg, source.machine);

strategies = { ...
    struct('name', 'ipopt_reference', 'wStrategy', 'ipopt', 'innerIterations', NaN), ...
    struct('name', 'pg_5', 'wStrategy', 'pg', 'innerIterations', 5), ...
    struct('name', 'pg_20', 'wStrategy', 'pg', 'innerIterations', 20), ...
    struct('name', 'lbfgs_5', 'wStrategy', 'lbfgs', 'innerIterations', 5), ...
    struct('name', 'lbfgs_20', 'wStrategy', 'lbfgs', 'innerIterations', 20)};
records = repmat(empty_record(), numel(strategies), 1);
setup = struct('classification', 'sequential reduced Variant-H prototype', ...
    'frozenStage1Hash', expectedHash, 'allowedFullSpots', nnz(allowedMask), ...
    'reducedSpots', nnz(reducedMask), 'representedAllowedLayers', numel(groupsReduced), ...
    'capPerControlPoint', capPerCP, 'rho', rho, 'lambda0', lambda0, ...
    'thresholdTau', tau, 'planningWeightThreshold', cfg.mmuThreshold, ...
    'planningMMURange', [min(planningMU), max(planningMU)], ...
    'machineMMURange', [min(units.machineMinMU), max(units.machineMinMU)], ...
    'machineMMUVerified', units.audit.machineMMUVerified, ...
    'planningThresholdEqualsMachineMMU', ...
        units.audit.planningThresholdIsPhysicalMMU, ...
    'roundTripMUError', roundTripError, 'maxAdmmIterations', 8, ...
    'variantIImplemented', false, 'fullProblemRun', false);
save(fullfile(outputDir, 'variant_h_reduced_setup.mat'), 'setup', ...
    'reducedMask', 'reducedToFull', 'fullToReduced', ...
    'selectionDiagnostics', 'minimumWeight', '-v7.3');

for j = 1:numel(strategies)
    spec = strategies{j};
    settings = struct('lambda0', lambda0, 'rho', rho, ...
        'maxIterations', 8, 'minimumIterations', 4, ...
        'absoluteTolerance', 1e-6, 'relativeTolerance', 1e-3, ...
        'wStrategy', spec.wStrategy);
    if isfinite(spec.innerIterations)
        settings.inner = struct('maxIterations', spec.innerIterations, ...
            'initialStep', 1e-3);
    else
        settings.ipopt = struct('ipopt', struct('max_iter', 300, ...
            'max_cpu_time', 3600, 'print_level', 3));
    end
    fprintf('\n=== Variant H reduced strategy: %s ===\n', spec.name);
    result = run_stage2_admm_variant_h(problem, settings);
    historyTable = struct2table(result.history);
    writetable(historyTable, fullfile(outputDir, ...
        sprintf('%s_trajectory.csv', spec.name)));
    supportIdentity = result.z > 0; %#ok<NASGU>
    save(fullfile(outputDir, sprintf('%s_result.mat', spec.name)), ...
        'result', 'supportIdentity', 'settings', '-v7.3');
    records(j) = summarize(spec.name, result);
    writetable(struct2table(records(1:j)), ...
        fullfile(outputDir, 'variant_h_strategy_summary.csv'));
end

save(fullfile(outputDir, 'variant_h_reduced_complete.mat'), ...
    'setup', 'records', '-v7.3');
outputs = struct('setup', setup, 'records', records, 'outputDir', outputDir);
end

function bdt = reduced_bdt(z, reducedToFull, nFull, map, groups, cfg, machine)
wFull = zeros(nFull, 1);
wFull(reducedToFull) = z;
delivery = compute_delivery_metrics(wFull, map, groups, cfg, machine);
bdt = delivery.T_total;
end

function r = summarize(name, result)
last = result.history(end);
r = empty_record();
r.strategy = string(name); r.iterations = result.iterations;
r.stopReason = string(result.stopReason); r.wallSeconds = result.wallSeconds;
r.objectiveEvaluations = result.callbackCounts.objective;
r.gradientEvaluations = result.callbackCounts.gradient;
r.primalResidual = last.primalResidual; r.dualResidual = last.dualResidual;
r.trueObjectiveZ = last.trueObjectiveZ; r.exactZeroSpots = last.exactZeroSpots;
r.activeSpots = last.activeSpots; r.activeLayers = last.activeLayers;
r.mmuFeasible = result.mmuFeasible;
r.minimumPositiveWeight = last.minimumPositiveZ;
r.constraintViolationW = last.constraintViolationW;
r.constraintViolationZ = last.constraintViolationZ;
r.provisionalBDT = last.provisionalBDT;
r.supportCyclePeriod = last.supportCyclePeriod;
end

function r = empty_record()
r = struct('strategy', "", 'iterations', NaN, 'stopReason', "", ...
    'wallSeconds', NaN, 'objectiveEvaluations', NaN, ...
    'gradientEvaluations', NaN, 'primalResidual', NaN, ...
    'dualResidual', NaN, 'trueObjectiveZ', NaN, ...
    'exactZeroSpots', NaN, 'activeSpots', NaN, 'activeLayers', NaN, ...
    'mmuFeasible', false, 'minimumPositiveWeight', NaN, ...
    'constraintViolationW', NaN, 'constraintViolationZ', NaN, ...
    'provisionalBDT', NaN, 'supportCyclePeriod', NaN);
end

function add_module_paths(root)
folders = {'config','planning','mapping','objective','delivery','reduction', ...
    'integrated','integrated/admm','evaluation','tests'};
for i = 1:numel(folders), addpath(fullfile(root, folders{i})); end
end

function verify_sha256(file, expected)
[status, output] = system(sprintf('sha256sum "%s"', file));
assert(status == 0, 'Could not hash frozen Stage-1 artifact.');
assert(strcmpi(strtok(strtrim(output)), expected), ...
    'Frozen Stage-1 artifact checksum mismatch.');
end

function merged = merge_structs_local(base, override)
merged = base;
names = fieldnames(override);
for i = 1:numel(names)
    key = names{i};
    if isfield(merged, key) && isstruct(merged.(key)) && ...
            isstruct(override.(key))
        merged.(key) = merge_structs_local(merged.(key), override.(key));
    else
        merged.(key) = override.(key);
    end
end
end

function create_empty_output_directory(path)
if exist(path, 'dir') == 7
    entries = dir(path);
    entries = entries(~ismember({entries.name}, {'.', '..'}));
    if ~isempty(entries)
        error('Output directory must be new or empty: %s', path);
    end
end
if exist(path, 'dir') ~= 7, mkdir(path); end
end
