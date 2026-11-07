function outputs = run_stage2_admm_rho_dynamics_audit(outputDir, auditFile)
%RUN_STAGE2_ADMM_RHO_DYNAMICS_AUDIT Reduced Variant-H dynamics only.
%
% Compares IPOPT-reference and PG-20 under fixed rho and bounded residual
% balancing. No full problem, Variant I, or other optimizer is reachable.

narginchk(2, 2);
root = fileparts(fileparts(mfilename('fullpath')));
add_module_paths(root);
create_empty_output_directory(outputDir);
proxTests = runtests(fullfile(root, 'tests', 'test_prox_l0_mmu.m'));
rhoTests = runtests(fullfile(root, 'tests', 'test_admm_rho_balance.m'));
assert(all([proxTests.Passed]) && all([rhoTests.Passed]), ...
    'ADMM operator/rho unit tests failed.');

artifactFile = fullfile(root, 'results', 'task25_27_pipeline_20deg.mat');
frozenFile = fullfile(root, 'results', 'full_integrated_layer_20deg_job22180', ...
    'alpha_0p3', 'alpha_0p3_final.mat');
expectedHash = '9c77e38f26f43515fc1515a29523b94ec2c52265df74e9686f0196efc7097f18';
verify_sha256(frozenFile, expectedHash);
source = load(artifactFile, 'cfg', 'machine', 'dij', 'pln', 'map');
layer = load(frozenFile, 'record', 'final');
units = load(auditFile, 'audit', 'planningMUPerSpot', 'nppm', 'machineMinMU');
patient = load(source.cfg.patientFile, 'cst');
assert(layer.record.planAccepted && abs(layer.record.alpha - 0.3) <= 1e-12);

cfg = merge_structs_local(default_config(), source.cfg);
cfg.lambdaB = 0; cfg.lambdaES = 0; cfg.lambdaT = 0; cfg.lambdaS = 0;
wStage1 = reshape(double(layer.final.w), [], 1);
allowedMask = wStage1 > 0;
assert(nnz(allowedMask) == 136409);
[candidateMask, selectionDiagnostics] = ...
    build_layer_stratified_support(wStage1, source.map, 100);
reducedMask = candidateMask & allowedMask;
[dijReduced, reducedToFull, fullToReduced, mapReduced] = ...
    reduce_dij_to_support(source.dij, reducedMask, source.map);
for i = 1:numel(mapReduced)
    mapReduced(i).originalGlobalIx = mapReduced(i).globalIx;
    mapReduced(i).globalIx = i;
end
[groupsReduced, ~] = build_group_index(mapReduced);
[groupsFull, ~] = build_group_index(source.map);
w0 = wStage1(reducedToFull);

planningMU = units.planningMUPerSpot(reducedToFull);
nppm = units.nppm(reducedToFull);
minimumWeight = nppm .* planningMU ./ 1e6;
assert(all(isfinite(minimumWeight) & minimumWeight >= 0));
assert(all(units.machineMinMU == 0), ...
    'This audit assumes the documented zero machine minMU artifact.');

doseModel = build_matrad_dose_model( ...
    dijReduced, patient.cst, source.pln, cfg.integrated.constraintPenalty);
positive = sort(w0(w0 > 0));
tau = positive(max(1, ceil(0.10 * numel(positive))));
rho0 = 1;
lambda0 = 0.5 * rho0 * tau^2;
problem = struct('doseModel', doseModel, 'w0', w0, ...
    'minimumWeight', minimumWeight, 'groups', groupsReduced, ...
    'reducedToFull', reducedToFull);
problem.bdtFunction = @(z) reduced_bdt(z, reducedToFull, ...
    numel(wStage1), source.map, groupsFull, cfg, source.machine);

balance = struct('mu', 1.5, 'increaseFactor', 2, 'decreaseFactor', 2, ...
    'minimum', 0.125, 'maximum', 8);
specs = { ...
    spec('ipopt_fixed', 'ipopt', 'fixed'), ...
    spec('ipopt_balanced', 'ipopt', 'residual_balancing'), ...
    spec('pg20_fixed', 'pg', 'fixed'), ...
    spec('pg20_balanced', 'pg', 'residual_balancing')};
records = repmat(empty_record(), numel(specs), 1);
setup = struct('classification', 'reduced Variant-H rho dynamics audit', ...
    'reducedSpots', numel(w0), 'representedAllowedLayers', numel(groupsReduced), ...
    'rho0', rho0, 'lambda0', lambda0, 'thresholdTauAtRho0', tau, ...
    'rhoBalance', balance, 'planningMinimumWeight', cfg.mmuThreshold, ...
    'machineMinimumMU', 0, 'maxAdmmIterations', 8, ...
    'fullProblemRun', false, 'variantIImplemented', false);
save(fullfile(outputDir, 'rho_audit_setup.mat'), 'setup', 'reducedMask', ...
    'reducedToFull', 'fullToReduced', 'selectionDiagnostics', '-v7.3');

for j = 1:numel(specs)
    s = specs{j};
    settings = struct('lambda0', lambda0, 'rho', rho0, ...
        'rhoStrategy', s.rhoStrategy, 'rhoBalance', balance, ...
        'maxIterations', 8, 'minimumIterations', 4, ...
        'absoluteTolerance', 1e-6, 'relativeTolerance', 1e-3, ...
        'wStrategy', s.wStrategy);
    if strcmp(s.wStrategy, 'pg')
        settings.inner = struct('maxIterations', 20, 'initialStep', 1e-3);
    else
        settings.ipopt = struct('ipopt', struct('max_iter', 300, ...
            'max_cpu_time', 3600, 'print_level', 3));
    end
    fprintf('\n=== Rho dynamics audit: %s ===\n', s.name);
    result = run_stage2_admm_variant_h(problem, settings);
    trajectory = struct2table(result.history);
    events = struct2table(result.supportEvents);
    writetable(trajectory, fullfile(outputDir, [s.name '_trajectory.csv']));
    writetable(events, fullfile(outputDir, [s.name '_support_events.csv']));
    cycleEvents = extract_cycle_events(result);
    writetable(cycleEvents, fullfile(outputDir, [s.name '_cycle_events.csv']));
    supportIdentity = result.z > 0; %#ok<NASGU>
    save(fullfile(outputDir, [s.name '_result.mat']), 'result', ...
        'supportIdentity', 'settings', 'cycleEvents', '-v7.3');
    records(j) = summarize(s, result, cycleEvents);
    writetable(struct2table(records(1:j)), ...
        fullfile(outputDir, 'rho_dynamics_summary.csv'));
end
save(fullfile(outputDir, 'rho_dynamics_complete.mat'), ...
    'setup', 'records', '-v7.3');
outputs = struct('setup', setup, 'records', records, 'outputDir', outputDir);
end

function s = spec(name, wStrategy, rhoStrategy)
s = struct('name', name, 'wStrategy', wStrategy, 'rhoStrategy', rhoStrategy);
end

function tableOut = extract_cycle_events(result)
emptyEvents = result.supportEvents;
emptyEvents = emptyEvents([]);
tableOut = struct2table(emptyEvents);
if isempty(result.history), return; end
p = result.history(end).supportCyclePeriod;
k = result.history(end).iteration;
if p <= 0 || isempty(result.supportEvents), return; end
events = struct2table(result.supportEvents);
tableOut = events(events.iteration >= k-p+1 & events.iteration <= k, :);
end

function r = summarize(specification, result, cycleEvents)
h = result.history; last = h(end); r = empty_record();
r.name = string(specification.name); r.wStrategy = string(specification.wStrategy);
r.rhoStrategy = string(specification.rhoStrategy); r.iterations = result.iterations;
r.stopReason = string(result.stopReason); r.finalRho = last.rhoNext;
r.rhoChanges = nnz([h.rho] ~= [h.rhoNext]);
r.primalResidual = last.primalResidual; r.dualResidual = last.dualResidual;
r.objective = last.trueObjectiveZ; r.exactZeroSpots = last.exactZeroSpots;
r.activeSpots = last.activeSpots; r.activeLayers = last.activeLayers;
r.provisionalBDT = last.provisionalBDT;
r.objectiveEvaluations = result.callbackCounts.objective;
r.gradientEvaluations = result.callbackCounts.gradient;
r.wallSeconds = result.wallSeconds; r.supportCyclePeriod = last.supportCyclePeriod;
r.supportEventRows = numel(result.supportEvents);
r.cycleEventRows = height(cycleEvents); r.cycleUniqueSpots = 0;
r.cycleMedianRelativeBoundaryDistance = NaN;
r.cycleMaximumRelativeBoundaryDistance = NaN;
r.cycleFractionWithin1eMinus6 = NaN;
if ~isempty(cycleEvents)
    r.cycleUniqueSpots = numel(unique(cycleEvents.reducedIndex));
    d = cycleEvents.relativeBoundaryDistance;
    r.cycleMedianRelativeBoundaryDistance = median(d);
    r.cycleMaximumRelativeBoundaryDistance = max(d);
    r.cycleFractionWithin1eMinus6 = mean(d <= 1e-6);
end
end

function r = empty_record()
r = struct('name', "", 'wStrategy', "", 'rhoStrategy', "", ...
    'iterations', NaN, 'stopReason', "", 'finalRho', NaN, ...
    'rhoChanges', NaN, 'primalResidual', NaN, 'dualResidual', NaN, ...
    'objective', NaN, 'exactZeroSpots', NaN, 'activeSpots', NaN, ...
    'activeLayers', NaN, 'provisionalBDT', NaN, ...
    'objectiveEvaluations', NaN, 'gradientEvaluations', NaN, ...
    'wallSeconds', NaN, 'supportCyclePeriod', NaN, ...
    'supportEventRows', NaN, 'cycleEventRows', NaN, ...
    'cycleUniqueSpots', NaN, 'cycleMedianRelativeBoundaryDistance', NaN, ...
    'cycleMaximumRelativeBoundaryDistance', NaN, ...
    'cycleFractionWithin1eMinus6', NaN);
end

function bdt = reduced_bdt(z, reducedToFull, nFull, map, groups, cfg, machine)
wFull = zeros(nFull, 1); wFull(reducedToFull) = z;
d = compute_delivery_metrics(wFull, map, groups, cfg, machine);
bdt = d.T_total;
end

function add_module_paths(root)
folders = {'config','planning','mapping','objective','delivery','reduction', ...
    'integrated','integrated/admm','evaluation','tests'};
for i = 1:numel(folders), addpath(fullfile(root, folders{i})); end
end

function verify_sha256(file, expected)
[status, output] = system(sprintf('sha256sum "%s"', file)); assert(status == 0);
assert(strcmpi(strtok(strtrim(output)), expected), 'Frozen checksum mismatch.');
end

function merged = merge_structs_local(base, override)
merged = base; names = fieldnames(override);
for i = 1:numel(names)
    key = names{i};
    if isfield(merged,key) && isstruct(merged.(key)) && isstruct(override.(key))
        merged.(key) = merge_structs_local(merged.(key), override.(key));
    else, merged.(key) = override.(key); end
end
end

function create_empty_output_directory(path)
if exist(path, 'dir') == 7
    e = dir(path); e = e(~ismember({e.name}, {'.','..'}));
    assert(isempty(e), 'Output directory must be empty: %s', path);
else, mkdir(path); end
end
