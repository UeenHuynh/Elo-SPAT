function outputs = run_convergence_sensitivity_analysis_20deg(outputDir)
% Analyze stopping-rule sensitivity from existing iterations 200--600.
% This runner performs no dose calculation and no optimization.

narginchk(0, 1);
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'delivery'));
addpath(fullfile(root, 'evaluation'));
if nargin < 1 || isempty(outputDir)
    outputDir = fullfile(root, 'results', ...
        ['convergence_sensitivity_20deg_' datestr(now, 'yyyymmdd_HHMMSS')]);
end
create_empty_output_directory(outputDir);

artifactFile = fullfile(root, 'results', 'task25_27_pipeline_20deg.mat');
supportFile = fullfile(root, 'results', ...
    'support_representativeness_20deg_job21861', ...
    'support_86400_checkpoint.mat');
job200 = fullfile(root, 'results', ...
    'revised_convergence_validation_20deg_job21894');
job600 = fullfile(root, 'results', ...
    'convergence_continuation_20deg_job22141');
source = load(artifactFile, 'cfg', 'machine', 'map', 'groups');
support = load(supportFile, 'supportMask');
reducedToFull = find(reshape(logical(support.supportMask), [], 1));
nFull = numel(source.map);
fullToReduced = zeros(nFull, 1);
fullToReduced(reducedToFull) = (1:numel(reducedToFull))';

alphas = [0.01; 1]; tokens = {'0p01'; '1'};
iterations = [200, 300, 400, 500, 600];
windows = [10, 25, 50, 100];
checkpointRows = repmat(empty_checkpoint_row(), 0, 1);
windowRows = repmat(empty_window_row(), 0, 1);
endpointRows = repmat(empty_endpoint_row(), 0, 1);
doseRows = repmat(empty_dose_row(), 0, 1);
identityRows = repmat(empty_identity_row(), 0, 1);

for a = 1:numel(alphas)
    token = tokens{a};
    seed = load(fullfile(job200, sprintf('alpha_%s_checkpoint.mat', token)), ...
        'optimizer', 'stats', 'supportFull', 'finalDelivery');
    final = load(fullfile(job600, sprintf('alpha_%s_final.mat', token)), ...
        'optimizer', 'stats', 'supportFull', 'finalDelivery');
    states = cell(numel(iterations), 1); states{1} = seed.optimizer;
    for j = 2:numel(iterations)-1
        tmp = load(fullfile(job600, sprintf( ...
            'alpha_%s_iter_%04d_checkpoint.mat', token, iterations(j))), ...
            'optimizer');
        states{j} = tmp.optimizer;
    end
    states{end} = final.optimizer;

    priorBDT = NaN; priorZeroMask = [];
    for j = 1:numel(iterations)
        state = states{j};
        wFull = zeros(nFull, 1); wFull(reducedToFull) = state.w;
        delivery = compute_delivery_metrics( ...
            wFull, source.map, source.groups, source.cfg, source.machine);
        row = empty_checkpoint_row();
        row.alpha = alphas(a); row.iteration = iterations(j);
        row.objective = state.objective;
        row.normalizedStationarity = state.finalNormalizedStationarity;
        row.absoluteProxGradientResidual = ...
            state.history.proxGradientResidualNorm(end);
        row.stepSize = state.history.stepSize(end);
        row.exactZeroLayers = state.history.exactZeroLayers(end);
        row.exactActiveLayers = state.history.exactActiveLayers(end);
        row.activeSpots = state.history.activeSpots(end);
        row.integratedBDT = delivery.T_total;
        row.integratedBDTChangeFromPreviousPercent = ...
            pct(delivery.T_total, priorBDT);
        checkpointRows(end+1, 1) = row; %#ok<AGROW>
        zeroMask = exact_zero_mask(state.w, source.groups, fullToReduced);
        identity = empty_identity_row();
        identity.alpha = alphas(a); identity.iteration = iterations(j);
        identity.exactZeroLayers = nnz(zeroMask);
        if isempty(priorZeroMask)
            identity.sameAsPreviousCheckpoint = NaN;
        else
            identity.sameAsPreviousCheckpoint = isequal(zeroMask, priorZeroMask);
            identity.newZeroLayers = nnz(~priorZeroMask & zeroMask);
            identity.reactivatedLayers = nnz(priorZeroMask & ~zeroMask);
            identity.symmetricDifference = nnz(xor(priorZeroMask, zeroMask));
        end
        identityRows(end+1, 1) = identity; %#ok<AGROW>
        priorZeroMask = zeroMask;
        priorBDT = delivery.T_total;
    end

    h = final.optimizer.history;
    for endpoint = [200, 600]
        q = find(h.iteration == endpoint, 1);
        for window = windows
            ix = (q-window+1):q;
            row = empty_window_row();
            row.alpha = alphas(a); row.endpoint = endpoint; row.window = window;
            row.normalizedStationarityEnd = h.normalizedStationarity(q);
            row.normalizedStationarityMax = max(h.normalizedStationarity(ix));
            row.relativeObjectiveWindow = ...
                (h.objective(ix(1))-h.objective(ix(end))) / ...
                max(1, abs(h.objective(ix(1))));
            row.maxPerIterationRelativeObjective = max(h.relativeObjective(ix));
            row.exactZeroLayerCountRange = ...
                max(h.exactZeroLayers(ix))-min(h.exactZeroLayers(ix));
            row.layerCountStable = row.exactZeroLayerCountRange == 0;
            row.stepSizeRange = max(h.stepSize(ix))-min(h.stepSize(ix));
            row.stepSizeStable = row.stepSizeRange <= 1e-12;
            windowRows(end+1, 1) = row; %#ok<AGROW>
        end
    end

    e = empty_endpoint_row();
    e.alpha = alphas(a);
    e.zeroLayers200 = nnz(seed.stats.exactZeroMask);
    e.zeroLayers600 = nnz(final.stats.exactZeroMask);
    e.zeroLayersRetained = nnz(seed.stats.exactZeroMask & final.stats.exactZeroMask);
    e.zeroLayersNew = nnz(~seed.stats.exactZeroMask & final.stats.exactZeroMask);
    e.zeroLayersReactivated = nnz(seed.stats.exactZeroMask & ~final.stats.exactZeroMask);
    e.supportSpots200 = nnz(seed.supportFull); e.supportSpots600 = nnz(final.supportFull);
    e.supportRemoved = nnz(seed.supportFull & ~final.supportFull);
    e.supportAdded = nnz(~seed.supportFull & final.supportFull);
    e.supportJaccard = nnz(seed.supportFull & final.supportFull) / ...
        nnz(seed.supportFull | final.supportFull);
    e.finalBDT200 = seed.finalDelivery.T_total;
    e.finalBDT600 = final.finalDelivery.T_total;
    e.finalBDTChangePercent = pct(e.finalBDT600, e.finalBDT200);
    endpointRows(end+1, 1) = e; %#ok<AGROW>

    dose200 = readtable(fullfile(job200, sprintf( ...
        'alpha_%s_final_ipopt_dose.csv', token)));
    dose600 = readtable(fullfile(job600, sprintf( ...
        'alpha_%s_final_ipopt_dose.csv', token)));
    assert(isequal(dose200.structure, dose600.structure));
    for s = 1:height(dose200)
        d = empty_dose_row(); d.alpha = alphas(a);
        d.structure = dose200.structure{s}; d.type = dose200.type{s};
        d.D98ChangeGy = dose600.finalD98(s)-dose200.finalD98(s);
        d.D95ChangeGy = dose600.finalD95(s)-dose200.finalD95(s);
        d.DmeanChangeGy = dose600.finalDmean(s)-dose200.finalDmean(s);
        d.DmaxChangeGy = dose600.finalDmax(s)-dose200.finalDmax(s);
        d.adverseDensePercent200 = dose200.adverseChangePercent(s);
        d.adverseDensePercent600 = dose600.adverseChangePercent(s);
        d.pass200 = dose200.structurePass(s); d.pass600 = dose600.structurePass(s);
        doseRows(end+1, 1) = d; %#ok<AGROW>
    end
end

checkpointTable = struct2table(checkpointRows);
windowTable = struct2table(windowRows);
endpointTable = struct2table(endpointRows);
doseTable = struct2table(doseRows);
identityTable = struct2table(identityRows);
ruleTable = candidate_rule_table(windowTable, endpointTable, job200, job600);
writetable(checkpointTable, fullfile(outputDir, 'checkpoint_bdt_trajectory.csv'));
writetable(windowTable, fullfile(outputDir, 'trailing_window_metrics.csv'));
writetable(endpointTable, fullfile(outputDir, 'endpoint_support_bdt_impact.csv'));
writetable(doseTable, fullfile(outputDir, 'endpoint_dose_impact.csv'));
writetable(identityTable, fullfile(outputDir, 'checkpoint_layer_identity.csv'));
writetable(ruleTable, fullfile(outputDir, 'candidate_stopping_rules.csv'));
save(fullfile(outputDir, 'convergence_sensitivity_complete.mat'), ...
    'checkpointTable', 'windowTable', 'endpointTable', 'doseTable', ...
    'identityTable', 'ruleTable', 'artifactFile', 'supportFile', ...
    'job200', 'job600', '-v7.3');
outputs = struct('checkpointTable', checkpointTable, 'windowTable', windowTable, ...
    'endpointTable', endpointTable, 'doseTable', doseTable, ...
    'identityTable', identityTable, 'ruleTable', ruleTable, ...
    'outputDir', outputDir);
end

function zeroMask = exact_zero_mask(wReduced, groupsFull, fullToReduced)
zeroMask = false(numel(groupsFull), 1);
for g = 1:numel(groupsFull)
    ix = fullToReduced(groupsFull(g).ix);
    ix = ix(ix > 0);
    if isempty(ix)
        error('Representative support unexpectedly omits full group %d.', g);
    end
    zeroMask(g) = all(wReduced(ix) == 0);
end
end

function tbl = candidate_rule_table(windows, endpoints, job200, job600)
stationarityGrid = [1e-7, 1e-6, 2e-6, 5e-6, 1e-5];
objectiveGrid = [1e-5, 5e-5, 1e-4, 2e-4];
window = 50;
rows = repmat(struct('alpha',NaN,'endpoint',NaN,'window',window, ...
    'stationarityTolerance',NaN,'relativeObjectiveWindowTolerance',NaN, ...
    'stationarityPass',false,'objectivePass',false,'supportStablePass',false, ...
    'stepStablePass',false,'numericalPass',false,'doseQualityPass',false, ...
    'supportPreservedPass',false,'BDTChange200To600Percent',NaN), 0, 1);
for a = [0.01, 1]
    token = alpha_token(a);
    s200 = readtable(fullfile(job200, 'revised_convergence_summary.csv'));
    s600 = readtable(fullfile(job600, 'continuation_summary.csv'));
    r200 = s200(abs(s200.alpha-a)<1e-12, :);
    r600 = s600(abs(s600.alpha-a)<1e-12, :);
    ep = endpoints(abs(endpoints.alpha-a)<1e-12, :);
    for endpoint = [200, 600]
        w = windows(abs(windows.alpha-a)<1e-12 & ...
            windows.endpoint==endpoint & windows.window==window, :);
        if endpoint == 200, clinical = r200; else, clinical = r600; end
        for st = stationarityGrid
            for ot = objectiveGrid
                row = rows_template(window);
                row.alpha = a; row.endpoint = endpoint;
                row.stationarityTolerance = st;
                row.relativeObjectiveWindowTolerance = ot;
                row.stationarityPass = w.normalizedStationarityMax <= st;
                row.objectivePass = w.relativeObjectiveWindow <= ot;
                row.supportStablePass = w.layerCountStable;
                row.stepStablePass = w.stepSizeStable;
                row.numericalPass = row.stationarityPass && row.objectivePass && ...
                    row.supportStablePass && row.stepStablePass;
                row.doseQualityPass = clinical.doseQualityPass;
                row.supportPreservedPass = clinical.supportPreserved;
                row.BDTChange200To600Percent = ep.finalBDTChangePercent;
                rows(end+1,1) = row; %#ok<AGROW>
            end
        end
    end
end
tbl = struct2table(rows);
end

function r = rows_template(window)
r = struct('alpha',NaN,'endpoint',NaN,'window',window, ...
    'stationarityTolerance',NaN,'relativeObjectiveWindowTolerance',NaN, ...
    'stationarityPass',false,'objectivePass',false,'supportStablePass',false, ...
    'stepStablePass',false,'numericalPass',false,'doseQualityPass',false, ...
    'supportPreservedPass',false,'BDTChange200To600Percent',NaN);
end

function r = empty_checkpoint_row()
r = struct('alpha',NaN,'iteration',NaN,'objective',NaN, ...
    'normalizedStationarity',NaN,'absoluteProxGradientResidual',NaN, ...
    'stepSize',NaN,'exactZeroLayers',NaN,'exactActiveLayers',NaN, ...
    'activeSpots',NaN,'integratedBDT',NaN, ...
    'integratedBDTChangeFromPreviousPercent',NaN);
end

function r = empty_window_row()
r = struct('alpha',NaN,'endpoint',NaN,'window',NaN, ...
    'normalizedStationarityEnd',NaN,'normalizedStationarityMax',NaN, ...
    'relativeObjectiveWindow',NaN,'maxPerIterationRelativeObjective',NaN, ...
    'exactZeroLayerCountRange',NaN,'layerCountStable',false, ...
    'stepSizeRange',NaN,'stepSizeStable',false);
end

function r = empty_endpoint_row()
r = struct('alpha',NaN,'zeroLayers200',NaN,'zeroLayers600',NaN, ...
    'zeroLayersRetained',NaN,'zeroLayersNew',NaN,'zeroLayersReactivated',NaN, ...
    'supportSpots200',NaN,'supportSpots600',NaN,'supportRemoved',NaN, ...
    'supportAdded',NaN,'supportJaccard',NaN,'finalBDT200',NaN, ...
    'finalBDT600',NaN,'finalBDTChangePercent',NaN);
end

function r = empty_dose_row()
r = struct('alpha',NaN,'structure','','type','','D98ChangeGy',NaN, ...
    'D95ChangeGy',NaN,'DmeanChangeGy',NaN,'DmaxChangeGy',NaN, ...
    'adverseDensePercent200',NaN,'adverseDensePercent600',NaN, ...
    'pass200',false,'pass600',false);
end

function r = empty_identity_row()
r = struct('alpha',NaN,'iteration',NaN,'exactZeroLayers',NaN, ...
    'sameAsPreviousCheckpoint',NaN,'newZeroLayers',NaN, ...
    'reactivatedLayers',NaN,'symmetricDifference',NaN);
end

function create_empty_output_directory(outputDir)
if exist(outputDir, 'dir') == 7
    entries = dir(outputDir); entries = entries(~ismember({entries.name},{'.','..'}));
    if ~isempty(entries), error('Output directory must be empty: %s', outputDir); end
else
    mkdir(outputDir);
end
end

function token = alpha_token(alpha)
token = strrep(sprintf('%.4g', alpha), '.', 'p');
end

function value = pct(x, ref)
if ~isfinite(ref), value = NaN; else, value = 100*(x-ref)/max(abs(ref),eps); end
end
