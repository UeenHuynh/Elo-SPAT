function outputs = run_final_stabilization_audit_20deg(outputDir)
% Retrospective Gate-C audit using existing artifacts only.
% No optimization, dose calculation, or final IPOPT is run here.

narginchk(0, 1);
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'integrated'));
if nargin < 1 || isempty(outputDir)
    outputDir = fullfile(root, 'results', ...
        ['final_stabilization_audit_20deg_' datestr(now,'yyyymmdd_HHMMSS')]);
end
if exist(outputDir, 'dir') ~= 7, mkdir(outputDir); end
entries = dir(outputDir);
entries = entries(~ismember({entries.name},{'.','..'}));
if ~isempty(entries)
    error('Output directory must be empty: %s', outputDir);
end

rule = default_stabilization_rule();
oldDir = fullfile(root, 'results', ...
    'integrated_layer_validation_20deg_job21835');
continuationDir = fullfile(root, 'results', ...
    'convergence_continuation_20deg_job22141');
sensitivityDir = fullfile(root, 'results', ...
    'convergence_sensitivity_20deg_job22160');
bdtTable = readtable(fullfile(sensitivityDir, 'checkpoint_bdt_trajectory.csv'));
identityTable = readtable(fullfile(sensitivityDir, 'checkpoint_layer_identity.csv'));
clinical = readtable(fullfile(continuationDir, 'continuation_summary.csv'));

rows = repmat(empty_row(), 3, 1);
zeroTrajectory = readtable(fullfile(oldDir, ...
    'zero_lambda_objective_trajectory.csv'));
zeroSummary = readtable(fullfile(oldDir, 'zero_lambda_summary.csv'));
zeroHistory = old_zero_history(zeroTrajectory);
zeroReport = evaluate_numerical_stabilization(zeroHistory, struct([]), rule);
zeroPlan = evaluate_plan_acceptance( ...
    logical(zeroSummary.representativenessPassed(1)), ...
    zeroSummary.maxConstraintViolation(1), ...
    logical(zeroSummary.supportEmbeddingPreserved(1)), ...
    contains(zeroSummary.finalIpoptStatus{1}, 'Solved'), rule);
rows(1) = make_row('zero_control_job21835', 0, zeroReport, zeroPlan, ...
    'none: historical run has only 80 updates and no BDT history', false);

alphas = [0.01, 1];
for i = 1:numel(alphas)
    alpha = alphas(i); token = strrep(sprintf('%.4g',alpha),'.','p');
    trajectory = readtable(fullfile(continuationDir, ...
        sprintf('alpha_%s_trajectory.csv', token)));
    history = nonzero_history(trajectory);
    b = bdtTable(abs(bdtTable.alpha-alpha)<1e-12, :);
    bdt = struct('iteration',b.iteration,'value',b.integratedBDT);
    ids = identityTable(abs(identityTable.alpha-alpha)<1e-12 & ...
        identityTable.iteration>=500, :);
    identityVerified = height(ids) >= 2 && ...
        all(ids.sameAsPreviousCheckpoint(2:end) == 1) && ...
        all(ids.symmetricDifference(2:end) == 0);
    report = evaluate_numerical_stabilization(history, bdt, rule);
    report.supportStable = report.supportStable && identityVerified;
    report.numerically_stabilized = report.numerically_stabilized && ...
        identityVerified;
    if ~identityVerified
        report.reason = 'layer_support_identity_not_verified';
    end
    c = clinical(abs(clinical.alpha-alpha)<1e-12, :);
    plan = evaluate_plan_acceptance(logical(c.doseQualityPass), ...
        c.maxConstraintViolation, logical(c.supportPreserved), ...
        logical(c.finalIpoptAcceptable), rule);
    rows(i+1) = make_row(sprintf('alpha_%g_iter600',alpha), alpha, ...
        report, plan, ...
        'per-iteration counts plus exact masks at iterations 500 and 600', ...
        identityVerified);
end

resultTable = struct2table(rows, 'AsArray', true);
writetable(resultTable, fullfile(outputDir, 'stabilization_audit.csv'));
save(fullfile(outputDir, 'stabilization_audit_complete.mat'), ...
    'rule','resultTable','oldDir','continuationDir','sensitivityDir','-v7.3');
outputs = struct('rule',rule,'resultTable',resultTable,'outputDir',outputDir);
disp(resultTable(:, {'caseName','numerically_stabilized', ...
    'stationary_converged','plan_accepted','reason'}));
end

function h = old_zero_history(t)
n = height(t)-1; objective = t.objective(2:end);
h = struct('iteration',t.iteration(2:end), 'objective',objective, ...
    'relativeObjective',abs(diff(t.objective))./max(1,abs(t.objective(1:end-1))), ...
    'stepSize',t.stepSize(2:end), ...
    'normalizedStationarity',NaN(n,1), ...
    'layerSupportChanged',true(n,1)); % identity was not recorded
end

function h = nonzero_history(t)
t = t(t.iteration>0,:);
h = struct('iteration',t.iteration,'objective',t.objective, ...
    'relativeObjective',t.relativeObjective,'stepSize',t.stepSize, ...
    'normalizedStationarity',t.normalizedStationarity, ...
    'layerSupportChanged',[true; diff(t.exactZeroLayers)~=0]);
end

function row = make_row(name, alpha, r, p, evidence, identityVerified)
row = empty_row(); row.caseName = name; row.alpha = alpha;
row.iterations = r.iterations; row.objectiveWindowChange = r.relativeObjectiveWindow;
row.maxPerIterationObjectiveChange = r.maxPerIterationRelativeObjective;
row.objectiveStable = r.objectiveStable; row.supportStable = r.supportStable;
row.supportIdentityVerified = identityVerified; row.stepRatio = r.stepRatio;
row.stepStable = r.stepStable; row.bdtRelativeRange = r.bdtRelativeRange;
row.bdtStable = r.bdtStable; row.numerically_stabilized = r.numerically_stabilized;
row.normalizedStationarity = r.finalNormalizedStationarity;
row.stationary_converged = r.stationary_converged;
row.doseQualityPass = p.doseQualityPass; row.constraintPass = p.constraintPass;
row.supportPreserved = p.supportPreserved;
row.finalIpoptAcceptable = p.finalIpoptAcceptable;
row.plan_accepted = p.plan_accepted; row.reason = r.reason;
row.supportEvidence = evidence;
end

function row = empty_row()
row = struct('caseName','','alpha',NaN,'iterations',NaN, ...
    'objectiveWindowChange',NaN,'maxPerIterationObjectiveChange',NaN, ...
    'objectiveStable',false,'supportStable',false, ...
    'supportIdentityVerified',false,'stepRatio',NaN,'stepStable',false, ...
    'bdtRelativeRange',NaN,'bdtStable',false, ...
    'numerically_stabilized',false,'normalizedStationarity',NaN, ...
    'stationary_converged',false,'doseQualityPass',false, ...
    'constraintPass',false,'supportPreserved',false, ...
    'finalIpoptAcceptable',false,'plan_accepted',false,'reason','', ...
    'supportEvidence','');
end
