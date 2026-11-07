function outputs = run_stage2_pdasc_reduced_study(outputDir)
%RUN_STAGE2_PDASC_REDUCED_STUDY Published PDASC support comparator.
%
% Scope is deliberately limited to the frozen 1,780-spot Stage-1 candidate
% set. PDASC selects support with its published primary-target LS/L0 model;
% every reported candidate is then polished with the same support-constrained
% matRad/IPOPT path and evaluated with the project dose/BDT pipeline.

narginchk(1, 1);
root = fileparts(fileparts(mfilename('fullpath')));
add_module_paths(root); create_empty_output_directory(outputDir);

testResults = runtests(fullfile(root, 'tests', 'test_pdasc_l0_ls.m'));
assert(all([testResults.Passed]), 'PDASC unit tests failed.');

artifactFile = fullfile(root, 'results', 'task25_27_pipeline_20deg.mat');
stage1File = fullfile(root, 'results', 'full_integrated_layer_20deg_job22180', ...
    'alpha_0p3', 'alpha_0p3_final.mat');
spotScaleFile = fullfile(root, 'results', ...
    'full_integrated_spot_setup_20deg_job22187', 'spot_setup_complete.mat');
admmDir = fullfile(root, 'results', 'stage2_admm_rho_dynamics_job22281');
admmFile = fullfile(admmDir, 'ipopt_fixed_result.mat');
admmSetupFile = fullfile(admmDir, 'rho_audit_setup.mat');
stage1Hash = '9c77e38f26f43515fc1515a29523b94ec2c52265df74e9686f0196efc7097f18';
admmHash = '3de84a2358cf2af89e364c332350c99d7ee4db39c3d2033e37baa715d013733c';
verify_sha256(stage1File, stage1Hash); verify_sha256(admmFile, admmHash);

source = load(artifactFile, 'cfg','machine','denseResult','dij','pln','map');
stage1 = load(stage1File, 'record','final','finalDose','finalDelivery','cfgTrial');
spotScale = load(spotScaleFile, 'setup');
patient = load(source.cfg.patientFile, 'ct','cst');
assert(stage1.record.planAccepted && abs(stage1.record.alpha - 0.3) <= 1e-12);

cfg = merge_structs_local(default_config(), source.cfg);
cfg.mode = 'integrated'; cfg.integrated.componentStage = 'spot';
cfg.lambdaEL = stage1.record.lambdaEL; cfg.lambdaB = 0;
cfg.lambdaES = 0; cfg.lambdaT = 0; cfg.enforceSelectedMinimumWeight = true;
alphaS = 0.03;
ixScale = find(abs(spotScale.setup.alphaSGrid - alphaS) <= 1e-12, 1);
assert(~isempty(ixScale)); cfg.lambdaS = spotScale.setup.lambdaS(ixScale);
cfg = configure_integrated_stage(cfg);

wStage1 = reshape(double(stage1.final.w), [], 1);
allowedMask = wStage1 > 0;
[candidateMask, selectionDiagnostics] = ...
    build_layer_stratified_support(wStage1, source.map, 100);
reducedMask = candidateMask & allowedMask;
[dijReduced, reducedToFull, fullToReduced, mapReduced] = ...
    reduce_dij_to_support(source.dij, reducedMask, source.map);
for i = 1:numel(mapReduced)
    mapReduced(i).originalGlobalIx = mapReduced(i).globalIx;
    mapReduced(i).globalIx = i;
end
[groupsReduced, cpGroupsReduced] = build_group_index(mapReduced);
[groupsFull, ~] = build_group_index(source.map);
assert(numel(reducedToFull) == 1780 && numel(groupsReduced) == 748);

doseModel = build_matrad_dose_model(dijReduced, patient.cst, source.pln, ...
    cfg.integrated.constraintPenalty);
context = build_integrated_context(groupsReduced, cpGroupsReduced, ...
    mapReduced, source.machine, cfg);
pdascModel = build_pdasc_primary_target_model(dijReduced, patient.cst, source.pln);

setup = struct('classification', ...
    'reduced published-form PDASC support comparator', ...
    'stage1Hash', stage1Hash, 'admmIpoptFixedHash', admmHash, ...
    'reducedSpots', numel(reducedToFull), 'allowedLayersRepresented', ...
    numel(groupsReduced), 'primaryTarget', pdascModel.targetName, ...
    'prescriptionPerFraction', pdascModel.prescriptionPerFraction, ...
    'targetVoxels', numel(pdascModel.targetVoxels), ...
    'pdascUsesClinicalObjective', false, 'pdascUsesOars', false, ...
    'commonPolishUsesClinicalObjective', true, ...
    'alphaSForPgMm', alphaS, 'lambdaSForPgMm', cfg.lambdaS, ...
    'lambdaELForPgMm', cfg.lambdaEL, 'fullSizePdascRun', false);
save(fullfile(outputDir, 'pdasc_study_setup.mat'), 'setup','reducedMask', ...
    'reducedToFull','fullToReduced','selectionDiagnostics','-v7.3');

%% Published-form PDASC continuation.
pdascSettings = struct('gamma',1,'lambdaMinRatio',1e-10, ...
    'maxContinuationLevels',41,'maxInnerIterations',50, ...
    'normalizeColumns',true,'restrictedSolver','normal_equations', ...
    'targetActiveSpots',Inf,'verbose',true);
pdasc = run_pdasc_l0_ls(pdascModel.U, pdascModel.p0, pdascSettings);
pdascTrajectory = augment_pdasc_trajectory(pdasc, groupsReduced, ...
    reducedToFull, numel(wStage1), source.map, groupsFull, cfg, source.machine);
writetable(pdascTrajectory, fullfile(outputDir, 'pdasc_continuation.csv'));
save(fullfile(outputDir, 'pdasc_selector_result.mat'), ...
    'pdasc','pdascModel','pdascSettings','pdascTrajectory','-v7.3');

%% Current period-1 smooth-L0 PG/MM on the identical reduced coordinates.
cfgPg = cfg; cfgPg.integrated.maxIterations = 100;
cfgPg.integrated.minimumIterations = 10;
cfgPg.integrated.stationarityStoppingEnabled = false;
cfgPg.integrated.spotMMUpdatePeriod = 1; cfgPg.integrated.verbose = true;
pgState = []; bdtHistory = struct('iteration',zeros(0,1),'value',zeros(0,1));
stabilization = struct('numerically_stabilized',false,'reason','not_started');
rule = default_stabilization_rule(); pgTimer = tic;
while isempty(pgState) || (pgState.iterations < 600 && ...
        ~stabilization.numerically_stabilized)
    if isempty(pgState)
        pgState = run_integrated_optimizer(doseModel, context, ...
            wStage1(reducedToFull), cfgPg);
    else
        pgState = run_integrated_optimizer(doseModel, context, ...
            pgState.w, cfgPg, pgState);
    end
    pgDelivery = reduced_delivery(pgState.w, reducedToFull, numel(wStage1), ...
        source.map, groupsFull, cfg, source.machine);
    bdtHistory.iteration(end+1,1) = pgState.iterations;
    bdtHistory.value(end+1,1) = pgDelivery.T_total;
    stabilization = evaluate_numerical_stabilization( ...
        pgState.history, bdtHistory, rule);
    save(fullfile(outputDir, sprintf('pgmm_iter_%04d.mat',pgState.iterations)), ...
        'pgState','bdtHistory','stabilization','cfgPg','-v7.3');
    if any(strcmp(pgState.status, {'line_search_failed','numerical_stagnation'}))
        break;
    end
end
pgWall = toc(pgTimer);
pgExtracted = extract_integrated_support(pgState.w, ...
    wStage1(reducedToFull), context, cfgPg);
save(fullfile(outputDir, 'pgmm_selector_result.mat'), 'pgState', ...
    'pgExtracted','bdtHistory','stabilization','pgWall','cfgPg','-v7.3');

%% Frozen reduced ADMM reference support; no ADMM continuation is run.
admm = load(admmFile, 'result','supportIdentity');
admmSetup = load(admmSetupFile, 'reducedToFull');
assert(isequal(double(admmSetup.reducedToFull(:)), double(reducedToFull(:))), ...
    'Frozen ADMM reduced coordinates differ from the PDASC study.');
admmSupport = reshape(logical(admm.supportIdentity), [], 1);
admmSeed = max(reshape(double(admm.result.z),[],1), 0);

% Immutable clinical references used in every post-polish comparison.
Dfull = nominal_matrix(source.dij);
denseDose = compute_dvh_metrics(full(Dfull*reshape(double(source.denseResult.w),[],1)), ...
    patient.ct,patient.cst,cfg);
denseDelivery = compute_delivery_metrics(source.denseResult.w,source.map, ...
    groupsFull,cfg,source.machine);
referenceRows = table(string({'Dense Arc';'Stage1 alphaEL=0.3'}), ...
    [denseDose.summary.targetD95Min;stage1.finalDose.summary.targetD95Min], ...
    [denseDose.summary.targetD98Min;stage1.finalDose.summary.targetD98Min], ...
    [denseDose.summary.oarDmaxMax;stage1.finalDose.summary.oarDmaxMax], ...
    [denseDelivery.T_total;stage1.finalDelivery.T_total], ...
    'VariableNames',{'reference','targetD95','targetD98','oarDmax','discreteBDT'});
writetable(referenceRows,fullfile(outputDir,'clinical_references.csv'));

%% Predeclared PDASC support sizes plus same-coordinate PG/MM and ADMM.
nReduced = numel(reducedToFull);
targetCounts = unique([round(0.50*nReduced), round(0.70*nReduced), ...
    round(0.80*nReduced), nnz(admmSupport)]);
pdascLevels = nearest_unique_levels(pdascTrajectory.projectedActiveSpots, targetCounts);
methods = repmat(empty_method(), 0, 1);
for j = 1:numel(pdascLevels)
    k = pdascLevels(j); state = pdasc.states{k};
    label = sprintf('PDASC_k%02d_n%d', k, nnz(state.support));
    methods(end+1) = polish_and_score(label, state.support, state.w, ...
        pdasc.history(k).elapsedSeconds, 'PDASC published LS/L0', ...
        source, patient, cfg, doseModel, reducedToFull, groupsFull, ...
        wStage1, denseDose, stage1, denseDelivery, outputDir, ...
        pdasc.history(k).innerStop); %#ok<AGROW>
end
methods(end+1) = polish_and_score('PGMM_period1_reduced', ...
    pgExtracted.mask, pgExtracted.w, pgWall, 'current smooth-L0 PG/MM', ...
    source, patient, cfg, doseModel, reducedToFull, groupsFull, ...
    wStage1, denseDose, stage1, denseDelivery, outputDir, ...
    stabilization.reason); %#ok<AGROW>
methods(end+1) = polish_and_score('ADMM_H_ipopt_fixed_frozen', ...
    admmSupport, admmSeed, admm.result.wallSeconds, ...
    'frozen job 22281 ADMM-H IPOPT fixed-rho', source, patient, cfg, ...
    doseModel, reducedToFull, groupsFull, wStage1, denseDose, stage1, ...
    denseDelivery, outputDir, admm.result.stopReason); %#ok<AGROW>

summary = struct2table(methods, 'AsArray', true);
writetable(summary, fullfile(outputDir, 'common_polish_benchmark.csv'));
save(fullfile(outputDir, 'pdasc_study_complete.mat'), 'setup','pdasc', ...
    'pdascTrajectory','pgState','pgExtracted','stabilization','admmSupport', ...
    'methods','targetCounts','pdascLevels','-v7.3');
outputs = struct('setup',setup,'pdasc',pdasc,'methods',methods, ...
    'outputDir',outputDir);
end

function tbl = augment_pdasc_trajectory(pdasc, groups, reducedToFull, ...
        nFull, map, groupsFull, cfg, machine)
h = pdasc.history; n = numel(h);
rows = repmat(struct('level',NaN,'lambda',NaN,'threshold',NaN, ...
    'innerIterations',NaN,'innerStop','','rawActiveSpots',NaN, ...
    'projectedActiveSpots',NaN,'exactZeroSpots',NaN,'activeLayers',NaN, ...
    'newActiveSpots',NaN,'removedActiveSpots',NaN,'projectedLeastSquares',NaN, ...
    'projectedObjective',NaN,'supportStableWithinLevel',false, ...
    'supportHash','','provisionalBDT',NaN,'elapsedSeconds',NaN),n,1);
for k=1:n
    support=pdasc.states{k}.support; d=reduced_delivery(pdasc.states{k}.w, ...
        reducedToFull,nFull,map,groupsFull,cfg,machine);
    rows(k).level=h(k).level; rows(k).lambda=h(k).lambda;
    rows(k).threshold=h(k).threshold; rows(k).innerIterations=h(k).innerIterations;
    rows(k).innerStop=h(k).innerStop; rows(k).rawActiveSpots=h(k).rawActiveSpots;
    rows(k).projectedActiveSpots=h(k).projectedActiveSpots;
    rows(k).exactZeroSpots=numel(support)-nnz(support);
    rows(k).activeLayers=count_active_layers(support,groups);
    rows(k).newActiveSpots=h(k).newActiveSpots;
    rows(k).removedActiveSpots=h(k).removedActiveSpots;
    rows(k).projectedLeastSquares=h(k).projectedLeastSquares;
    rows(k).projectedObjective=h(k).projectedObjective;
    rows(k).supportStableWithinLevel=strcmp(h(k).innerStop,'active_set_stable');
    rows(k).supportHash=h(k).projectedSupportHash;
    rows(k).provisionalBDT=d.T_total; rows(k).elapsedSeconds=h(k).elapsedSeconds;
end
tbl=struct2table(rows);
end

function levels = nearest_unique_levels(counts, targets)
levels=zeros(numel(targets),1);
for i=1:numel(targets), [~,levels(i)]=min(abs(counts-targets(i))); end
levels=unique(levels,'stable');
levels=levels(isfinite(counts(levels)) & counts(levels)>0);
end

function row = polish_and_score(name, supportReduced, seedReduced, selectorWall, ...
        selector, source, patient, cfg, doseModel, reducedToFull, groupsFull, ...
        wStage1, denseDose, stage1, denseDelivery, outputDir, selectorStatus)
supportReduced=reshape(logical(supportReduced),[],1);
seedReduced=reshape(double(seedReduced),[],1);
assert(numel(supportReduced)==numel(reducedToFull));
supportFull=false(numel(wStage1),1); supportFull(reducedToFull(supportReduced))=true;
wInit=zeros(size(wStage1));
seed=max(seedReduced(supportReduced),cfg.mmuThreshold);
wInit(reducedToFull(supportReduced))=seed;
timer=tic;
final=run_final_reoptimization(source.dij,patient.cst,source.pln, ...
    wInit,cfg,supportFull,[]);
polishWall=toc(timer);
tol=1e-8*max(1,max(abs(final.w)));
supportPreserved=~any(abs(final.w(~supportFull))>tol);
[clinicalObjective,~,diagnostics]=evaluate_matrad_dose_model( ...
    doseModel,final.w(reducedToFull),false);
D=nominal_matrix(source.dij);
dose=compute_dvh_metrics(full(D*final.w),patient.ct,patient.cst,cfg);
delivery=compute_delivery_metrics(final.w,source.map,groupsFull,cfg,source.machine);
solver=compute_solver_metrics(final.resultGUI,struct(),final.optimizer);
doseTable=dose_comparison_table(dose,denseDose,stage1.finalDose);
writetable(doseTable,fullfile(outputDir,[name '_dose_metrics.csv']));
selectedWeights=final.w(supportFull);
row=empty_method(); row.method=name; row.selector=selector;
row.selectorStatus=selectorStatus;
row.selectorWallSeconds=selectorWall; row.polishWallSeconds=polishWall;
row.totalWallSeconds=selectorWall+polishWall;
row.selectedSpots=nnz(supportReduced); row.exactZeroReducedSpots=numel(supportReduced)-nnz(supportReduced);
row.selectedLayers=count_active_layers(supportReduced,build_reduced_groups(groupsFull,reducedToFull));
row.finalActiveSpots=delivery.nActiveSpots; row.finalActiveLayers=delivery.nActiveLayers;
row.supportHash=support_hash(supportReduced); row.supportPreserved=supportPreserved;
row.finalIpoptStatus=solver.statusMessage; row.finalIpoptConverged=solver.converged;
row.clinicalObjective=clinicalObjective; row.doseObjective=diagnostics.doseObjective;
row.maxConstraintViolation=diagnostics.maxConstraintViolation;
row.targetD95=dose.summary.targetD95Min; row.targetD98=dose.summary.targetD98Min;
row.oarDmax=dose.summary.oarDmaxMax; row.discreteBDT=delivery.T_total;
row.d95ChangeFromDensePercent=pct(row.targetD95,denseDose.summary.targetD95Min);
row.d98ChangeFromDensePercent=pct(row.targetD98,denseDose.summary.targetD98Min);
row.bdtChangeFromDensePercent=pct(row.discreteBDT,denseDelivery.T_total);
row.bdtChangeFromStage1Percent=pct(row.discreteBDT,stage1.finalDelivery.T_total);
row.minimumSelectedWeight=min(selectedWeights); row.subMmuSpots= ...
    nnz(selectedWeights<cfg.mmuThreshold-tol);
record=row; %#ok<NASGU>
save(fullfile(outputDir,[name '_final.mat']),'record','final','dose', ...
    'delivery','solver','supportFull','supportReduced','selector','-v7.3');
end

function groupsReduced = build_reduced_groups(groupsFull,reducedToFull)
fullToReduced=zeros(max(reducedToFull),1); fullToReduced(reducedToFull)=1:numel(reducedToFull);
groupsReduced=groupsFull;
keep=false(numel(groupsFull),1);
for g=1:numel(groupsFull)
    ix=groupsFull(g).ix; ix=ix(ix<=numel(fullToReduced));
    local=fullToReduced(ix); local=local(local>0);
    groupsReduced(g).ix=local; keep(g)=~isempty(local);
end
groupsReduced=groupsReduced(keep);
end

function n=count_active_layers(support,groups)
n=0; for g=1:numel(groups), n=n+any(support(groups(g).ix)); end
end
function delivery=reduced_delivery(w,reducedToFull,nFull,map,groups,cfg,machine)
wFull=zeros(nFull,1); wFull(reducedToFull)=max(reshape(double(w),[],1),0);
delivery=compute_delivery_metrics(wFull,map,groups,cfg,machine);
end
function tbl=dose_comparison_table(candidate,dense,stage1)
n=numel(candidate.structures); rows=repmat(struct('structure','','type','', ...
    'isTarget',false,'denseD98',NaN,'stage1D98',NaN,'finalD98',NaN, ...
    'denseD95',NaN,'stage1D95',NaN,'finalD95',NaN, ...
    'denseDmean',NaN,'stage1Dmean',NaN,'finalDmean',NaN, ...
    'denseDmax',NaN,'stage1Dmax',NaN,'finalDmax',NaN),n,1);
for i=1:n
    a=candidate.structures(i); d=dense.structures(i); s=stage1.structures(i);
    rows(i).structure=a.name; rows(i).type=a.type; rows(i).isTarget=a.isTarget;
    for key={'D98','D95','Dmean','Dmax'}
        q=key{1}; rows(i).(['dense' q])=d.(q);
        rows(i).(['stage1' q])=s.(q); rows(i).(['final' q])=a.(q);
    end
end
tbl=struct2table(rows);
end
function row=empty_method()
row=struct('method','','selector','','selectorStatus','', ...
    'selectorWallSeconds',NaN, ...
    'polishWallSeconds',NaN,'totalWallSeconds',NaN,'selectedSpots',NaN, ...
    'exactZeroReducedSpots',NaN,'selectedLayers',NaN,'finalActiveSpots',NaN, ...
    'finalActiveLayers',NaN,'supportHash','','supportPreserved',false, ...
    'finalIpoptStatus','','finalIpoptConverged',false,'clinicalObjective',NaN, ...
    'doseObjective',NaN,'maxConstraintViolation',NaN,'targetD95',NaN, ...
    'targetD98',NaN,'oarDmax',NaN,'discreteBDT',NaN, ...
    'd95ChangeFromDensePercent',NaN,'d98ChangeFromDensePercent',NaN, ...
    'bdtChangeFromDensePercent',NaN,'bdtChangeFromStage1Percent',NaN, ...
    'minimumSelectedWeight',NaN,'subMmuSpots',NaN);
end
function D=nominal_matrix(dij), if iscell(dij.physicalDose), D=dij.physicalDose{1}; else, D=dij.physicalDose; end, end
function value=pct(x,ref), value=100*(double(x)-double(ref))/max(abs(double(ref)),eps); end
function value=support_hash(mask), md=java.security.MessageDigest.getInstance('SHA-256'); md.update(uint8(mask(:))); bytes=typecast(md.digest(),'uint8'); value=lower(reshape(dec2hex(bytes,2).',1,[])); end
function verify_sha256(file,expected), [s,o]=system(sprintf('sha256sum "%s"',file)); assert(s==0); assert(strcmpi(strtok(strtrim(o)),expected),'Frozen artifact checksum mismatch.'); end
function add_module_paths(root), f={'config','planning','mapping','objective','delivery','reduction','robustness','evaluation','integrated','integrated/pdasc'}; for i=1:numel(f), addpath(fullfile(root,f{i})); end, end
function create_empty_output_directory(d), if exist(d,'dir')==7, e=dir(d); e=e(~ismember({e.name},{'.','..'})); assert(isempty(e),'Output directory must be empty: %s',d); else, mkdir(d); end, end
function merged=merge_structs_local(base,override), merged=base; n=fieldnames(override); for i=1:numel(n), k=n{i}; if isfield(merged,k)&&isstruct(merged.(k))&&isstruct(override.(k)), merged.(k)=merge_structs_local(merged.(k),override.(k)); else, merged.(k)=override.(k); end, end, end
