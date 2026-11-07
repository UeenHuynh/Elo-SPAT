function outputs=run_working_set_reduced_benchmark(outputDir,pgComparatorFile,checkpointDir)
%RUN_WORKING_SET_REDUCED_BENCHMARK Independent exact-matRad working-set study.
narginchk(2,3); here=fileparts(mfilename('fullpath'));
if nargin < 3, checkpointDir=''; end
root=fileparts(fileparts(here)); add_paths(root,here);
create_empty_output_directory(outputDir);
if isempty(checkpointDir)
    pricingTests=runtests(fullfile(here,'tests','test_working_set_pricing.m'));
    quadraticTests=runtests(fullfile(here,'tests','test_working_set_quadratic.m'));
    assert(all([pricingTests.Passed])&&all([quadraticTests.Passed]), ...
        'Working-set tests failed.');
else
    assert(exist(checkpointDir,'dir')==7,'Checkpoint directory is missing.');
end

artifactFile=fullfile(root,'results','task25_27_pipeline_20deg.mat');
stage1File=fullfile(root,'results','full_integrated_layer_20deg_job22180', ...
    'alpha_0p3','alpha_0p3_final.mat');
admmFile=fullfile(root,'results','stage2_admm_rho_dynamics_job22281', ...
    'ipopt_fixed_result.mat');
admmSetupFile=fullfile(root,'results','stage2_admm_rho_dynamics_job22281', ...
    'rho_audit_setup.mat');
verify_sha256(stage1File,'9c77e38f26f43515fc1515a29523b94ec2c52265df74e9686f0196efc7097f18');
verify_sha256(admmFile,'3de84a2358cf2af89e364c332350c99d7ee4db39c3d2033e37baa715d013733c');
assert(exist(pgComparatorFile,'file')==2,'PG/MM comparator file is missing.');

source=load(artifactFile,'cfg','machine','denseResult','dij','pln','map');
stage1=load(stage1File,'record','final','finalDose','finalDelivery');
patient=load(source.cfg.patientFile,'ct','cst');
cfg=merge_structs_local(default_config(),source.cfg);
cfg.lambdaB=0;cfg.lambdaES=0;cfg.lambdaT=0;cfg.enforceSelectedMinimumWeight=true;
wStage1=reshape(double(stage1.final.w),[],1);allowed=wStage1>0;
[candidate,selectionDiagnostics]=build_layer_stratified_support(wStage1,source.map,100);
reducedMask=candidate&allowed;
[dijReduced,reducedToFull,fullToReduced,mapReduced]= ...
    reduce_dij_to_support(source.dij,reducedMask,source.map);
for i=1:numel(mapReduced),mapReduced(i).originalGlobalIx=mapReduced(i).globalIx;mapReduced(i).globalIx=i;end
[groupsReduced,~]=build_group_index(mapReduced);[groupsFull,~]=build_group_index(source.map);
assert(numel(reducedToFull)==1780&&numel(groupsReduced)==748);
fullModel=build_matrad_dose_model(dijReduced,patient.cst,source.pln,0);
if isempty(checkpointDir)
    callbackValidation=validate_working_set_callbacks(dijReduced,patient.cst, ...
        source.pln,wStage1(reducedToFull), ...
        fullfile(outputDir,'callback_gradient_validation.csv'));
else
    oldSetup=load(fullfile(checkpointDir,'working_set_setup.mat'), ...
        'reducedToFull','callbackValidation');
    assert(isequal(double(oldSetup.reducedToFull(:)),double(reducedToFull(:))), ...
        'Checkpoint reduced/full mapping mismatch.');
    callbackValidation=oldSetup.callbackValidation;
    copyfile(fullfile(checkpointDir,'callback_gradient_validation.csv'),outputDir);
end
problem=struct('dij',dijReduced,'cst',{patient.cst},'pln',source.pln, ...
    'fullModel',fullModel,'groups',groupsReduced);
problem.bdtFunction=@(w) reduced_bdt(w,reducedToFull,numel(wStage1), ...
    source.map,groupsFull,cfg,source.machine);
inner=struct('ipopt',struct('max_iter',300,'max_cpu_time',3600, ...
    'print_level',3,'hessian_approximation','limited-memory'));
specs={struct('name','working_set_spot','pricingMode','spot', ...
        'spotBatchSize',100,'layerBatchSize',0), ...
    struct('name','working_set_layer','pricingMode','layer', ...
        'spotBatchSize',50,'layerBatchSize',20)};
ws=cell(numel(specs),1);
setup=struct('classification','independent reduced exact-matRad working-set study', ...
    'candidateSpots',numel(reducedToFull),'candidateLayers',numel(groupsReduced), ...
    'forbiddenStage1Layers',20,'pgComparatorFile',pgComparatorFile, ...
    'deliveryPenaltyEnabled',false,'admmContinued',false,'pdascModified',false);
save(fullfile(outputDir,'working_set_setup.mat'),'setup','reducedMask', ...
    'reducedToFull','fullToReduced','selectionDiagnostics', ...
    'callbackValidation','-v7.3');
for j=1:numel(specs)
    s=specs{j};
    if isempty(checkpointDir)
        settings=struct('pricingMode',s.pricingMode, ...
            'spotBatchSize',s.spotBatchSize,'layerBatchSize',s.layerBatchSize, ...
            'maxOuterIterations',30,'absoluteTolerance',1e-7, ...
            'relativeTolerance',1e-6,'removeZeroNonbasic',true, ...
            'zeroTolerance',1e-10,'innerSettings',inner,'verbose',true);
        ws{j}=run_matrad_working_set(problem,settings);
        trajectory=struct2table(ws{j}.history); %#ok<NASGU>
        writetable(trajectory,fullfile(outputDir,[s.name '_trajectory.csv']));
        result=ws{j}; %#ok<NASGU>
        save(fullfile(outputDir,[s.name '_selector.mat']),'result','settings','-v7.3');
    else
        saved=load(fullfile(checkpointDir,[s.name '_selector.mat']),'result','settings');
        ws{j}=saved.result;
        settings=saved.settings; %#ok<NASGU>
        copyfile(fullfile(checkpointDir,[s.name '_trajectory.csv']),outputDir);
        result=ws{j}; %#ok<NASGU>
        save(fullfile(outputDir,[s.name '_selector.mat']),'result','settings','-v7.3');
    end
end

pg=load(pgComparatorFile,'pgState','pgExtracted','bdtHistory','stabilization','pgWall');
assert(numel(pg.pgExtracted.mask)==numel(reducedToFull));
admm=load(admmFile,'result','supportIdentity');
admmSetup=load(admmSetupFile,'reducedToFull');
assert(isequal(double(admmSetup.reducedToFull(:)),double(reducedToFull(:))));

methods=repmat(empty_method(),0,1);
for j=1:numel(specs)
    support=ws{j}.w>1e-10;
    methods(end+1)=polish(specs{j}.name,support,ws{j}.w, ...
        ws{j}.wallSeconds,ws{j}.objectiveEvaluations,ws{j}.gradientEvaluations, ...
        'exact',ws{j}.stopReason,source,patient,cfg,fullModel,reducedToFull, ...
        groupsReduced,groupsFull,wStage1,stage1,outputDir); %#ok<AGROW>
    writetable(struct2table(methods,'AsArray',true), ...
        fullfile(outputDir,'working_set_benchmark.csv'));
end
[pgObj,pgGrad]=estimate_pg_work(pg.pgState,pg.bdtHistory);
methods(end+1)=polish('PGMM_period1_reduced',pg.pgExtracted.mask, ...
    pg.pgExtracted.w,pg.pgWall,pgObj,pgGrad,'estimated_from_call_structure', ...
    pg.stabilization.reason, ...
    source,patient,cfg,fullModel,reducedToFull,groupsReduced,groupsFull, ...
    wStage1,stage1,outputDir); %#ok<AGROW>
methods(end+1)=polish('ADMM_H_ipopt_fixed_frozen',admm.supportIdentity, ...
    max(admm.result.z,0),admm.result.wallSeconds, ...
    admm.result.callbackCounts.objective,admm.result.callbackCounts.gradient, ...
    'exact',admm.result.stopReason,source,patient,cfg,fullModel,reducedToFull, ...
    groupsReduced,groupsFull,wStage1,stage1,outputDir); %#ok<AGROW>
writetable(struct2table(methods,'AsArray',true), ...
    fullfile(outputDir,'working_set_benchmark.csv'));
save(fullfile(outputDir,'working_set_benchmark_complete.mat'), ...
    'setup','ws','methods','specs','-v7.3');
outputs=struct('setup',setup,'workingSet',{ws},'methods',methods,'outputDir',outputDir);
end

function row=polish(name,supportReduced,seedReduced,selectorWall,objEvals,gradEvals, ...
        evaluationCountType,selectorStatus,source,patient,cfg,fullModel,reducedToFull,groupsReduced, ...
        groupsFull,wStage1,stage1,outputDir)
supportReduced=reshape(logical(supportReduced),[],1);seedReduced=reshape(double(seedReduced),[],1);
supportFull=false(size(wStage1));supportFull(reducedToFull(supportReduced))=true;
wInit=zeros(size(wStage1));wInit(reducedToFull(supportReduced))= ...
    max(seedReduced(supportReduced),cfg.mmuThreshold);
timer=tic;final=run_final_reoptimization(source.dij,patient.cst,source.pln, ...
    wInit,cfg,supportFull,[]);polishWall=toc(timer);
[objective,~,diagnostics]=evaluate_matrad_dose_model(fullModel, ...
    final.w(reducedToFull),false);
D=nominal_matrix(source.dij);dose=compute_dvh_metrics(full(D*final.w), ...
    patient.ct,patient.cst,cfg);
delivery=compute_delivery_metrics(final.w,source.map,groupsFull,cfg,source.machine);
solver=compute_solver_metrics(final.resultGUI,struct(),final.optimizer);
tol=1e-8*max(1,max(abs(final.w)));
row=empty_method();row.method=name;row.selectorStatus=selectorStatus;
row.evaluationCountType=evaluationCountType;
row.selectorWallSeconds=selectorWall;row.polishWallSeconds=polishWall;
row.objectiveEvaluations=objEvals;row.gradientEvaluations=gradEvals;
row.selectedSpots=nnz(supportReduced);row.selectedLayers=count_layers(supportReduced,groupsReduced);
row.finalActiveSpots=delivery.nActiveSpots;row.finalActiveLayers=delivery.nActiveLayers;
row.clinicalObjective=objective;row.maxConstraintViolation=diagnostics.maxConstraintViolation;
row.targetD95=dose.summary.targetD95Min;row.targetD98=dose.summary.targetD98Min;
row.oarDmax=dose.summary.oarDmaxMax;row.discreteBDT=delivery.T_total;
row.supportPreserved=~any(abs(final.w(~supportFull))>tol);
row.finalIpoptStatus=solver.statusMessage;
tbl=dose_table(dose,stage1.finalDose);writetable(tbl,fullfile(outputDir,[name '_dose.csv']));
record=row; %#ok<NASGU>
save(fullfile(outputDir,[name '_final.mat']),'record','final','dose', ...
    'delivery','solver','supportReduced','supportFull','-v7.3');
end
function [obj,grad]=estimate_pg_work(state,bdtHistory)
bt=sum(state.history.backtrackingTrials);n=state.iterations;segments=numel(bdtHistory.iteration);
obj=segments+2*n+bt;grad=2*n;
end
function row=empty_method()
row=struct('method','','selectorStatus','','evaluationCountType','', ...
    'selectorWallSeconds',NaN, ...
    'polishWallSeconds',NaN,'objectiveEvaluations',NaN,'gradientEvaluations',NaN, ...
    'selectedSpots',NaN,'selectedLayers',NaN,'finalActiveSpots',NaN, ...
    'finalActiveLayers',NaN,'clinicalObjective',NaN,'maxConstraintViolation',NaN, ...
    'targetD95',NaN,'targetD98',NaN,'oarDmax',NaN,'discreteBDT',NaN, ...
    'supportPreserved',false,'finalIpoptStatus','');
end
function tbl=dose_table(a,b)
n=numel(a.structures);r=repmat(struct('structure','','type','','isTarget',false, ...
    'stage1D98',NaN,'finalD98',NaN,'stage1D95',NaN,'finalD95',NaN, ...
    'stage1Dmean',NaN,'finalDmean',NaN,'stage1Dmax',NaN,'finalDmax',NaN),n,1);
for i=1:n,r(i).structure=a.structures(i).name;r(i).type=a.structures(i).type;r(i).isTarget=a.structures(i).isTarget;for c={'D98','D95','Dmean','Dmax'},q=c{1};r(i).(['stage1' q])=b.structures(i).(q);r(i).(['final' q])=a.structures(i).(q);end,end
tbl=struct2table(r);
end
function n=count_layers(mask,groups),n=0;for g=1:numel(groups),n=n+any(mask(groups(g).ix));end,end
function bdt=reduced_bdt(w,ix,nFull,map,groups,cfg,machine),x=zeros(nFull,1);x(ix)=w;d=compute_delivery_metrics(x,map,groups,cfg,machine);bdt=d.T_total;end
function D=nominal_matrix(dij),if iscell(dij.physicalDose),D=dij.physicalDose{1};else,D=dij.physicalDose;end,end
function verify_sha256(f,e),[s,o]=system(sprintf('sha256sum "%s"',f));assert(s==0&&strcmpi(strtok(strtrim(o)),e),'Frozen checksum mismatch.');end
function create_empty_output_directory(d),if exist(d,'dir')==7,e=dir(d);e=e(~ismember({e.name},{'.','..'}));assert(isempty(e),'Output directory not empty.');else,mkdir(d);end,end
function add_paths(root,here),f={'config','planning','mapping','objective','delivery','reduction','evaluation','integrated','robustness'};for i=1:numel(f),addpath(fullfile(root,f{i}));end,addpath(here);end
function m=merge_structs_local(a,b),m=a;n=fieldnames(b);for i=1:numel(n),k=n{i};if isfield(m,k)&&isstruct(m.(k))&&isstruct(b.(k)),m.(k)=merge_structs_local(m.(k),b.(k));else,m.(k)=b.(k);end,end,end
