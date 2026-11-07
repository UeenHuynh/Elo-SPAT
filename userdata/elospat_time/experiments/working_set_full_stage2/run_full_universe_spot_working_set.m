function output=run_full_universe_spot_working_set(outputDir)
%RUN_FULL_UNIVERSE_SPOT_WORKING_SET Stage-2 exact-clinical spot pricing.
narginchk(1,1);here=fileparts(mfilename('fullpath'));root=fileparts(fileparts(here));
add_paths(root,here);create_empty(outputDir);
artifactFile=fullfile(root,'results','task25_27_pipeline_20deg.mat');
stage1File=fullfile(root,'results','full_integrated_layer_20deg_job22180', ...
    'alpha_0p3','alpha_0p3_final.mat');
verify_sha256(stage1File,'9c77e38f26f43515fc1515a29523b94ec2c52265df74e9686f0196efc7097f18');
source=load(artifactFile,'cfg','machine','denseResult','dij','pln','map');
stage1=load(stage1File,'record','final','finalDose','finalDelivery');
patient=load(source.cfg.patientFile,'ct','cst');cfg=merge_structs(default_config(),source.cfg);
cfg.lambdaB=0;cfg.lambdaES=0;cfg.lambdaT=0;
wStage1=reshape(double(stage1.final.w),[],1);allowedFull=wStage1>0;
assert(nnz(allowedFull)==136409,'Frozen Stage-1 allowed universe changed.');
[dijAllowed,allowedToFull,fullToAllowed,mapAllowed]= ...
    reduce_dij_to_support(source.dij,allowedFull,source.map);
for i=1:numel(mapAllowed),mapAllowed(i).originalGlobalIx=mapAllowed(i).globalIx;mapAllowed(i).globalIx=i;end
[groupsAllowed,cpGroupsAllowed]=build_group_index(mapAllowed);
assert(numel(groupsAllowed)==748&&numel(cpGroupsAllowed)==18);
[fullGroups,~]=build_group_index(source.map);
layerAllowed=false(numel(fullGroups),1);
for q=1:numel(fullGroups),layerAllowed(q)=any(allowedFull(fullGroups(q).ix));end
assert(nnz(~layerAllowed)==20&&nnz(layerAllowed)==748);

fullModel=build_matrad_dose_model(dijAllowed,patient.cst,source.pln,0);
wAllowed=wStage1(allowedToFull);
[stage1Objective,~,stage1Diagnostics]=evaluate_matrad_dose_model(fullModel,wAllowed,false);
g0=reshape(double(fullModel.optiProb.matRad_objectiveGradient( ...
    zeros(size(wAllowed)),fullModel.dij,fullModel.cst)),[],1);
[initialMask,seedDiagnostics]=build_cp_gradient_seed(g0,wAllowed,mapAllowed,500);
assert(nnz(initialMask)==9000&&numel(unique([mapAllowed(initialMask).cpIndex]))==18);
initialW=zeros(size(wAllowed));initialW(initialMask)=wAllowed(initialMask);
D=nominal_matrix(source.dij);
problem=struct('dij',dijAllowed,'cst',{patient.cst},'pln',source.pln, ...
    'fullModel',fullModel,'groups',groupsAllowed,'initialW',initialW);
problem.doseFunction=@dose_metrics_allowed;
problem.bdtFunction=@delivery_allowed;
problem.meaningfulFunction=@meaningful_gate;
inner=struct('desiredDualTolerance',1e-3,'acceptableDualTolerance',1e-2, ...
    'acceptableObjectiveWindowTolerance',1e-6,'stabilityWindow',8, ...
    'ipopt',struct('max_iter',500,'max_cpu_time',14400,'print_level',3, ...
    'hessian_approximation','limited-memory','limited_memory_max_history',50));
settings=struct('spotBatchSize',1000,'maxOuterIterations',20, ...
    'absoluteTolerance',1e-7,'relativeTolerance',1e-6, ...
    'zeroTolerance',1e-10,'removeZeroNonbasic',true, ...
    'initialGradientEvaluations',1,'innerSettings',inner,'outputDir',outputDir);
setup=struct('classification','full Stage-1-universe exact-clinical spot working set', ...
    'allowedSpots',nnz(allowedFull),'allowedLayers',nnz(layerAllowed), ...
    'permanentlyExcludedLayers',nnz(~layerAllowed),'initialSpots',nnz(initialMask), ...
    'initialPerControlPoint',500,'spotBatchSize',settings.spotBatchSize, ...
    'stage1Objective',stage1Objective,'stage1ConstraintViolation', ...
    stage1Diagnostics.maxConstraintViolation,'stage1Checksum', ...
    '9c77e38f26f43515fc1515a29523b94ec2c52265df74e9686f0196efc7097f18', ...
    'otherOptimizersEnabled',false,'deliveryPenaltyEnabled',false);
save(fullfile(outputDir,'full_ws_setup.mat'),'setup','settings','allowedFull', ...
    'allowedToFull','fullToAllowed','layerAllowed','initialMask','seedDiagnostics','-v7.3');
writetable(struct2table(seedDiagnostics),fullfile(outputDir,'initial_seed_by_cp.csv'));
result=run_full_spot_working_set_solver(problem,settings);

last=result.history(end);final=[];finalDose=[];finalDelivery=[];finalSolver=[];
planAccepted=false;supportPreserved=false;doseQualityPass=false;
if last.clinicallyMeaningful
    supportAllowed=result.w>settings.zeroTolerance;
    supportFull=false(size(wStage1));supportFull(allowedToFull(supportAllowed))=true;
    assert(~any(supportFull(~allowedFull)),'Excluded Stage-1 spot became available.');
    wInit=zeros(size(wStage1));wInit(allowedToFull(supportAllowed))= ...
        max(result.w(supportAllowed),cfg.mmuThreshold);
    final=run_final_reoptimization(source.dij,patient.cst,source.pln,wInit,cfg,supportFull,[]);
    tol=1e-8*max(1,max(abs(final.w)));supportPreserved=~any(abs(final.w(~supportFull))>tol);
    finalDose=compute_dvh_metrics(full(D*final.w),patient.ct,patient.cst,cfg);
    finalDelivery=compute_delivery_metrics(final.w,source.map,fullGroups,cfg,source.machine);
    finalSolver=compute_solver_metrics(final.resultGUI,struct(),final.optimizer);
    [doseQualityPass,doseTable,worstAdverse]=dose_quality(finalDose,stage1.finalDose, ...
        cfg.evaluation.ctvD95TolerancePercent); %#ok<ASGLU>
    ipoptAccepted=finalSolver.converged||finalSolver.acceptableButNotConverged|| ...
        contains(lower(finalSolver.statusMessage),'acceptable');
    masterConverged=strcmp(result.stopReason,'restricted_and_pricing_converged');
    planAccepted=masterConverged&&doseQualityPass&&supportPreserved&&ipoptAccepted;
    writetable(doseTable,fullfile(outputDir,'full_ws_final_dose.csv'));
    save(fullfile(outputDir,'full_ws_final_polish.mat'),'final','finalDose', ...
        'finalDelivery','finalSolver','supportAllowed','supportFull', ...
        'supportPreserved','doseQualityPass','planAccepted','worstAdverse','-v7.3');
end
summary=struct('stopReason',result.stopReason,'outerIterations',result.outerIterations, ...
    'activeSpots',last.positiveSpots,'activeLayers',last.activeLayers, ...
    'restrictedSolveAccepted',last.restrictedSolveAccepted, ...
    'pricingViolation',last.maximumInactiveSpotViolation, ...
    'pricingTolerance',last.pricingTolerance,'objective',last.objective, ...
    'targetD95',last.targetD95,'targetD98',last.targetD98,'oarDmax',last.oarDmax, ...
    'provisionalBDT',last.provisionalBDT,'clinicallyMeaningful',last.clinicallyMeaningful, ...
    'masterConverged',strcmp(result.stopReason,'restricted_and_pricing_converged'), ...
    'finalPolishRun',~isempty(final),'supportPreserved',supportPreserved, ...
    'doseQualityPass',doseQualityPass,'planAccepted',planAccepted);
writetable(struct2table(summary,'AsArray',true),fullfile(outputDir,'full_ws_summary.csv'));
save(fullfile(outputDir,'full_ws_complete.mat'),'setup','result','summary','-v7.3');
output=struct('setup',setup,'result',result,'summary',summary,'outputDir',outputDir);

    function dose=dose_metrics_allowed(x)
        z=zeros(size(wStage1));z(allowedToFull)=x;
        dose=compute_dvh_metrics(full(D*z),patient.ct,patient.cst,cfg);
    end
    function delivery=delivery_allowed(x)
        z=zeros(size(wStage1));z(allowedToFull)=x;
        delivery=compute_delivery_metrics(z,source.map,fullGroups,cfg,source.machine);
    end
    function pass=meaningful_gate(dose,diagnostics,objective)
        pass=dose.summary.targetD95Min>=0.95*stage1.finalDose.summary.targetD95Min&& ...
            dose.summary.targetD98Min>=0.95*stage1.finalDose.summary.targetD98Min&& ...
            dose.summary.oarDmaxMax<=1.10*stage1.finalDose.summary.oarDmaxMax&& ...
            diagnostics.maxConstraintViolation<=1e-6&&objective<=1.10*stage1Objective;
    end
end

function [pass,tbl,worst]=dose_quality(candidate,reference,tolerance)
n=numel(candidate.structures);r=repmat(struct('structure','','type','', ...
    'referenceD98',NaN,'finalD98',NaN,'referenceD95',NaN,'finalD95',NaN, ...
    'referenceDmean',NaN,'finalDmean',NaN,'referenceDmax',NaN,'finalDmax',NaN, ...
    'adverseChangePercent',NaN,'structurePass',false),n,1);worst=0;
for i=1:n
    a=reference.structures(i);b=candidate.structures(i);r(i).structure=a.name;r(i).type=a.type;
    for q={'D98','D95','Dmean','Dmax'},x=q{1};r(i).(['reference' x])=a.(x);r(i).(['final' x])=b.(x);end
    if a.isTarget,adverse=max([0,-pct(b.D95,a.D95),-pct(b.D98,a.D98)]);
    elseif strcmpi(a.type,'OAR'),adverse=max([0,pct(b.Dmean,a.Dmean),pct(b.Dmax,a.Dmax)]);
    else,adverse=0;end
    r(i).adverseChangePercent=adverse;r(i).structurePass=adverse<=tolerance;worst=max(worst,adverse);
end
tbl=struct2table(r);pass=all([r.structurePass]);
end
function value=pct(x,ref),value=100*(double(x)-double(ref))/max(abs(double(ref)),eps);end
function D=nominal_matrix(dij),if iscell(dij.physicalDose),D=dij.physicalDose{1};else,D=dij.physicalDose;end,end
function verify_sha256(f,e),[s,o]=system(sprintf('sha256sum "%s"',f));assert(s==0&&strcmpi(strtok(strtrim(o)),e));end
function create_empty(d),if exist(d,'dir')==7,e=dir(d);e=e(~ismember({e.name},{'.','..'}));assert(isempty(e));else,mkdir(d);end,end
function add_paths(root,here),f={'config','planning','mapping','objective','delivery','reduction','evaluation','integrated','robustness'};for i=1:numel(f),addpath(fullfile(root,f{i}));end,addpath(fullfile(root,'experiments','working_set_stage2'));addpath(here);end
function m=merge_structs(a,b),m=a;f=fieldnames(b);for i=1:numel(f),k=f{i};if isfield(m,k)&&isstruct(m.(k))&&isstruct(b.(k)),m.(k)=merge_structs(m.(k),b.(k));else,m.(k)=b.(k);end,end,end
