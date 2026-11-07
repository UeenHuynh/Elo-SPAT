function output=run_coreset_initialized_full_working_set(outputDir)
%RUN_CORESET_INITIALIZED_FULL_WORKING_SET Coreset gate, then clinical pricing.
narginchk(1,1);here=fileparts(mfilename('fullpath'));root=fileparts(fileparts(here));
add_paths(root,here);create_empty(outputDir);
artifact=load(fullfile(root,'results','task25_27_pipeline_20deg.mat'), ...
    'cfg','machine','dij','pln','map');
stage1File=fullfile(root,'results','full_integrated_layer_20deg_job22180', ...
    'alpha_0p3','alpha_0p3_final.mat');
verify_sha256(stage1File,'9c77e38f26f43515fc1515a29523b94ec2c52265df74e9686f0196efc7097f18');
stage1=load(stage1File,'record','final','finalDose','finalDelivery');
patient=load(artifact.cfg.patientFile,'ct','cst');cfg=merge_structs(default_config(),artifact.cfg);
cfg.lambdaB=0;cfg.lambdaES=0;cfg.lambdaT=0;
wFullRef=reshape(double(stage1.final.w),[],1);allowedFull=wFullRef>0;
assert(nnz(allowedFull)==136409);
[dijAllowed,allowedToFull,fullToAllowed,mapAllowed]= ...
    reduce_dij_to_support(artifact.dij,allowedFull,artifact.map);
for i=1:numel(mapAllowed),mapAllowed(i).originalGlobalIx=mapAllowed(i).globalIx;mapAllowed(i).globalIx=i;end
[groupsAllowed,cpGroups]=build_group_index(mapAllowed);assert(numel(groupsAllowed)==748&&numel(cpGroups)==18);
[groupsFull,~]=build_group_index(artifact.map);
Dfull=nominal_matrix(artifact.dij);D=nominal_matrix(dijAllowed);wRef=wFullRef(allowedToFull);
dRef=D*wRef;nVox=size(D,1);
[voxelWeights,weightTable]=build_clinical_reconstruction_weights(patient.cst,nVox);
clinicalVoxels=find(voxelWeights>0);sqrtW=sqrt(voxelWeights(clinicalVoxels));
weightedD=spdiags(sqrtW,0,numel(sqrtW),numel(sqrtW))*D(clinicalVoxels,:);
weightedReference=sqrtW.*dRef(clinicalVoxels);
columnNormSquared=full(sum(weightedD.^2,1))';
[initialMask,initialTargetTable]=build_target_cp_seed(D,dRef,patient.cst,mapAllowed,50);
fullModel=build_matrad_dose_model(dijAllowed,patient.cst,artifact.pln,0);
[referenceObjective,~,referenceDiagnostics]=evaluate_matrad_dose_model(fullModel,wRef,false);

problem=struct('D',D,'weightedD',weightedD,'weightedReference',weightedReference, ...
    'columnNormSquared',columnNormSquared,'referenceW',wRef,'initialMask',initialMask, ...
    'groups',groupsAllowed,'doseFunction',@dose_allowed, ...
    'objectiveFunction',@objective_allowed,'gateFunction',@gate);
coresetSettings=struct('maxOuterIterations',30,'batchSize',2000, ...
    'maximumSupportSize',60000,'minimumAdmissionScore',1e-14,'zeroTolerance',1e-10, ...
    'outputDir',outputDir,'nnls',struct('maxIterations',150, ...
    'stationarityTolerance',1e-5,'objectiveTolerance',1e-7, ...
    'stabilityWindow',6,'powerIterations',10));
setup=struct('classification','Stage-1 dose-preserving coreset initializer', ...
    'allowedSpots',numel(wRef),'allowedLayers',numel(groupsAllowed), ...
    'excludedLayers',20,'clinicalVoxels',numel(clinicalVoxels), ...
    'initialSpots',nnz(initialMask),'targetAggregateWeight',10, ...
    'oarAggregateWeight',1,'bodyAggregateWeight',0, ...
    'weightedErrorTolerance',0.02,'targetMetricTolerancePercent',2, ...
    'oarMetricTolerancePercent',2,'oarAbsoluteToleranceGy',0.2, ...
    'objectiveTolerancePercent',2,'referenceObjective',referenceObjective, ...
    'referenceConstraintViolation',referenceDiagnostics.maxConstraintViolation, ...
    'clinicalPricingConditionalOnCoresetGate',true,'otherOptimizersEnabled',false);
save(fullfile(outputDir,'coreset_setup.mat'),'setup','coresetSettings','allowedFull', ...
    'allowedToFull','fullToAllowed','clinicalVoxels','voxelWeights','initialMask','-v7.3');
writetable(weightTable,fullfile(outputDir,'clinical_reconstruction_weights.csv'));
writetable(initialTargetTable,fullfile(outputDir,'initial_target_cp_seed.csv'));
coreset=run_stage1_dose_coreset(problem,coresetSettings);
coresetPassed=strcmp(coreset.stopReason,'representativeness_gate_passed');
clinical=[];polish=[];polishDose=[];polishDelivery=[];polishSolver=[];
supportPreserved=false;dosePass=false;
if coresetPassed
    clinicalProblem=struct('dij',dijAllowed,'cst',{patient.cst},'pln',artifact.pln, ...
        'fullModel',fullModel,'groups',groupsAllowed,'initialW',coreset.w, ...
        'doseFunction',@dose_allowed,'bdtFunction',@delivery_allowed, ...
        'meaningfulFunction',@clinical_meaningful);
    inner=struct('ipopt',struct('max_iter',500,'max_cpu_time',14400, ...
        'print_level',3,'hessian_approximation','limited-memory', ...
        'limited_memory_max_history',50));
    clinicalSettings=struct('spotBatchSize',1000,'maxOuterIterations',15, ...
        'absoluteTolerance',1e-7,'relativeTolerance',1e-6,'zeroTolerance',1e-10, ...
        'removeZeroNonbasic',true,'initialGradientEvaluations',0, ...
        'innerSettings',inner,'outputDir',outputDir);
    clinical=run_full_spot_working_set_solver(clinicalProblem,clinicalSettings);
    if clinical.history(end).clinicallyMeaningful,polishW=clinical.w;polishSource='clinical_working_set';
    else,polishW=coreset.w;polishSource='representative_coreset_fallback';end
    supportAllowed=polishW>1e-10;supportFull=false(size(wFullRef));
    supportFull(allowedToFull(supportAllowed))=true;wInit=zeros(size(wFullRef));
    wInit(allowedToFull(supportAllowed))=max(polishW(supportAllowed),cfg.mmuThreshold);
    polish=run_final_reoptimization(artifact.dij,patient.cst,artifact.pln,wInit,cfg,supportFull,[]);
    tol=1e-8*max(1,max(abs(polish.w)));supportPreserved=~any(abs(polish.w(~supportFull))>tol);
    polishDose=compute_dvh_metrics(full(Dfull*polish.w),patient.ct,patient.cst,cfg);
    polishDelivery=compute_delivery_metrics(polish.w,artifact.map,groupsFull,cfg,artifact.machine);
    polishSolver=compute_solver_metrics(polish.resultGUI,struct(),polish.optimizer);
    [dosePass,doseTable,worstAdverse]=final_dose_quality(polishDose,stage1.finalDose,2); %#ok<ASGLU>
    writetable(doseTable,fullfile(outputDir,'final_polish_dose.csv'));
    save(fullfile(outputDir,'final_polish.mat'),'polish','polishDose','polishDelivery', ...
        'polishSolver','polishSource','supportAllowed','supportFull', ...
        'supportPreserved','dosePass','worstAdverse','-v7.3');
end
last=coreset.history(end);summary=struct('coresetPassed',coresetPassed, ...
    'coresetStopReason',coreset.stopReason,'coresetIterations',coreset.iterations, ...
    'coresetSpots',last.activeSpots,'coresetLayers',last.activeLayers, ...
    'weightedRelativeError',last.weightedRelativeError,'targetD95',last.targetD95, ...
    'targetD98',last.targetD98,'oarDmax',last.oarDmax, ...
    'clinicalPricingRun',~isempty(clinical),'clinicalMasterConverged', ...
    ~isempty(clinical)&&strcmp(clinical.stopReason,'restricted_and_pricing_converged'), ...
    'finalPolishRun',~isempty(polish),'supportPreserved',supportPreserved, ...
    'finalDosePass',dosePass);
writetable(struct2table(summary,'AsArray',true),fullfile(outputDir,'coreset_working_set_summary.csv'));
save(fullfile(outputDir,'coreset_working_set_complete.mat'),'setup','coreset','clinical','summary','-v7.3');
output=struct('setup',setup,'coreset',coreset,'clinical',clinical,'summary',summary);

    function dose=dose_allowed(x),dose=compute_dvh_metrics(full(D*x),patient.ct,patient.cst,cfg);end
    function objective=objective_allowed(x),objective=evaluate_matrad_dose_model(fullModel,x,false);end
    function delivery=delivery_allowed(x)
        z=zeros(size(wFullRef));z(allowedToFull)=x;
        delivery=compute_delivery_metrics(z,artifact.map,groupsFull,cfg,artifact.machine);
    end
    function result=gate(relError,dose,objective)
        [targetOK,oarOK,reason]=representativeness_metrics(dose,stage1.finalDose,2,0.2);
        errorOK=relError<=0.02;objectiveOK=objective<=1.02*referenceObjective;
        pass=errorOK&&targetOK&&oarOK&&objectiveOK;
        result=struct('pass',pass,'weightedErrorPass',errorOK, ...
            'targetPass',targetOK,'oarPass',oarOK,'objectivePass',objectiveOK, ...
            'reason',sprintf('error=%d,target=%d,oar=%d,objective=%d;%s', ...
            errorOK,targetOK,oarOK,objectiveOK,reason));
    end
    function pass=clinical_meaningful(dose,diagnostics,objective)
        [targetOK,oarOK]=representativeness_metrics(dose,stage1.finalDose,2,0.2);
        pass=targetOK&&oarOK&&objective<=1.02*referenceObjective&& ...
            diagnostics.maxConstraintViolation<=1e-6;
    end
end

function [targetOK,oarOK,reason]=representativeness_metrics(c,r,pctTol,absTol)
targetOK=true;oarOK=true;failed={};
for i=1:numel(c.structures)
    a=r.structures(i);b=c.structures(i);
    if a.isTarget
        ok=b.D95>=a.D95*(1-pctTol/100)&&b.D98>=a.D98*(1-pctTol/100);targetOK=targetOK&&ok;
    elseif strcmpi(a.type,'OAR')&&~strcmpi(a.name,'BODY')&&~strcmpi(a.name,'EXTERNAL')
        meanTol=max(absTol,pctTol/100*max(abs(a.Dmean),1));
        maxTol=max(absTol,pctTol/100*max(abs(a.Dmax),1));
        ok=b.Dmean<=a.Dmean+meanTol&&b.Dmax<=a.Dmax+maxTol;oarOK=oarOK&&ok;
    else,ok=true;end
    if ~ok,failed{end+1}=a.name;end %#ok<AGROW>
end
reason=strjoin(failed,'|');
end
function [pass,tbl,worst]=final_dose_quality(c,r,tol)
n=numel(c.structures);rows=repmat(struct('structure','','type','','referenceD95',NaN, ...
    'finalD95',NaN,'referenceD98',NaN,'finalD98',NaN,'referenceDmean',NaN, ...
    'finalDmean',NaN,'referenceDmax',NaN,'finalDmax',NaN,'adversePercent',NaN, ...
    'pass',false),n,1);worst=0;
for i=1:n,a=r.structures(i);b=c.structures(i);rows(i).structure=a.name;rows(i).type=a.type;
    for q={'D95','D98','Dmean','Dmax'},x=q{1};rows(i).(['reference' x])=a.(x);rows(i).(['final' x])=b.(x);end
    if a.isTarget,ad=max([0,-pct(b.D95,a.D95),-pct(b.D98,a.D98)]);
    elseif strcmpi(a.type,'OAR'),ad=max([0,pct(b.Dmean,a.Dmean),pct(b.Dmax,a.Dmax)]);else,ad=0;end
    rows(i).adversePercent=ad;rows(i).pass=ad<=tol;worst=max(worst,ad);
end
tbl=struct2table(rows);pass=all([rows.pass]);
end
function x=pct(a,b),x=100*(a-b)/max(abs(b),eps);end
function D=nominal_matrix(dij),if iscell(dij.physicalDose),D=dij.physicalDose{1};else,D=dij.physicalDose;end,end
function verify_sha256(f,e),[s,o]=system(sprintf('sha256sum "%s"',f));assert(s==0&&strcmpi(strtok(strtrim(o)),e));end
function create_empty(d),if exist(d,'dir')==7,e=dir(d);e=e(~ismember({e.name},{'.','..'}));assert(isempty(e));else,mkdir(d);end,end
function add_paths(root,here),f={'config','planning','mapping','objective','delivery','reduction','evaluation','integrated','robustness'};for i=1:numel(f),addpath(fullfile(root,f{i}));end,addpath(fullfile(root,'experiments','working_set_stage2'));addpath(fullfile(root,'experiments','working_set_full_stage2'));addpath(here);end
function m=merge_structs(a,b),m=a;f=fieldnames(b);for i=1:numel(f),k=f{i};if isfield(m,k)&&isstruct(m.(k))&&isstruct(b.(k)),m.(k)=merge_structs(m.(k),b.(k));else,m.(k)=b.(k);end,end,end
