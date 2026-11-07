function outputs = run_stage2_mm_schedule_benchmark_20deg(outputDir, mmPeriod)
% Compare smooth-L0 MM update schedules on the frozen full 20-degree case.
% This is an algorithmic screen only: no final IPOPT or operating-point
% selection is performed.

narginchk(2,2);
validPeriods=[1 5 20 50];
assert(isscalar(mmPeriod) && any(mmPeriod==validPeriods), ...
    'mmPeriod must be one of [1 5 20 50].');
root=fileparts(fileparts(mfilename('fullpath'))); add_module_paths(root);
create_empty_output_directory(outputDir);

artifactFile=fullfile(root,'results','task25_27_pipeline_20deg.mat');
frozenFile=fullfile(root,'results','full_integrated_layer_20deg_job22180', ...
    'alpha_0p3','alpha_0p3_final.mat');
scalingFile=fullfile(root,'results','full_integrated_spot_setup_20deg_job22187', ...
    'spot_setup_complete.mat');
expectedHash='9c77e38f26f43515fc1515a29523b94ec2c52265df74e9686f0196efc7097f18';
verify_sha256(frozenFile,expectedHash);
source=load(artifactFile,'cfg','machine','dij','pln','map');
layer=load(frozenFile,'record','final');
scaleSource=load(scalingFile,'setup');
patient=load(source.cfg.patientFile,'cst');
assert(layer.record.planAccepted && abs(layer.record.alpha-0.3)<=1e-12);

alphaS=0.03;
q=find(abs(scaleSource.setup.alphaSGrid-alphaS)<=1e-12,1);
assert(~isempty(q),'alphaS=0.03 is absent from saved calibration.');
cfg=merge_structs_local(default_config(),source.cfg); cfg.mode='integrated';
cfg.integrated.componentStage='spot';
cfg.lambdaEL=layer.record.lambdaEL;
cfg.lambdaS=scaleSource.setup.lambdaS(q);
cfg.lambdaB=0; cfg.lambdaES=0; cfg.lambdaT=0;
cfg.integrated.spotMMUpdatePeriod=mmPeriod;
cfg.integrated.maxIterations=25;
cfg.integrated.stationarityStoppingEnabled=false;
cfg.integrated.verbose=true; cfg.integrated.printEvery=25;
cfg.captureConsoleOutput=false;
cfg=configure_integrated_stage(cfg);
validate_config(cfg,source.machine);

w0=reshape(double(layer.final.w),[],1);
assert(numel(w0)==double(source.dij.totalNumOfBixels));
[groups,cpGroups]=build_group_index(source.map);
doseModel=build_matrad_dose_model(source.dij,patient.cst,source.pln, ...
    cfg.integrated.constraintPenalty);
context=build_integrated_context(groups,cpGroups,source.map,source.machine,cfg);
[initialObjective,initialParts]=evaluate_integrated_objective( ...
    w0,doseModel,context,cfg);

setup=struct('purpose','MM schedule diagnostic; no final IPOPT', ...
    'artifactFile',artifactFile,'frozenFile',frozenFile, ...
    'frozenSha256',expectedHash,'scalingFile',scalingFile, ...
    'alphaEL',layer.record.alpha,'lambdaEL',cfg.lambdaEL, ...
    'alphaS',alphaS,'lambdaS',cfg.lambdaS,'epsS',cfg.epsS, ...
    'mmPeriod',mmPeriod,'blockSize',25,'maximumUpdates',100, ...
    'nWeights',numel(w0),'initialObjective',initialObjective, ...
    'initialParts',initialParts,'penalties',struct( ...
    'lambdaB',0,'lambdaES',0,'lambdaT',0));
save(fullfile(outputDir,'benchmark_setup.mat'),'setup','cfg','-v7.3');

state=[]; previousW=w0; cumulativeSeconds=0;
rows=repmat(empty_row(),0,1);
for block=1:4
    timer=tic;
    if isempty(state)
        state=run_integrated_optimizer(doseModel,context,w0,cfg);
    else
        state=run_integrated_optimizer(doseModel,context,state.w,cfg,state);
    end
    segmentSeconds=toc(timer); cumulativeSeconds=cumulativeSeconds+segmentSeconds;
    delivery=compute_delivery_metrics( ...
        state.w,source.map,groups,cfg,source.machine);
    row=make_row(state,previousW,initialObjective,delivery,mmPeriod, ...
        segmentSeconds,cumulativeSeconds,cfg);
    rows(end+1,1)=row; %#ok<AGROW>
    optimizer=state; %#ok<NASGU>
    save(fullfile(outputDir,sprintf('mm_period_%02d_iter_%04d.mat', ...
        mmPeriod,state.iterations)),'optimizer','delivery','row','setup','cfg','-v7.3');
    writetable(struct2table(rows,'AsArray',true), ...
        fullfile(outputDir,'mm_schedule_blocks.csv'));
    fprintf(['MM BENCH period=%d iter=%d qUpdates=%d newZero=%d active=%d ' ...
        'J=%.10g dJ=%.6g stationarity=%.3e BDT=%.6f sec=%.1f\n'], ...
        mmPeriod,row.iteration,row.cumulativeMMUpdates,row.newExactZeroSpots, ...
        row.activeSpots,row.objective,row.objectiveDecreaseFromInitial, ...
        row.normalizedStationarity,row.provisionalBDT,row.cumulativeSeconds);
    previousW=state.w;
    if ~strcmp(state.status,'maximum_iterations')
        break;
    end
end

summary=rows(end); summaryStruct=summary; %#ok<NASGU>
save(fullfile(outputDir,'benchmark_complete.mat'),'summaryStruct','rows', ...
    'optimizer','setup','cfg','-v7.3');
writetable(struct2table(summary,'AsArray',true), ...
    fullfile(outputDir,'mm_schedule_summary.csv'));
outputs=struct('summary',summary,'rows',rows,'outputDir',outputDir);
end

function row=make_row(state,previousW,initialObjective,delivery,period,segSec,cumSec,cfg)
h=state.history; segmentStart=state.iterations-state.segmentIterations+1;
ix=h.iteration>=segmentStart;
oldZero=previousW==0; newZero=state.w==0;
row=empty_row(); row.mmPeriod=period; row.iteration=state.iterations;
row.segmentIterations=state.segmentIterations;
row.cumulativeMMUpdates=sum(h.spotMMUpdated~=0);
row.segmentMMUpdates=sum(h.spotMMUpdated(ix)~=0);
row.objective=state.objective;
row.objectiveDecreaseFromInitial=initialObjective-state.objective;
row.doseObjective=state.parts.dose;
row.layerPenalty=state.parts.layerL21half;
row.spotPenalty=state.parts.spotSmoothL0;
row.exactZeroLayers=h.exactZeroLayers(end);
row.activeLayers=h.exactActiveLayers(end);
row.exactZeroSpots=nnz(newZero); row.activeSpots=nnz(state.w>cfg.spotThreshold);
row.newExactZeroSpots=nnz(~oldZero & newZero);
row.reactivatedSpots=nnz(oldZero & ~newZero);
row.supportIdentity=support_hash(newZero);
row.normalizedStationarity=state.finalNormalizedStationarity;
row.minimumStep=min(h.stepSize(ix)); row.maximumStep=max(h.stepSize(ix));
row.backtrackingTrials=sum(h.backtrackingTrials(ix));
row.maxConstraintViolation=h.maxConstraintViolation(end);
row.provisionalBDT=delivery.T_total;
row.subMMUSpots=nnz(state.w>0 & state.w<cfg.mmuThreshold);
row.minimumPositiveWeight=min_positive(state.w);
row.segmentSeconds=segSec; row.cumulativeSeconds=cumSec;
row.optimizerStatus=state.status; row.stoppingReason=state.stoppingReason;
end

function r=empty_row()
r=struct('mmPeriod',NaN,'iteration',NaN,'segmentIterations',NaN, ...
    'cumulativeMMUpdates',NaN,'segmentMMUpdates',NaN,'objective',NaN, ...
    'objectiveDecreaseFromInitial',NaN,'doseObjective',NaN, ...
    'layerPenalty',NaN,'spotPenalty',NaN,'exactZeroLayers',NaN, ...
    'activeLayers',NaN,'exactZeroSpots',NaN,'activeSpots',NaN, ...
    'newExactZeroSpots',NaN,'reactivatedSpots',NaN,'supportIdentity','', ...
    'normalizedStationarity',NaN,'minimumStep',NaN,'maximumStep',NaN, ...
    'backtrackingTrials',NaN,'maxConstraintViolation',NaN, ...
    'provisionalBDT',NaN,'subMMUSpots',NaN,'minimumPositiveWeight',NaN, ...
    'segmentSeconds',NaN,'cumulativeSeconds',NaN, ...
    'optimizerStatus','','stoppingReason','');
end

function value=min_positive(w)
q=w(w>0); if isempty(q), value=NaN; else, value=min(q); end
end
function value=support_hash(mask)
md=java.security.MessageDigest.getInstance('SHA-256');
md.update(uint8(mask(:))); bytes=typecast(md.digest(),'uint8');
value=lower(reshape(dec2hex(bytes,2).',1,[]));
end
function verify_sha256(file,expected)
[s,o]=system(sprintf('sha256sum "%s"',file)); assert(s==0);
assert(strcmpi(strtok(strtrim(o)),expected),'Frozen Stage-1 checksum mismatch.');
end
function create_empty_output_directory(d)
if exist(d,'dir')==7
    e=dir(d); e=e(~ismember({e.name},{'.','..'}));
    assert(isempty(e),'Output directory must be empty: %s',d);
else
    mkdir(d);
end
end
function add_module_paths(root)
f={'config','planning','mapping','objective','delivery','reduction', ...
    'robustness','evaluation','integrated'};
for i=1:numel(f), addpath(fullfile(root,f{i})); end
end
function merged=merge_structs_local(base,override)
merged=base; n=fieldnames(override);
for i=1:numel(n)
    k=n{i};
    if isfield(merged,k)&&isstruct(merged.(k))&&isstruct(override.(k))
        merged.(k)=merge_structs_local(merged.(k),override.(k));
    else
        merged.(k)=override.(k);
    end
end
end
