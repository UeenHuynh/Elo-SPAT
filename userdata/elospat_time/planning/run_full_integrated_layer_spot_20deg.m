function outputs = run_full_integrated_layer_spot_20deg(outputDir, alphaSGrid)
% Full 20-degree Stage-2 layer + smooth-L0/MM spot optimization.

narginchk(0,2);
root=fileparts(fileparts(mfilename('fullpath'))); add_module_paths(root);
if nargin<1 || isempty(outputDir), outputDir=fullfile(root,'results', ...
        ['full_integrated_layer_spot_20deg_' datestr(now,'yyyymmdd_HHMMSS')]); end
if nargin<2 || isempty(alphaSGrid), alphaSGrid=[0.01;0.03;0.1]; end
alphaSGrid=reshape(double(alphaSGrid),[],1); create_or_reuse_output_directory(outputDir);

artifactFile=fullfile(root,'results','task25_27_pipeline_20deg.mat');
frozenFile=fullfile(root,'results','full_integrated_layer_20deg_job22180', ...
    'alpha_0p3','alpha_0p3_final.mat');
scalingFile=fullfile(root,'results','full_integrated_spot_setup_20deg_job22187', ...
    'spot_setup_complete.mat');
expectedHash='9c77e38f26f43515fc1515a29523b94ec2c52265df74e9686f0196efc7097f18';
verify_sha256(frozenFile,expectedHash);
assert(exist(scalingFile,'file')==2,'Missing Stage-2 calibration job 22187.');
scaleSource=load(scalingFile,'setup');
source=load(artifactFile,'cfg','machine','denseResult','dij','pln','map','groups');
layer=load(frozenFile,'record','optimizer','final','finalDose', ...
    'finalDelivery','cfgTrial');
patient=load(source.cfg.patientFile,'ct','cst');
assert(abs(layer.record.alpha-0.3)<=1e-12 && layer.record.planAccepted);
assert(strcmpi(scaleSource.setup.frozenSha256,expectedHash));

cfg=merge_structs_local(default_config(),source.cfg); cfg.mode='integrated';
cfg.integrated.componentStage='spot'; cfg.lambdaEL=layer.record.lambdaEL;
cfg.lambdaB=0; cfg.lambdaES=0; cfg.lambdaT=0;
cfg.integrated.stationarityStoppingEnabled=false;
cfg.integrated.maxIterations=50; cfg.integrated.minimumIterations=10;
cfg.integrated.verbose=true; cfg.captureConsoleOutput=false;
cfg.enforceSelectedMinimumWeight=true; cfg=configure_integrated_stage(cfg);
assert_stage2(cfg,layer.record.lambdaEL);

nFull=double(source.dij.totalNumOfBixels);
w0=reshape(double(layer.final.w),[],1);
denseW=reshape(double(source.denseResult.w),[],1);
assert(numel(w0)==nFull && numel(denseW)==nFull && numel(source.map)==nFull);
[groups,cpGroups]=build_group_index(source.map);
doseModel=build_matrad_dose_model(source.dij,patient.cst,source.pln, ...
    cfg.integrated.constraintPenalty);
context=build_integrated_context(groups,cpGroups,source.map,source.machine,cfg);
cfgScale=cfg; cfgScale.lambdaS=0;
[~,parts0]=evaluate_integrated_objective(w0,doseModel,context,cfgScale);
spotScale=parts0.dose/max(parts0.spotSmoothL0,eps);
scaleRelativeError=abs(spotScale-scaleSource.setup.spotScale)/ ...
    max(1,abs(scaleSource.setup.spotScale));
assert(scaleRelativeError<=1e-11,'Spot scale disagrees with calibration by %.3e.', ...
    scaleRelativeError);
lambdaSGrid=zeros(size(alphaSGrid));
for i=1:numel(alphaSGrid)
    q=find(abs(scaleSource.setup.alphaSGrid-alphaSGrid(i))<=1e-12,1);
    assert(~isempty(q),'alphaS %.12g absent from calibration.',alphaSGrid(i));
    lambdaSGrid(i)=scaleSource.setup.lambdaS(q);
end

rule=default_stabilization_rule(); maxTotalIterations=600;
D=nominal_matrix(source.dij);
denseDose=dose_metrics(D,denseW,patient,cfg);
denseDelivery=compute_delivery_metrics(denseW,source.map,groups,cfg,source.machine);
layerDose=layer.finalDose; layerDelivery=layer.finalDelivery;
setup=struct('artifactFile',artifactFile,'frozenFile',frozenFile, ...
    'frozenSha256',expectedHash,'scalingFile',scalingFile, ...
    'nWeights',nFull,'alphaEL',layer.record.alpha,'lambdaEL',cfg.lambdaEL, ...
    'spotScale',spotScale,'alphaSGrid',alphaSGrid,'lambdaSGrid',lambdaSGrid, ...
    'epsS',cfg.epsS,'mmuThreshold',cfg.mmuThreshold, ...
    'checkpointInterval',50,'maxTotalIterations',maxTotalIterations, ...
    'rule',rule,'penalties',struct('lambdaB',0,'lambdaES',0,'lambdaT',0));
save(fullfile(outputDir,'stage2_setup.mat'),'setup','-v7.3');
writetable(table(alphaSGrid,lambdaSGrid), ...
    fullfile(outputDir,'stage2_lambda_scaling.csv'));

records=repmat(empty_record(),numel(alphaSGrid),1);
for a=1:numel(alphaSGrid)
    alphaS=alphaSGrid(a); token=alpha_token(alphaS);
    cfgTrial=cfg; cfgTrial.lambdaS=lambdaSGrid(a);
    [state,bdtHistory,stabilization,resumeFile]= ...
        load_latest_checkpoint(outputDir,token,alphaS,cfgTrial.lambdaS,cfgTrial.lambdaEL);
    try
        while isempty(state) || (state.iterations<maxTotalIterations && ...
                ~stabilization.numerically_stabilized)
            if isempty(state)
                state=run_integrated_optimizer(doseModel,context,w0,cfgTrial);
            else
                state=run_integrated_optimizer(doseModel,context,state.w,cfgTrial,state);
            end
            delivery=compute_delivery_metrics(state.w,source.map,groups,cfgTrial,source.machine);
            bdtHistory.iteration(end+1,1)=state.iterations;
            bdtHistory.value(end+1,1)=delivery.T_total;
            stabilization=evaluate_numerical_stabilization(state.history,bdtHistory,rule);
            optimizer=state; %#ok<NASGU>
            save(fullfile(outputDir,sprintf('alphaS_%s_iter_%04d.mat', ...
                token,state.iterations)),'optimizer','bdtHistory','stabilization', ...
                'cfgTrial','alphaS','scalingFile','frozenFile','-v7.3');
            if strcmp(state.status,'line_search_failed') || ...
                    strcmp(state.status,'numerical_stagnation'), break; end
        end

        r=empty_record(); r.alphaS=alphaS; r.lambdaS=cfgTrial.lambdaS;
        r.lambdaEL=cfgTrial.lambdaEL; r.iterations=state.iterations;
        r.optimizerStatus=state.status; r.stoppingReason=state.stoppingReason;
        r.objective=state.objective;
        r.normalizedStationarity=state.finalNormalizedStationarity;
        r.numericallyStabilized=stabilization.numerically_stabilized;
        r.stabilizationReason=stabilization.reason;
        r.layerSupportStable=stabilization.supportStable;
        r.spotSupportStable=stabilization.spotSupportStable;
        r.exactZeroLayers=state.history.exactZeroLayers(end);
        r.activeLayers=state.history.exactActiveLayers(end);
        r.activeSpots=nnz(state.w>cfg.spotThreshold);
        r.exactZeroSpots=nnz(state.w==0);
        r.positiveSpots=nnz(state.w>0);
        r.subMmuSpots=nnz(state.w>0 & state.w<cfg.mmuThreshold);
        r.minimumPositiveWeight=min_positive(state.w);
        r.integratedBDT=bdtHistory.value(end);
        r.maxConstraintViolation=state.history.maxConstraintViolation(end);
        r.bdtChangeFromDensePercent=pct(r.integratedBDT,denseDelivery.T_total);
        r.bdtChangeFromLayerPercent=pct(r.integratedBDT,layerDelivery.T_total);

        if stabilization.numerically_stabilized
            extracted=extract_integrated_support(state.w,w0,context,cfgTrial);
            supportMask=reshape(logical(extracted.mask),[],1);
            final=run_final_reoptimization(source.dij,patient.cst,source.pln, ...
                extracted.w,cfgTrial,supportMask,extracted.layerMask);
            supportPreserved=~any(abs(final.w(~supportMask))>rule.supportLeakageTolerance);
            finalDose=dose_metrics(D,final.w,patient,cfgTrial);
            finalDelivery=compute_delivery_metrics(final.w,source.map,groups,cfgTrial,source.machine);
            finalSolver=compute_solver_metrics(final.resultGUI,struct(),final.optimizer);
            [dosePass,doseTable,worstDense,worstLayer]=dose_comparison( ...
                finalDose,denseDose,layerDose,cfg.evaluation.ctvD95TolerancePercent);
            ipoptAcceptable=finalSolver.converged || ...
                finalSolver.acceptableButNotConverged || ...
                contains(lower(finalSolver.statusMessage),'acceptable');
            acceptance=evaluate_plan_acceptance(dosePass,r.maxConstraintViolation, ...
                supportPreserved,ipoptAcceptable,rule);
            r.planAccepted=acceptance.plan_accepted; r.doseQualityPass=dosePass;
            r.supportPreserved=supportPreserved; r.finalIpoptStatus=finalSolver.statusMessage;
            r.finalActiveLayers=finalDelivery.nActiveLayers;
            r.finalActiveSpots=finalDelivery.nActiveSpots;
            r.finalExactZeroSpots=nnz(final.w==0);
            r.finalSubMmuSpots=nnz(final.w>0 & final.w<cfg.mmuThreshold-rule.supportLeakageTolerance);
            r.finalMinimumPositiveWeight=min_positive(final.w);
            r.minimumWeightPass=r.finalSubMmuSpots==0 && ...
                r.finalMinimumPositiveWeight>=cfg.mmuThreshold-rule.supportLeakageTolerance;
            r.finalBDT=finalDelivery.T_total;
            r.bdtChangeFromDensePercent=pct(r.finalBDT,denseDelivery.T_total);
            r.bdtChangeFromLayerPercent=pct(r.finalBDT,layerDelivery.T_total);
            r.worstAdverseVsDensePercent=worstDense;
            r.worstAdverseVsLayerPercent=worstLayer;
            writetable(doseTable,fullfile(outputDir,sprintf( ...
                'alphaS_%s_final_ipopt_dose.csv',token)));
            record=r; %#ok<NASGU>
            save(fullfile(outputDir,sprintf('alphaS_%s_final.mat',token)), ...
                'record','optimizer','bdtHistory','stabilization','extracted', ...
                'supportMask','final','finalDose','finalDelivery','finalSolver', ...
                'acceptance','doseTable','cfgTrial','scalingFile','frozenFile', ...
                'resumeFile','-v7.3');
        end
        records(a)=r;
    catch err
        records(a)=empty_record(); records(a).alphaS=alphaS;
        records(a).lambdaS=cfgTrial.lambdaS; records(a).lambdaEL=cfgTrial.lambdaEL;
        records(a).runStatus=['failed: ' err.identifier];
        errorInfo=struct('identifier',err.identifier,'message',err.message,'stack',err.stack); %#ok<NASGU>
        save(fullfile(outputDir,sprintf('alphaS_%s_failure.mat',token)), ...
            'errorInfo','cfgTrial','-v7.3');
        warning('Stage-2 alphaS=%g failed: %s',alphaS,err.message);
    end
    writetable(struct2table(records(1:a),'AsArray',true), ...
        fullfile(outputDir,'stage2_summary.csv'));
end
save(fullfile(outputDir,'stage2_complete.mat'),'records','setup', ...
    'denseDose','denseDelivery','layerDose','layerDelivery','-v7.3');
outputs=struct('records',records,'setup',setup,'outputDir',outputDir);
end

function assert_stage2(cfg,lambdaEL)
assert(strcmpi(cfg.integrated.componentStage,'spot'));
assert(abs(cfg.lambdaEL-lambdaEL)<=eps(max(1,abs(lambdaEL))));
assert(cfg.lambdaB==0 && cfg.lambdaES==0 && cfg.lambdaT==0);
assert(~cfg.integrated.stationarityStoppingEnabled);
end

function r=empty_record()
r=struct('alphaS',NaN,'lambdaS',NaN,'lambdaEL',NaN,'iterations',NaN, ...
    'optimizerStatus','','stoppingReason','','objective',NaN, ...
    'normalizedStationarity',NaN,'numericallyStabilized',false, ...
    'stabilizationReason','','layerSupportStable',false,'spotSupportStable',false, ...
    'exactZeroLayers',NaN,'activeLayers',NaN,'activeSpots',NaN, ...
    'exactZeroSpots',NaN,'positiveSpots',NaN,'subMmuSpots',NaN, ...
    'minimumPositiveWeight',NaN,'integratedBDT',NaN,'maxConstraintViolation',NaN, ...
    'finalActiveLayers',NaN,'finalActiveSpots',NaN,'finalExactZeroSpots',NaN, ...
    'finalSubMmuSpots',NaN,'finalMinimumPositiveWeight',NaN, ...
    'minimumWeightPass',false,'supportPreserved',false,'finalIpoptStatus','', ...
    'finalBDT',NaN,'bdtChangeFromDensePercent',NaN,'bdtChangeFromLayerPercent',NaN, ...
    'worstAdverseVsDensePercent',NaN,'worstAdverseVsLayerPercent',NaN, ...
    'doseQualityPass',false,'planAccepted',false,'runStatus','completed');
end

function [pass,tbl,worstDense,worstLayer]=dose_comparison(candidate,dense,layer,tolerance)
n=numel(dense.structures); rows=repmat(struct('structure','','type','', ...
    'isTarget',false,'denseD98',NaN,'layerD98',NaN,'finalD98',NaN, ...
    'denseD95',NaN,'layerD95',NaN,'finalD95',NaN, ...
    'denseDmean',NaN,'layerDmean',NaN,'finalDmean',NaN, ...
    'denseDmax',NaN,'layerDmax',NaN,'finalDmax',NaN, ...
    'adverseVsDensePercent',NaN,'adverseVsLayerPercent',NaN, ...
    'structurePass',false),n,1); worstDense=0; worstLayer=0;
for i=1:n
    d=dense.structures(i); l=layer.structures(i); f=candidate.structures(i);
    rows(i).structure=d.name; rows(i).type=d.type; rows(i).isTarget=d.isTarget;
    names={'D98','D95','Dmean','Dmax'};
    for j=1:numel(names), k=names{j}; rows(i).(['dense' k])=d.(k); ...
            rows(i).(['layer' k])=l.(k); rows(i).(['final' k])=f.(k); end
    rows(i).adverseVsDensePercent=adverse_change(f,d);
    rows(i).adverseVsLayerPercent=adverse_change(f,l);
    rows(i).structurePass=rows(i).adverseVsDensePercent<=tolerance;
    worstDense=max(worstDense,rows(i).adverseVsDensePercent);
    worstLayer=max(worstLayer,rows(i).adverseVsLayerPercent);
end
tbl=struct2table(rows); pass=all([rows.structurePass]);
end
function v=adverse_change(f,r), if r.isTarget, v=max([0,-pct(f.D95,r.D95),-pct(f.D98,r.D98)]); elseif strcmpi(r.type,'OAR'), v=max([0,pct(f.Dmean,r.Dmean),pct(f.Dmax,r.Dmax)]); else, v=0; end, end
function value=min_positive(w), q=w(w>0); if isempty(q), value=NaN; else, value=min(q); end, end
function verify_sha256(file,expected), [s,o]=system(sprintf('sha256sum "%s"',file)); assert(s==0); assert(strcmpi(strtok(strtrim(o)),expected),'Frozen Stage-1 checksum mismatch.'); end
function [state,bdtHistory,stabilization,resumeFile]=load_latest_checkpoint(d,token,alphaS,lambdaS,lambdaEL), state=[]; bdtHistory=struct('iteration',zeros(0,1),'value',zeros(0,1)); stabilization=struct('numerically_stabilized',false,'reason','not_started'); resumeFile=''; fs=dir(fullfile(d,sprintf('alphaS_%s_iter_*.mat',token))); if isempty(fs), return; end; it=zeros(numel(fs),1); for i=1:numel(fs), p=regexp(fs(i).name,'_iter_(\d+)\.mat$','tokens','once'); if ~isempty(p), it(i)=str2double(p{1}); end, end; [~,q]=max(it); resumeFile=fullfile(fs(q).folder,fs(q).name); x=load(resumeFile,'optimizer','bdtHistory','stabilization','cfgTrial','alphaS'); assert(abs(x.alphaS-alphaS)<=1e-12); assert(abs(x.cfgTrial.lambdaS-lambdaS)<=1e-12*max(1,abs(lambdaS))); assert(abs(x.cfgTrial.lambdaEL-lambdaEL)<=1e-12*max(1,abs(lambdaEL))); state=x.optimizer; bdtHistory=x.bdtHistory; stabilization=x.stabilization; fprintf('RESUME alphaS=%g from iteration %d: %s\n',alphaS,state.iterations,resumeFile); end
function value=pct(x,ref), value=100*(double(x)-double(ref))/max(abs(double(ref)),eps); end
function D=nominal_matrix(dij), if iscell(dij.physicalDose), D=dij.physicalDose{1}; else, D=dij.physicalDose; end, end
function m=dose_metrics(D,w,patient,cfg), m=compute_dvh_metrics(full(D*reshape(double(w),[],1)),patient.ct,patient.cst,cfg); end
function token=alpha_token(alpha), token=strrep(sprintf('%.4g',alpha),'.','p'); end
function add_module_paths(root), f={'config','planning','mapping','objective','delivery','reduction','robustness','evaluation','integrated'}; for i=1:numel(f), addpath(fullfile(root,f{i})); end, end
function create_or_reuse_output_directory(d), if exist(d,'dir')~=7, mkdir(d); end, end
function merged=merge_structs_local(base,override), merged=base; n=fieldnames(override); for i=1:numel(n), k=n{i}; if isfield(merged,k)&&isstruct(merged.(k))&&isstruct(override.(k)), merged.(k)=merge_structs_local(merged.(k),override.(k)); else, merged.(k)=override.(k); end, end, end
