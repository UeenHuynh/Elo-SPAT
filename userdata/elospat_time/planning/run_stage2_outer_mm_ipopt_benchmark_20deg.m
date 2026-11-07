function outputs = run_stage2_outer_mm_ipopt_benchmark_20deg(outputDir)
% First true sequential outer-MM/IPOPT benchmark for Stage 2.
% Stage-1 allowed support is fixed; R_EL is diagnostic only, not optimized.
% No final support polishing is performed.

narginchk(1,1);
root=fileparts(fileparts(mfilename('fullpath'))); add_module_paths(root);
create_empty_output_directory(outputDir);

artifactFile=fullfile(root,'results','task25_27_pipeline_20deg.mat');
frozenFile=fullfile(root,'results','full_integrated_layer_20deg_job22180', ...
    'alpha_0p3','alpha_0p3_final.mat');
scalingFile=fullfile(root,'results','full_integrated_spot_setup_20deg_job22187', ...
    'spot_setup_complete.mat');
pgRoot=fullfile(root,'results','stage2_mm_schedule_benchmark_20deg_job22249','period_1');
expectedHash='9c77e38f26f43515fc1515a29523b94ec2c52265df74e9686f0196efc7097f18';
verify_sha256(frozenFile,expectedHash);
source=load(artifactFile,'cfg','machine','dij','pln','map');
layer=load(frozenFile,'record','final');
scaleSource=load(scalingFile,'setup');
patient=load(source.cfg.patientFile,'cst');
assert(layer.record.planAccepted && abs(layer.record.alpha-0.3)<=1e-12);

cfg=merge_structs_local(default_config(),source.cfg); cfg.mode='integrated';
cfg.integrated.componentStage='spot'; cfg.lambdaEL=layer.record.lambdaEL;
cfg.lambdaB=0; cfg.lambdaES=0; cfg.lambdaT=0;
alphaS=0.03; ixAlpha=find(abs(scaleSource.setup.alphaSGrid-alphaS)<=1e-12,1);
assert(~isempty(ixAlpha)); cfg.lambdaS=scaleSource.setup.lambdaS(ixAlpha);
cfg=configure_integrated_stage(cfg);

wStage1=reshape(double(layer.final.w),[],1);
allowedMask=wStage1>0;
assert(nnz(allowedMask)==136409,'Unexpected frozen Stage-1 allowed support size.');
[dijReduced,reducedToFull,fullToReduced,mapReduced]= ...
    reduce_dij_to_support(source.dij,allowedMask,source.map);
% cut_dij_by_mask deliberately preserves original full-space globalIx for
% traceability. Group operators on reduced vectors require reduced-space
% indices, so retain the original index explicitly and remap globalIx.
for i=1:numel(mapReduced)
    mapReduced(i).originalGlobalIx=mapReduced(i).globalIx;
    mapReduced(i).globalIx=i;
end
wAnchor=wStage1(reducedToFull);
doseModel=build_matrad_dose_model( ...
    dijReduced,patient.cst,source.pln,cfg.integrated.constraintPenalty);
[groupsFull,~]=build_group_index(source.map);
[groupsReduced,~]=build_group_index(mapReduced);

innerSettings=struct();
innerSettings.ipopt=struct('max_iter',300,'max_cpu_time',21600, ...
    'tol',1e-6,'dual_inf_tol',1e-3,'constr_viol_tol',1e-6, ...
    'compl_inf_tol',1e-3,'acceptable_iter',8,'acceptable_tol',1e-4, ...
    'acceptable_constr_viol_tol',1e-5,'acceptable_dual_inf_tol',1e-2, ...
    'acceptable_compl_inf_tol',1e-2,'acceptable_obj_change_tol',1e-6, ...
    'print_level',5);
maxOuterCycles=3;
majorizerTolerance=1e-8;
projectionTolerance=1e-10;

pgReference=build_pg_reference(pgRoot,doseModel,reducedToFull,groupsReduced,cfg);
writetable(pgReference,fullfile(outputDir,'period1_pg_reference.csv'));
setup=struct('classification','sequential outer-MM/IPOPT', ...
    'innerObjective','F_matRad(w)+q''*w; Stage-1 R_EL omitted', ...
    'artifactFile',artifactFile,'frozenFile',frozenFile, ...
    'frozenSha256',expectedHash,'scalingFile',scalingFile, ...
    'period1ReferenceRoot',pgRoot,'alphaEL',layer.record.alpha, ...
    'lambdaELDiagnosticOnly',cfg.lambdaEL,'alphaS',alphaS, ...
    'lambdaS',cfg.lambdaS,'epsS',cfg.epsS,'allowedSpots',nnz(allowedMask), ...
    'fullSpots',numel(allowedMask),'nLayers',numel(groupsFull), ...
    'maxOuterCycles',maxOuterCycles,'innerSettings',innerSettings, ...
    'majorizerTolerance',majorizerTolerance, ...
    'projectionTolerance',projectionTolerance,'warmStart','primal_only', ...
    'finalIpoptPolishingRun',false,'penalties',struct( ...
    'lambdaELInInner',0,'lambdaB',0,'lambdaES',0,'lambdaT',0));
save(fullfile(outputDir,'outer_mm_setup.mat'),'setup','allowedMask', ...
    'reducedToFull','fullToReduced','pgReference','-v7.3');

rows=repmat(empty_cycle_row(),0,1); cumulativeWall=0;
cumulativeObjectiveCalls=0; cumulativeGradientCalls=0;
stopReason='maximum_outer_cycles_reached';
for outer=1:maxOuterCycles
    anchor=wAnchor;
    phiAnchor=anchor./(anchor+cfg.epsS);
    derivative=cfg.epsS./(anchor+cfg.epsS).^2;
    q=cfg.lambdaS*derivative;
    tangentConstant=cfg.lambdaS*sum(phiAnchor-derivative.*anchor);
    [fAnchor,~,anchorDiag]=evaluate_matrad_dose_model(doseModel,anchor,false);
    trueAnchor=fAnchor+cfg.lambdaS*sum(phiAnchor);
    majorizerAnchor=fAnchor+q'*anchor+tangentConstant;
    anchorEquality=abs(majorizerAnchor-trueAnchor);

    inner=solve_weighted_matrad_ipopt(doseModel,anchor,q,innerSettings);
    candidate=inner.w;
    [fCandidate,~,candidateDiag]=evaluate_matrad_dose_model( ...
        doseModel,candidate,false);
    phiCandidate=candidate./(candidate+cfg.epsS);
    trueCandidate=fCandidate+cfg.lambdaS*sum(phiCandidate);
    majorizerCandidate=fCandidate+q'*candidate+tangentConstant;
    dominanceGap=majorizerCandidate-trueCandidate;
    majorizerDecrease=majorizerAnchor-majorizerCandidate;
    trueDecrease=trueAnchor-trueCandidate;
    majorizerPass=anchorEquality<=majorizerTolerance*max(1,abs(trueAnchor)) && ...
        majorizerDecrease>=-majorizerTolerance*max(1,abs(majorizerAnchor)) && ...
        dominanceGap>=-majorizerTolerance*max(1,abs(trueCandidate));
    projectionPass=inner.projectionRelativeNorm<=projectionTolerance;
    cycleAccepted=inner.innerAccepted && majorizerPass && projectionPass;

    wFull=zeros(size(wStage1)); wFull(reducedToFull)=candidate;
    delivery=compute_delivery_metrics( ...
        wFull,source.map,groupsFull,cfg,source.machine);
    cumulativeWall=cumulativeWall+inner.wallSeconds;
    cumulativeObjectiveCalls=cumulativeObjectiveCalls+inner.callbackCounts.objective;
    cumulativeGradientCalls=cumulativeGradientCalls+inner.callbackCounts.gradient;
    row=make_cycle_row(outer,anchor,candidate,wFull,groupsFull, ...
        fAnchor,fCandidate,trueAnchor,trueCandidate,majorizerAnchor, ...
        majorizerCandidate,anchorEquality,dominanceGap,majorizerDecrease, ...
        trueDecrease,majorizerPass,projectionPass,cycleAccepted,inner, ...
        delivery,anchorDiag,candidateDiag,cumulativeWall, ...
        cumulativeObjectiveCalls,cumulativeGradientCalls,cfg,pgReference);
    rows(end+1,1)=row; %#ok<AGROW>
    save(fullfile(outputDir,sprintf('outer_cycle_%02d.mat',outer)), ...
        'inner','row','anchor','candidate','q','tangentConstant', ...
        'wFull','delivery','setup','-v7.3');
    writetable(struct2table(rows,'AsArray',true), ...
        fullfile(outputDir,'outer_mm_cycles.csv'));
    fprintf(['OUTER-MM cycle=%d accepted=%d status=%s iter=%g wall=%.1f ' ...
        'true=%.10g dTrue=%.6g dQ=%.6g primal=%.3e dual=%.3e ' ...
        'compl=%.3e exactZeroAllowed=%d active=%d BDT=%.6f\n'], ...
        outer,cycleAccepted,inner.statusMessage,inner.iterations,inner.wallSeconds, ...
        trueCandidate,trueDecrease,majorizerDecrease, ...
        inner.postKKT.primalInfeasibility,inner.postKKT.dualInfeasibility, ...
        inner.postKKT.complementarity,row.exactZeroAllowedSpots, ...
        row.activeFullSpots,row.provisionalBDT);
    if ~cycleAccepted
        stopReason='inner_or_majorizer_gate_failed'; break;
    end
    wAnchor=candidate;
end

summary=struct('completedOuterCycles',numel(rows), ...
    'acceptedOuterCycles',sum([rows.cycleAccepted]), ...
    'stopReason',stopReason,'finalTrueObjective',rows(end).trueObjectiveCandidate, ...
    'finalActiveSpots',rows(end).activeFullSpots, ...
    'finalExactZeroAllowedSpots',rows(end).exactZeroAllowedSpots, ...
    'finalActiveLayers',rows(end).activeLayers, ...
    'finalProvisionalBDT',rows(end).provisionalBDT, ...
    'cumulativeWallSeconds',rows(end).cumulativeWallSeconds, ...
    'cumulativeObjectiveCalls',rows(end).cumulativeObjectiveCalls, ...
    'cumulativeGradientCalls',rows(end).cumulativeGradientCalls, ...
    'finalIpoptPolishingRun',false,'classification',setup.classification);
save(fullfile(outputDir,'outer_mm_complete.mat'),'summary','rows','setup', ...
    'allowedMask','reducedToFull','fullToReduced','pgReference','-v7.3');
writetable(struct2table(summary,'AsArray',true), ...
    fullfile(outputDir,'outer_mm_summary.csv'));
outputs=struct('summary',summary,'rows',rows,'outputDir',outputDir);
end

function row=make_cycle_row(k,anchor,w,wFull,groupsFull,fA,fW,tA,tW,qA,qW, ...
    eqErr,domGap,dQ,dTrue,qPass,projectionPass,accepted,inner,delivery, ...
    anchorDiag,candidateDiag,cumWall,cumObj,cumGrad,cfg,pg)
row=empty_cycle_row(); row.outerCycle=k; row.trueObjectiveAnchor=tA;
row.trueObjectiveCandidate=tW; row.trueObjectiveDecrease=dTrue;
row.doseObjectiveAnchor=fA; row.doseObjectiveCandidate=fW;
row.majorizerAnchor=qA; row.majorizerCandidate=qW;
row.majorizerDecrease=dQ; row.anchorEqualityAbsoluteError=eqErr;
row.majorizerDominanceGap=domGap; row.majorizerPass=qPass;
row.innerAccepted=inner.innerAccepted; row.desiredConverged=inner.desiredConverged;
row.acceptableStabilized=inner.acceptableStabilized;
row.projectionPass=projectionPass; row.cycleAccepted=accepted;
row.ipoptStatus=inner.status; row.ipoptStatusMessage=inner.statusMessage;
row.ipoptIterations=inner.iterations; row.ipoptCpuSeconds=inner.cpuSeconds;
row.wallSeconds=inner.wallSeconds; row.cumulativeWallSeconds=cumWall;
row.primalInfeasibility=inner.postKKT.primalInfeasibility;
row.dualInfeasibility=inner.postKKT.dualInfeasibility;
row.complementarity=inner.postKKT.complementarity;
row.ipoptFinalInfPr=last_finite(inner.iterationHistory.infPr);
row.ipoptFinalInfDu=last_finite(inner.iterationHistory.infDu);
row.objectiveWindowRelativeRange=inner.stability.relativeRange;
row.objectiveCalls=inner.callbackCounts.objective;
row.gradientCalls=inner.callbackCounts.gradient;
row.constraintCalls=inner.callbackCounts.constraints;
row.jacobianCalls=inner.callbackCounts.jacobian;
row.cumulativeObjectiveCalls=cumObj; row.cumulativeGradientCalls=cumGrad;
row.exactZeroAllowedSpots=nnz(w==0);
row.exactZeroFullSpots=nnz(wFull==0);
row.activeFullSpots=nnz(wFull>cfg.spotThreshold);
row.newExactZeroAllowedSpots=nnz(anchor>0 & w==0);
row.reactivatedAllowedSpots=nnz(anchor==0 & w>0);
[row.exactZeroLayers,row.activeLayers]=layer_counts(wFull,groupsFull);
row.provisionalBDT=delivery.T_total;
row.subMMUSpots=nnz(w>0 & w<cfg.mmuThreshold);
row.minimumPositiveWeight=min_positive(w);
row.anchorConstraintViolation=anchorDiag.maxConstraintViolation;
row.candidateConstraintViolation=candidateDiag.maxConstraintViolation;
row.projectionRelativeNorm=inner.projectionRelativeNorm;
[row.nearestPGIteration,row.nearestPGWallSeconds,row.nearestPGTrueObjective, ...
    row.nearestPGActiveSpots,row.nearestPGExactZeroSpots, ...
    row.nearestPGObjectiveCalls,row.nearestPGGradientCalls]=nearest_pg(pg,cumWall);
end

function pg=build_pg_reference(pgRoot,model,reducedToFull,groups,cfg)
blockFile=fullfile(pgRoot,'mm_schedule_blocks.csv'); assert(exist(blockFile,'file')==2);
blocks=readtable(blockFile); n=height(blocks);
iteration=zeros(n,1); wallSeconds=zeros(n,1); trueObjective=zeros(n,1);
jointDiagnostic=zeros(n,1); doseObjective=zeros(n,1); activeSpots=zeros(n,1);
exactZeroSpots=zeros(n,1); objectiveCalls=zeros(n,1); gradientCalls=zeros(n,1);
for i=1:n
    iteration(i)=blocks.iteration(i); wallSeconds(i)=blocks.cumulativeSeconds(i);
    file=fullfile(pgRoot,sprintf('mm_period_01_iter_%04d.mat',iteration(i)));
    s=load(file,'optimizer'); wFull=reshape(double(s.optimizer.w),[],1);
    w=wFull(reducedToFull); [f,~,~]=evaluate_matrad_dose_model(model,w,false);
    trueObjective(i)=f+cfg.lambdaS*sum(w./(w+cfg.epsS));
    jointDiagnostic(i)=trueObjective(i)+cfg.lambdaEL*layer_penalty(w,groups);
    doseObjective(i)=f; activeSpots(i)=nnz(wFull>cfg.spotThreshold);
    exactZeroSpots(i)=nnz(wFull==0);
    backtracking=sum(s.optimizer.history.backtrackingTrials(1:iteration(i)));
    objectiveCalls(i)=1+2*iteration(i)+backtracking;
    gradientCalls(i)=2*iteration(i);
end
pg=table(iteration,wallSeconds,trueObjective,jointDiagnostic,doseObjective, ...
    activeSpots,exactZeroSpots,objectiveCalls,gradientCalls);
end

function r=empty_cycle_row()
r=struct('outerCycle',NaN,'trueObjectiveAnchor',NaN, ...
    'trueObjectiveCandidate',NaN,'trueObjectiveDecrease',NaN, ...
    'doseObjectiveAnchor',NaN,'doseObjectiveCandidate',NaN, ...
    'majorizerAnchor',NaN,'majorizerCandidate',NaN,'majorizerDecrease',NaN, ...
    'anchorEqualityAbsoluteError',NaN,'majorizerDominanceGap',NaN, ...
    'majorizerPass',false,'innerAccepted',false,'desiredConverged',false, ...
    'acceptableStabilized',false,'projectionPass',false,'cycleAccepted',false, ...
    'ipoptStatus',NaN,'ipoptStatusMessage','','ipoptIterations',NaN, ...
    'ipoptCpuSeconds',NaN,'wallSeconds',NaN,'cumulativeWallSeconds',NaN, ...
    'primalInfeasibility',NaN,'dualInfeasibility',NaN,'complementarity',NaN, ...
    'ipoptFinalInfPr',NaN,'ipoptFinalInfDu',NaN, ...
    'objectiveWindowRelativeRange',NaN,'objectiveCalls',NaN, ...
    'gradientCalls',NaN,'constraintCalls',NaN,'jacobianCalls',NaN, ...
    'cumulativeObjectiveCalls',NaN,'cumulativeGradientCalls',NaN, ...
    'exactZeroAllowedSpots',NaN,'exactZeroFullSpots',NaN, ...
    'activeFullSpots',NaN,'newExactZeroAllowedSpots',NaN, ...
    'reactivatedAllowedSpots',NaN,'exactZeroLayers',NaN,'activeLayers',NaN, ...
    'provisionalBDT',NaN,'subMMUSpots',NaN,'minimumPositiveWeight',NaN, ...
    'anchorConstraintViolation',NaN,'candidateConstraintViolation',NaN, ...
    'projectionRelativeNorm',NaN,'nearestPGIteration',NaN, ...
    'nearestPGWallSeconds',NaN,'nearestPGTrueObjective',NaN, ...
    'nearestPGActiveSpots',NaN,'nearestPGExactZeroSpots',NaN, ...
    'nearestPGObjectiveCalls',NaN,'nearestPGGradientCalls',NaN);
end

function [iter,wall,obj,active,zero,objCalls,gradCalls]=nearest_pg(pg,target)
[~,i]=min(abs(pg.wallSeconds-target)); iter=pg.iteration(i); wall=pg.wallSeconds(i);
obj=pg.trueObjective(i); active=pg.activeSpots(i); zero=pg.exactZeroSpots(i);
objCalls=pg.objectiveCalls(i); gradCalls=pg.gradientCalls(i);
end
function value=layer_penalty(w,groups)
value=0; for g=1:numel(groups), value=value+sqrt(norm(w(groups(g).ix),2)); end
end
function [nZero,nActive]=layer_counts(w,groups)
z=false(numel(groups),1); for g=1:numel(groups), z(g)=all(w(groups(g).ix)==0); end
nZero=nnz(z); nActive=numel(groups)-nZero;
end
function value=last_finite(x)
x=x(isfinite(x)); if isempty(x), value=NaN; else, value=x(end); end
end
function value=min_positive(w)
x=w(w>0); if isempty(x), value=NaN; else, value=min(x); end
end
function verify_sha256(file,expected)
[s,o]=system(sprintf('sha256sum "%s"',file)); assert(s==0);
assert(strcmpi(strtok(strtrim(o)),expected),'Frozen artifact checksum mismatch.');
end
function create_empty_output_directory(d)
if exist(d,'dir')==7
    e=dir(d); e=e(~ismember({e.name},{'.','..'}));
    assert(isempty(e),'Output directory must be empty: %s',d);
else, mkdir(d); end
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
    else, merged.(k)=override.(k); end
end
end
