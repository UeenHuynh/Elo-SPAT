function outputs = run_stage2_spot_continuation_20deg(outputDir, alphaS)
% Continue Stage-2 exactly from a job-22208 iteration-600 checkpoint.
% Final IPOPT is delegated to the validated Stage-2 runner only after Gate C.

narginchk(2,2);
root=fileparts(fileparts(mfilename('fullpath'))); add_module_paths(root);
assert(any(abs(alphaS-[0.01 0.03])<=1e-12), ...
    'Only alphaS=0.01 and 0.03 are authorized for this continuation.');
if exist(outputDir,'dir')~=7, mkdir(outputDir); end

artifactFile=fullfile(root,'results','task25_27_pipeline_20deg.mat');
frozenFile=fullfile(root,'results','full_integrated_layer_20deg_job22180', ...
    'alpha_0p3','alpha_0p3_final.mat');
sourceRoot=fullfile(root,'results','full_integrated_layer_spot_20deg_job22208');
token=alpha_token(alphaS);
seedFile=fullfile(sourceRoot,['alphaS_' token], ...
    sprintf('alphaS_%s_iter_0600.mat',token));
assert(exist(seedFile,'file')==2,'Missing iteration-600 seed: %s',seedFile);

source=load(artifactFile,'cfg','machine','dij','pln','map');
layer=load(frozenFile,'record','final');
seed=load(seedFile,'optimizer','bdtHistory','stabilization','cfgTrial','alphaS');
assert(abs(seed.alphaS-alphaS)<=1e-12);
assert(seed.optimizer.iterations==600,'Seed must be iteration 600.');
assert(abs(layer.record.alpha-0.3)<=1e-12 && layer.record.planAccepted);
cfgTrial=seed.cfgTrial;
assert(strcmpi(cfgTrial.integrated.componentStage,'spot'));
assert(cfgTrial.lambdaS>0 && cfgTrial.lambdaEL==layer.record.lambdaEL);
assert(cfgTrial.lambdaB==0 && cfgTrial.lambdaES==0 && cfgTrial.lambdaT==0);
cfgTrial.integrated.maxIterations=50;
cfgTrial.integrated.stationarityStoppingEnabled=false;
cfgTrial.integrated.verbose=true;

patient=load(source.cfg.patientFile,'ct','cst');
[groups,cpGroups]=build_group_index(source.map);
doseModel=build_matrad_dose_model(source.dij,patient.cst,source.pln, ...
    cfgTrial.integrated.constraintPenalty);
context=build_integrated_context(groups,cpGroups,source.map,source.machine,cfgTrial);
verify_seed(seed.optimizer,doseModel,context,cfgTrial);

rule=default_stabilization_rule();
blockSize=50; maxAdditionalIterations=600; maxTotal=1200;
stagnationBlocks=3; stagnationMaximumNewZeros=1;
stagnationObjectiveTolerance=1e-6;
stagnation=struct('blocks',stagnationBlocks, ...
    'maximumNewZerosPerBlock',stagnationMaximumNewZeros, ...
    'relativeObjectiveWindowTolerance',stagnationObjectiveTolerance, ...
    'maximumStepRatio',rule.maximumStepRatio);
setup=struct('alphaS',alphaS,'alphaEL',layer.record.alpha, ...
    'lambdaS',cfgTrial.lambdaS,'lambdaEL',cfgTrial.lambdaEL, ...
    'seedFile',seedFile,'seedIteration',600,'blockSize',blockSize, ...
    'maxAdditionalIterations',maxAdditionalIterations,'maxTotalIterations',maxTotal, ...
    'rule',rule,'stagnation',stagnation,'frozenLayerFile',frozenFile, ...
    'penalties',struct('lambdaB',0,'lambdaES',0,'lambdaT',0));
save(fullfile(outputDir,'continuation_setup.mat'),'setup','-v7.3');

state=seed.optimizer; bdtHistory=seed.bdtHistory;
blockRows=repmat(empty_block(),0,1); stopReason='additional_budget_exhausted';
while state.iterations<maxTotal
    blockStart=state.iterations;
    startMask=(state.w==0);
    startObjective=state.objective;
    state=run_integrated_optimizer(doseModel,context,state.w,cfgTrial,state);
    delivery=compute_delivery_metrics(state.w,source.map,groups,cfgTrial,source.machine);
    bdtHistory.iteration(end+1,1)=state.iterations;
    bdtHistory.value(end+1,1)=delivery.T_total;
    stabilization=evaluate_numerical_stabilization(state.history,bdtHistory,rule);
    endMask=(state.w==0);
    row=make_block(alphaS,blockStart,state,startMask,endMask,startObjective, ...
        delivery.T_total,stabilization,rule);
    blockRows(end+1,1)=row; %#ok<AGROW>
    exactZeroSpotMask=endMask; %#ok<NASGU>
    optimizer=state; %#ok<NASGU>
    supportIdentity=row.supportIdentity; %#ok<NASGU>
    save(fullfile(outputDir,sprintf('alphaS_%s_iter_%04d.mat',token,state.iterations)), ...
        'optimizer','bdtHistory','stabilization','cfgTrial','alphaS', ...
        'exactZeroSpotMask','supportIdentity','seedFile','frozenFile','-v7.3');
    writetable(struct2table(blockRows,'AsArray',true), ...
        fullfile(outputDir,'continuation_blocks.csv'));
    fprintf(['BLOCK alphaS=%g %d->%d newZero=%d active=%d hash=%s ' ...
        'objWindow=%.6e step=[%.3e,%.3e] stationarity=%.6e BDT=%.6f ' ...
        'GateC=%d reason=%s\n'],alphaS,row.blockStart,row.blockEnd, ...
        row.newExactZeroSpots,row.activeSpots,row.supportIdentity, ...
        row.objectiveWindowChange,row.minimumStep,row.maximumStep, ...
        row.normalizedStationarity,row.provisionalBDT,row.gateCPass,row.gateCReason);

    if stabilization.numerically_stabilized
        stopReason='gate_c_satisfied'; break;
    end
    if clear_stagnation(blockRows,stagnation)
        stopReason='clear_numerical_stagnation'; break;
    end
    if ~strcmp(state.status,'maximum_iterations')
        stopReason=state.stoppingReason; break;
    end
end

summary=struct('alphaS',alphaS,'lambdaS',cfgTrial.lambdaS, ...
    'startIteration',600,'finalIteration',state.iterations, ...
    'additionalIterations',state.iterations-600,'stopReason',stopReason, ...
    'gateCPass',stabilization.numerically_stabilized, ...
    'gateCReason',stabilization.reason,'activeSpots',nnz(state.w>cfgTrial.spotThreshold), ...
    'exactZeroSpots',nnz(state.w==0),'activeLayers',state.history.exactActiveLayers(end), ...
    'exactZeroLayers',state.history.exactZeroLayers(end), ...
    'objective',state.objective,'normalizedStationarity',state.finalNormalizedStationarity, ...
    'provisionalBDT',bdtHistory.value(end),'finalIpoptRun',false);
save(fullfile(outputDir,'continuation_complete.mat'),'summary','blockRows', ...
    'optimizer','bdtHistory','stabilization','cfgTrial','setup','-v7.3');
writetable(struct2table(summary,'AsArray',true), ...
    fullfile(outputDir,'continuation_summary.csv'));

% Reuse the validated final-IPOPT/acceptance path only after Gate C. The
% checkpoint naming and variables above are intentionally compatible with it.
if stabilization.numerically_stabilized
    fprintf('Gate C passed; invoking support-constrained final IPOPT.\n');
    finalOutputs=run_full_integrated_layer_spot_20deg(outputDir,alphaS); %#ok<NASGU>
    summary.finalIpoptRun=true;
    save(fullfile(outputDir,'continuation_complete.mat'),'summary','blockRows', ...
        'optimizer','bdtHistory','stabilization','cfgTrial','setup', ...
        'finalOutputs','-v7.3');
    writetable(struct2table(summary,'AsArray',true), ...
        fullfile(outputDir,'continuation_summary.csv'));
end
outputs=struct('summary',summary,'blocks',blockRows,'outputDir',outputDir);
end

function row=make_block(alphaS,startIter,state,startMask,endMask,startObjective,bdt,stab,rule)
h=state.history; ix=h.iteration>startIter; steps=h.stepSize(ix);
row=empty_block(); row.alphaS=alphaS; row.blockStart=startIter;
row.blockEnd=state.iterations; row.newExactZeroSpots=nnz(~startMask & endMask);
row.reactivatedSpots=nnz(startMask & ~endMask);
row.activeSpots=nnz(state.w>0); row.exactZeroSpots=nnz(endMask);
row.supportIdenticalToPrevious=isequal(startMask,endMask);
row.supportIdentity=support_hash(endMask);
row.objectiveStart=startObjective; row.objectiveEnd=state.objective;
row.objectiveWindowChange=stab.relativeObjectiveWindow;
row.minimumStep=min(steps); row.maximumStep=max(steps);
row.stepRatio=row.maximumStep/max(row.minimumStep,realmin);
row.normalizedStationarity=state.finalNormalizedStationarity;
row.provisionalBDT=bdt; row.gateCPass=stab.numerically_stabilized;
row.gateCReason=stab.reason; row.objectiveStable=stab.objectiveStable;
row.spotSupportStable=stab.spotSupportStable;
row.noObjectiveIncrease=stab.maximumObjectiveIncrease<= ...
    rule.objectiveIncreaseTolerance*max(1,max(abs(h.objective(ix))));
end

function tf=clear_stagnation(rows,rule)
tf=false; if numel(rows)<rule.blocks, return; end
q=rows(end-rule.blocks+1:end);
tf=all([q.newExactZeroSpots]<=rule.maximumNewZerosPerBlock) && ...
    all([q.reactivatedSpots]==0) && ...
    all([q.objectiveWindowChange]<=rule.relativeObjectiveWindowTolerance) && ...
    all([q.stepRatio]<=rule.maximumStepRatio) && all([q.noObjectiveIncrease]);
end

function row=empty_block()
row=struct('alphaS',NaN,'blockStart',NaN,'blockEnd',NaN, ...
    'newExactZeroSpots',NaN,'reactivatedSpots',NaN,'activeSpots',NaN, ...
    'exactZeroSpots',NaN,'supportIdenticalToPrevious',false, ...
    'supportIdentity','','objectiveStart',NaN,'objectiveEnd',NaN, ...
    'objectiveWindowChange',NaN,'minimumStep',NaN,'maximumStep',NaN, ...
    'stepRatio',NaN,'normalizedStationarity',NaN,'provisionalBDT',NaN, ...
    'gateCPass',false,'gateCReason','','objectiveStable',false, ...
    'spotSupportStable',false,'noObjectiveIncrease',false);
end

function verify_seed(state,doseModel,context,cfg)
assert(numel(state.w)==context.nWeights && numel(state.history.iteration)==600);
assert(state.history.iteration(end)==600 && state.acceptedStep>0);
[objective,~]=evaluate_integrated_objective(state.w,doseModel,context,cfg);
err=abs(objective-state.objective)/max(1,abs(state.objective));
assert(err<=1e-11,'Seed objective mismatch: %.3e',err);
end

function value=support_hash(mask)
md=java.security.MessageDigest.getInstance('SHA-256');
md.update(uint8(mask(:))); bytes=typecast(md.digest(),'uint8');
value=lower(reshape(dec2hex(bytes,2).',1,[]));
end

function token=alpha_token(alpha), token=strrep(sprintf('%.4g',alpha),'.','p'); end
function add_module_paths(root), f={'config','planning','mapping','objective','delivery','reduction','robustness','evaluation','integrated'}; for i=1:numel(f), addpath(fullfile(root,f{i})); end, end
