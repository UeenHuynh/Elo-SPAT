function output=run_bidirectional_support_experiment(outputDir)
%RUN_BIDIRECTIONAL_SUPPORT_EXPERIMENT Stage-1 deletion/repair discovery.
narginchk(1,1);here=fileparts(mfilename('fullpath'));root=fileparts(fileparts(here));
add_paths(root,here);create_empty(outputDir);
artifactFile=fullfile(root,'results','task25_27_pipeline_20deg.mat');
stage1File=fullfile(root,'results','full_integrated_layer_20deg_job22180','alpha_0p3','alpha_0p3_final.mat');
source=load(artifactFile,'cfg','machine','dij','pln','map');
stage1=load(stage1File,'final','finalDose','finalDelivery','supportMask');
patient=load(source.cfg.patientFile,'ct','cst');cfg=merge_structs(default_config(),source.cfg);
wFullStage1=reshape(double(stage1.final.w),[],1);allowedFull=reshape(logical(stage1.supportMask),[],1);
assert(nnz(allowedFull)==136409);[dij,allowedToFull,fullToAllowed,map]= ...
    reduce_dij_to_support(source.dij,allowedFull,source.map); %#ok<ASGLU>
for i=1:numel(map),map(i).originalGlobalIx=map(i).globalIx;map(i).globalIx=i;end
[groups,~]=build_group_index(map);assert(numel(groups)==748);
[groupsFull,~]=build_group_index(source.map);w0=wFullStage1(allowedToFull);
model=build_matrad_dose_model(dij,patient.cst,source.pln,0);assert(~model.hasConstraints);
D=nominal_matrix(dij);A=state_matrix(model);n=numel(w0);zeroTol=1e-10;
reference=evaluate_plan(w0,model,D,patient,cfg,source,map,groupsFull,allowedToFull,wFullStage1);
assert(abs(reference.objective-3139.658625654965)<=1e-8*max(1,reference.objective));
assert(reference.maxConstraintViolation<=1e-6);

settings=struct('initialDeathBatch',8,'minimumDeathBatch',8,'maximumDeathBatch',64, ...
    'maximumTransitions',7,'maximumFailuresAtMinimum',2,'rescueFactor',4, ...
    'maximumCorrectionSize',256,'maximumPerLayerFraction',0.25, ...
    'objectiveToleranceFraction',0.02,'doseTolerancePercent',2, ...
    'oarAbsoluteToleranceGy',0.2,'zeroTolerance',zeroTol);
inner=struct('desiredDualTolerance',1e-3,'acceptableDualTolerance',2e-2, ...
    'acceptableObjectiveWindowTolerance',1e-6,'stabilityWindow',6, ...
    'ipopt',struct('max_iter',80,'max_cpu_time',3600,'print_level',3, ...
    'hessian_approximation','limited-memory','limited_memory_max_history',30, ...
    'acceptable_iter',6,'acceptable_tol',1e-4));

w=w0;batch=settings.initialDeathBatch;failAtMinimum=0;acceptedCount=0;
rows=repmat(empty_transition(),0,1);components=repmat(empty_component(),0,1);
mapping=repmat(empty_mapping(),0,1);fullGradientEvaluations=0;totalTimer=tic;
for k=1:settings.maximumTransitions
    incumbent=evaluate_plan(w,model,D,patient,cfg,source,map,groupsFull,allowedToFull,wFullStage1);
    [gInc,compInc]=objective_information(model,A,w);fullGradientEvaluations=fullGradientEvaluations+1+numel(compInc);
    [deathIx,predictedRemovalLoss,proposalTable]=propose_deaths( ...
        w,gInc,compInc,map,groups,batch,settings); %#ok<ASGLU>
    if isempty(deathIx),break;end
    wDeleted=w;wDeleted(deathIx)=0;
    defect=evaluate_plan(wDeleted,model,D,patient,cfg,source,map,groupsFull,allowedToFull,wFullStage1);
    rawActualLoss=defect.objective-incumbent.objective;
    [gDefect,compDefect]=objective_information(model,A,wDeleted);
    fullGradientEvaluations=fullGradientEvaluations+1+numel(compDefect);
    [correctionIx,rescueTable,curvatureScale]=price_rescue(w,wDeleted,gInc,gDefect,A, ...
        deathIx,batch,map,groups,settings,rawActualLoss);
    repairTimer=tic;repair=solve_bidirectional_repair(model,wDeleted,correctionIx,inner);
    repairSeconds=toc(repairTimer);wTrial=repair.wFull;wTrial(wTrial<=zeroTol)=0;
    repaired=evaluate_plan(wTrial,model,D,patient,cfg,source,map,groupsFull,allowedToFull,wFullStage1);
    gate=acceptance_gate(repaired,reference,settings,repair.innerAccepted);

    incumbentActive=w>zeroTol;trialActive=wTrial>zeroTol;
    successfulDeaths=deathIx(~trialActive(deathIx));returnedDeaths=deathIx(trialActive(deathIx));
    births=find(~incumbentActive&trialActive);otherDeaths=find(incumbentActive&~trialActive&~ismember((1:n)',deathIx));
    netChange=nnz(trialActive)-nnz(incumbentActive);
    row=empty_transition();row.transition=k;row.deathBatchRequested=batch;
    row.incumbentSpots=nnz(incumbentActive);row.incumbentLayers=count_layers(w,groups,zeroTol);
    row.proposedDeaths=numel(deathIx);row.successfulDeaths=numel(successfulDeaths);
    row.returnedDeaths=numel(returnedDeaths);row.rescueBirths=numel(births);
    row.otherDeaths=numel(otherDeaths);row.netSupportChange=netChange;
    correctionMask=false(n,1);correctionMask(correctionIx)=true;
    row.repairSetSize=numel(correctionIx);row.restrictedSupportSize=nnz((wDeleted>zeroTol)|correctionMask);
    row.objectiveBefore=incumbent.objective;row.objectiveAfterDeletion=defect.objective;
    row.objectiveAfterRepair=repaired.objective;row.predictedDeletionLoss=predictedRemovalLoss;
    row.actualDeletionLoss=rawActualLoss;row.actualRepairedLoss=repaired.objective-incumbent.objective;
    row.predictionRatio=rawActualLoss/max(abs(predictedRemovalLoss),eps);
    row.maxConstraintViolation=repaired.maxConstraintViolation;row.worstAdversePercent=gate.worstAdverse;
    row.objectivePass=gate.objectivePass;row.dosePass=gate.dosePass;
    row.solverAccepted=repair.innerAccepted;row.accepted=gate.pass;
    row.failureClass=gate.failureClass;row.ipoptStatus=repair.statusMessage;
    row.ipoptIterations=repair.iterations;row.repairWallSeconds=repairSeconds;
    row.objectiveEvaluations=repair.callbackCounts.objective;
    row.gradientEvaluations=repair.callbackCounts.gradient;row.fullGradientEvaluations=fullGradientEvaluations;
    row.bdtBefore=incumbent.delivery.T_total;row.bdtAfterRepair=repaired.delivery.T_total;
    row.ptvD95=repaired.dose.summary.targetD95Min;row.ptvD98=repaired.dose.summary.targetD98Min;
    row.oarDmax=repaired.dose.summary.oarDmaxMax;row.affectedLayers=count_affected_layers(deathIx,correctionIx,groups,n);
    row.deathGlobalIx=index_list(allowedToFull(deathIx));
    row.successfulDeathGlobalIx=index_list(allowedToFull(successfulDeaths));
    row.returnedDeathGlobalIx=index_list(allowedToFull(returnedDeaths));
    row.birthGlobalIx=index_list(allowedToFull(births));
    row.wallSeconds=toc(totalTimer);rows(end+1,1)=row; %#ok<AGROW>

    components=[components;component_rows(k,'defect',compInc,compDefect); ...
        component_rows(k,'repair',compInc,repaired.components)]; %#ok<AGROW>
    mapping=[mapping;mapping_rows(k,deathIx,correctionIx,w,wTrial,map,groups, ...
        allowedToFull,gDefect,compDefect)]; %#ok<AGROW>
    writetable(struct2table(rows,'AsArray',true),fullfile(outputDir,'bidirectional_transitions.csv'));
    writetable(struct2table(components,'AsArray',true),fullfile(outputDir,'bidirectional_components.csv'));
    writetable(struct2table(mapping,'AsArray',true),fullfile(outputDir,'death_rescue_mapping.csv'));
    checkpoint=struct('transition',k,'wIncumbent',w,'wDeleted',wDeleted,'wTrial',wTrial, ...
        'deathIx',deathIx,'correctionIx',correctionIx,'repair',repair,'gate',gate, ...
        'proposalTable',proposalTable,'rescueTable',rescueTable,'curvatureScale',curvatureScale,'row',row); %#ok<NASGU>
    save(fullfile(outputDir,sprintf('transition_%02d.mat',k)),'checkpoint','-v7.3');
    fprintf(['BIDIR k=%d batch=%d kill=%d return=%d birth=%d net=%d objRaw=%.4f ' ...
        'objRepair=%.4f accepted=%d solver=%s\n'],k,batch,numel(successfulDeaths), ...
        numel(returnedDeaths),numel(births),netChange,rawActualLoss, ...
        row.actualRepairedLoss,gate.pass,repair.statusMessage);
    if gate.pass&&netChange<0
        w=wTrial;acceptedCount=acceptedCount+1;failAtMinimum=0;batch=min(2*batch,settings.maximumDeathBatch);
    else
        if batch==settings.minimumDeathBatch,failAtMinimum=failAtMinimum+1;end
        batch=max(settings.minimumDeathBatch,floor(batch/2));
    end
    if failAtMinimum>=settings.maximumFailuresAtMinimum,break;end
end

[final,finalSummary]=finalize_result(w,w0,reference,rows,model,D,patient,cfg,source, ...
    groupsFull,allowedToFull,wFullStage1,allowedFull,settings,acceptedCount);
status=classify_status(rows,finalSummary,settings);
output=struct('status',status,'settings',settings,'reference',reference,'rows',rows, ...
    'final',final,'finalSummary',finalSummary,'acceptedTransitions',acceptedCount, ...
    'wallSeconds',toc(totalTimer),'stage1File',stage1File,'artifactFile',artifactFile);
save(fullfile(outputDir,'bidirectional_support_complete.mat'),'output','-v7.3');
write_report(outputDir,output);fprintf('FINAL_STATUS: %s\n',status);
end

function [g,comps]=objective_information(model,A,w)
[~,g]=evaluate_matrad_dose_model(model,w,true);x=A*w;nVox=numel(x);
comps=repmat(struct('label','','value',NaN,'doseGradient',zeros(nVox,1), ...
    'weightGradient',zeros(numel(w),1)),0,1);
for i=1:size(model.cst,1)
    vox=reshape(double(model.cst{i,4}{1}),[],1);if isempty(vox),continue;end
    for j=1:numel(model.cst{i,6})
        o=model.cst{i,6}{j};if ~isa(o,'DoseObjectives.matRad_DoseObjective'),continue;end
        assert(strcmp(o.robustness,'none'));
        o=model.optiProb.BP.setBiologicalDosePrescriptions(o,model.cst{i,5}.alphaX,model.cst{i,5}.betaX);
        c=struct('label',sprintf('%s:%s',char(string(model.cst{i,2})),class(o)), ...
            'value',double(o.penalty*o.computeDoseObjectiveFunction(x(vox))), ...
            'doseGradient',zeros(nVox,1),'weightGradient',zeros(numel(w),1));
        c.doseGradient(vox)=o.penalty*o.computeDoseObjectiveGradient(x(vox));
        c.weightGradient=A'*c.doseGradient;comps(end+1,1)=c; %#ok<AGROW>
    end
end
assert(norm(sum(reshape([comps.weightGradient],numel(w),[]),2)-g)<=1e-8*max(1,norm(g)));
end

function [deathIx,predictedLoss,tbl]=propose_deaths(w,g,comps,map,groups,batch,settings)
n=numel(w);risk=zeros(n,1);dominant=zeros(n,1);dominantValue=zeros(n,1);
for c=1:numel(comps)
    response=abs(w.*comps(c).weightGradient)/max(1,abs(comps(c).value));
    risk=risk+response.^2;take=response>dominantValue;dominant(take)=c;dominantValue(take)=response(take);
end
risk=sqrt(risk);active=w>settings.zeroTolerance;risk(~active)=Inf;
[~,order]=sort(risk,'ascend');groupOf=zeros(n,1);for q=1:numel(groups),groupOf(groups(q).ix)=q;end
cp=[map.cpIndex]';layerCap=max(2,ceil(settings.maximumPerLayerFraction*batch));
cpCap=max(2,ceil(batch/4));layerCount=zeros(numel(groups),1);cpCount=zeros(max(cp),1);deathIx=zeros(0,1);
for q=1:numel(order)
    i=order(q);if ~isfinite(risk(i)),break;end;gi=groupOf(i);ci=cp(i);
    if layerCount(gi)>=layerCap||cpCount(ci)>=cpCap,continue;end
    deathIx(end+1,1)=i;layerCount(gi)=layerCount(gi)+1;cpCount(ci)=cpCount(ci)+1; %#ok<AGROW>
    if numel(deathIx)>=batch,break;end
end
predictedLoss=-g(deathIx)'*w(deathIx);
labels=cell(numel(deathIx),1);for i=1:numel(deathIx),if dominant(deathIx(i))>0,labels{i}=comps(dominant(deathIx(i))).label;end,end
tbl=table(deathIx,reshape([map(deathIx).originalGlobalIx],[],1),w(deathIx),risk(deathIx), ...
    g(deathIx),labels,'VariableNames',{'allowedIx','globalIx','weight','clinicalRisk','gradient','dominantComponent'});
end

function [ix,tbl,gamma]=price_rescue(wInc,wDef,gInc,gDef,A,deathIx,batch,map,groups,settings,rawLoss)
colNorm2=reshape(full(sum(A.^2,1)),[],1);den=sum(colNorm2(deathIx).*wInc(deathIx).^2);
gamma=max(1e-12,2*max(rawLoss+gInc(deathIx)'*wInc(deathIx),1e-12)/max(den,eps));
L=max(gamma*colNorm2,1e-12);active=wDef>settings.zeroTolerance;score=zeros(size(wDef));step=zeros(size(wDef));
down=active&gDef>0;step(down)=min(wDef(down),gDef(down)./L(down));
score(down)=gDef(down).*step(down)-0.5*L(down).*step(down).^2;
up=gDef<0;step(up)=-gDef(up)./L(up);score(up)=0.5*gDef(up).^2./L(up);
target=min(settings.maximumCorrectionSize,max(2*batch,settings.rescueFactor*batch));
[~,order]=sort(score,'descend');order=order(score(order)>0);ix=order(1:min(target,numel(order)));
if isempty(ix),[~,order]=sort(abs(gDef),'descend');ix=order(1:min(target,numel(order)));end
groupOf=zeros(numel(wDef),1);for q=1:numel(groups),groupOf(groups(q).ix)=q;end
tbl=table(ix,reshape([map(ix).originalGlobalIx],[],1),score(ix),gDef(ix),step(ix), ...
    active(ix),ismember(ix,deathIx),groupOf(ix),'VariableNames', ...
    {'allowedIx','globalIx','rescueScore','defectGradient','predictedStep','wasActive','wasDeleted','layerId'});
end

function p=evaluate_plan(w,model,D,patient,cfg,source,mapAllowed,groupsFull,allowedToFull,wFullTemplate)
[f,~,d]=evaluate_matrad_dose_model(model,w,false);components=component_values(model,state_matrix(model)*w);
z=zeros(size(wFullTemplate));z(allowedToFull)=w;
p=struct('objective',f,'components',components,'maxConstraintViolation',d.maxConstraintViolation, ...
    'dose',compute_dvh_metrics(full(D*w),patient.ct,patient.cst,cfg), ...
    'delivery',compute_delivery_metrics(z,source.map,groupsFull,cfg,source.machine), ...
    'activeSpots',nnz(w>1e-10)); %#ok<NASGU>
end

function comps=component_values(model,x)
comps=repmat(struct('label','','value',NaN),0,1);
for i=1:size(model.cst,1),vox=reshape(double(model.cst{i,4}{1}),[],1);if isempty(vox),continue;end
    for j=1:numel(model.cst{i,6}),o=model.cst{i,6}{j};if ~isa(o,'DoseObjectives.matRad_DoseObjective'),continue;end
        o=model.optiProb.BP.setBiologicalDosePrescriptions(o,model.cst{i,5}.alphaX,model.cst{i,5}.betaX);
        comps(end+1,1)=struct('label',sprintf('%s:%s',char(string(model.cst{i,2})),class(o)), ...
            'value',double(o.penalty*o.computeDoseObjectiveFunction(x(vox)))); %#ok<AGROW>
    end
end
end

function A=state_matrix(model)
D=model.dij.physicalDose;if iscell(D),D=D{1};end
if isa(model.optiProb.BP,'matRad_DoseProjection'),A=D;
elseif isa(model.optiProb.BP,'matRad_ConstantRBEProjection'),A=double(model.dij.RBE)*D;
else,error('Bidirectional prototype requires a linear exact clinical state.');end
end

function gate=acceptance_gate(p,reference,s,solverAccepted)
[dosePass,worst]=dose_quality(p.dose,reference.dose,s.doseTolerancePercent,s.oarAbsoluteToleranceGy);
objectivePass=p.objective<=(1+s.objectiveToleranceFraction)*reference.objective;
constraintPass=p.maxConstraintViolation<=1e-6;
if ~solverAccepted,failure='numerical_solver_failure';
elseif ~(dosePass&&objectivePass&&constraintPass),failure='support_clinically_insufficient';
else,failure='none';end
gate=struct('pass',solverAccepted&&dosePass&&objectivePass&&constraintPass, ...
    'dosePass',dosePass,'objectivePass',objectivePass,'constraintPass',constraintPass, ...
    'worstAdverse',worst,'failureClass',failure);
end

function [pass,worst]=dose_quality(c,r,pctTol,absTol)
pass=true;worst=0;
for i=1:numel(c.structures),a=r.structures(i);b=c.structures(i);
    if a.isTarget
        adverse=max([0,-pct(b.D95,a.D95),-pct(b.D98,a.D98)]);ok=adverse<=pctTol;
    elseif strcmpi(a.type,'OAR')&&~strcmpi(a.name,'BODY')&&~strcmpi(a.name,'EXTERNAL')
        dm=b.Dmean-a.Dmean;dx=b.Dmax-a.Dmax;
        ok=dm<=max(absTol,pctTol/100*max(abs(a.Dmean),1))&& ...
            dx<=max(absTol,pctTol/100*max(abs(a.Dmax),1));
        adverse=max([0,pct(b.Dmean,a.Dmean),pct(b.Dmax,a.Dmax)]);
    else,adverse=0;ok=true;end
    pass=pass&&ok;worst=max(worst,adverse);
end
end

function rows=component_rows(k,phase,before,after)
rows=repmat(empty_component(),numel(before),1);
for i=1:numel(before)
    rows(i).transition=k;rows(i).phase=phase;rows(i).component=before(i).label;
    rows(i).incumbentValue=before(i).value;rows(i).trialValue=after(i).value;
    rows(i).change=after(i).value-before(i).value;
end
end

function rows=mapping_rows(k,deaths,rescue,wBefore,wAfter,map,groups,toFull,gDefect,components)
rows=repmat(empty_mapping(),numel(rescue),1);groupOf=zeros(numel(wBefore),1);
for q=1:numel(groups),groupOf(groups(q).ix)=q;end
deathCp=unique([map(deaths).cpIndex]);deathLayers=unique(groupOf(deaths));
for i=1:numel(rescue),r=rescue(i);rows(i).transition=k;rows(i).rescueAllowedIx=r;
    rows(i).rescueGlobalIx=toFull(r);rows(i).wasDeleted=ismember(r,deaths);
    rows(i).wasInactiveBefore=wBefore(r)<=1e-10;rows(i).weightBefore=wBefore(r);
    rows(i).weightAfter=wAfter(r);rows(i).becameActive=wBefore(r)<=1e-10&&wAfter(r)>1e-10;
    rows(i).defectGradient=gDefect(r);[rows(i).dominantComponentGradient,q]= ...
        max(abs(arrayfun(@(c)c.weightGradient(r)/max(1,abs(c.value)),components)));
    rows(i).dominantDefectComponent=components(q).label;rows(i).cpIndex=map(r).cpIndex;
    rows(i).layerId=groupOf(r);rows(i).sameDeathCP=ismember(map(r).cpIndex,deathCp);
    rows(i).sameDeathLayer=ismember(groupOf(r),deathLayers);
    rows(i).deathGlobalIx=index_list(toFull(deaths));
end
end

function [final,summary]=finalize_result(w,w0,reference,rows,model,D,patient,cfg,source, ...
    groupsFull,toFull,wFullTemplate,allowedFull,settings,acceptedCount)
netReduction=nnz(w0>settings.zeroTolerance)-nnz(w>settings.zeroTolerance);
runPolish=acceptedCount>=2&&netReduction>=8;final=[];
summary=struct('polishRun',runPolish,'netReductionBeforePolish',netReduction, ...
    'supportPreserved',false,'solverAccepted',false,'planAccepted',false, ...
    'objective',NaN,'activeSpots',nnz(w>settings.zeroTolerance),'activeLayers',NaN, ...
    'bdt',NaN,'worstAdverse',NaN,'failureClass','not_run');
if ~runPolish,return;end
supportAllowed=w>settings.zeroTolerance;supportFull=false(size(wFullTemplate));supportFull(toFull(supportAllowed))=true;
wInit=zeros(size(wFullTemplate));wInit(toFull(supportAllowed))=w(supportAllowed);
assert(~any(supportFull(~allowedFull)));
timer=tic;final=run_final_reoptimization(source.dij,patient.cst,source.pln,wInit,cfg,supportFull,[]);
summary.polishWallSeconds=toc(timer);tol=1e-8*max(1,max(abs(final.w)));
summary.supportPreserved=~any(abs(final.w(~supportFull))>tol);
solver=compute_solver_metrics(final.resultGUI,struct(),final.optimizer);
summary.solverAccepted=solver.converged||solver.acceptableButNotConverged||contains(lower(solver.statusMessage),'acceptable');
wAllowed=final.w(toFull);plan=evaluate_plan(wAllowed,model,D,patient,cfg,source,[],groupsFull,toFull,wFullTemplate);
gate=acceptance_gate(plan,reference,settings,summary.solverAccepted&&summary.supportPreserved);
summary.planAccepted=gate.pass;summary.objective=plan.objective;summary.activeSpots=nnz(wAllowed>settings.zeroTolerance);
summary.activeLayers=count_layers(wAllowed,build_group_index_local(source.map,toFull),settings.zeroTolerance);
summary.bdt=plan.delivery.T_total;summary.worstAdverse=gate.worstAdverse;summary.failureClass=gate.failureClass;
summary.solverStatus=solver.statusMessage;summary.netReductionAfterPolish=nnz(w0>settings.zeroTolerance)-summary.activeSpots;
end

function groups=build_group_index_local(fullMap,toFull)
m=fullMap(toFull);for i=1:numel(m),m(i).globalIx=i;end;groups=build_group_index(m);
end

function status=classify_status(rows,finalSummary,settings)
if isempty(rows),status='INCONCLUSIVE';return;end
accepted=[rows.accepted];solverAccepted=[rows.solverAccepted];
reduction=-sum([rows(accepted).netSupportChange]);
if sum(accepted)>=2&&reduction>=16&&finalSummary.polishRun&&finalSummary.planAccepted
    status='BIDIRECTIONAL_SUPPORT_MECHANISM_SUPPORTED';
elseif mean(solverAccepted)<0.5
    status='INCONCLUSIVE';
elseif sum(accepted)<2||reduction<8|| ...
        (finalSummary.polishRun&&~finalSummary.planAccepted)
    status='BIDIRECTIONAL_SUPPORT_MECHANISM_NOT_SUPPORTED';
else
    status='INCONCLUSIVE';
end
end

function write_report(out,o)
fid=fopen(fullfile(out,'BIDIRECTIONAL_SUPPORT_REPORT.md'),'w');assert(fid>0);c=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'# Bidirectional Stage-2 support-compression experiment\n\nStatus: **%s**\n\n',o.status);
fprintf(fid,'Stage-1 objective %.12g; BDT %.6f s; allowed support 136409 spots / 748 layers.\n\n', ...
    o.reference.objective,o.reference.delivery.T_total);
fprintf(fid,'| k | deaths | killed | returned | births | net | objective raw loss | repaired loss | accepted | failure | repair vars | BDT |\n');
fprintf(fid,'|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---:|---:|\n');
for i=1:numel(o.rows),r=o.rows(i);fprintf(fid,'| %d | %d | %d | %d | %d | %d | %.6g | %.6g | %d | %s | %d | %.3f |\n', ...
        r.transition,r.proposedDeaths,r.successfulDeaths,r.returnedDeaths,r.rescueBirths, ...
        r.netSupportChange,r.actualDeletionLoss,r.actualRepairedLoss,r.accepted, ...
        r.failureClass,r.repairSetSize,r.bdtAfterRepair);end
fprintf(fid,'\nAccepted transitions: %d. Total wall time %.1f s.\n\n',o.acceptedTransitions,o.wallSeconds);
fprintf(fid,'Final polish run=%d, accepted=%d, active spots=%d, objective=%.9g, BDT=%.3f s.\n\n', ...
    o.finalSummary.polishRun,o.finalSummary.planAccepted,o.finalSummary.activeSpots, ...
    o.finalSummary.objective,o.finalSummary.bdt);
fprintf(fid,'Comparator job 22292: 35204 spots, objective about 4253.25 (+35%%), BDT about 4752.56 s; it is not used as initialization.\n\n');
fprintf(fid,'Artifacts: `bidirectional_transitions.csv`, `bidirectional_components.csv`, `death_rescue_mapping.csv`, and per-transition MAT checkpoints.\n');
end

function n=count_layers(w,groups,tol),n=0;for q=1:numel(groups),n=n+any(w(groups(q).ix)>tol);end,end
function n=count_affected_layers(a,b,groups,N),lookup=zeros(N,1);for q=1:numel(groups),lookup(groups(q).ix)=q;end;n=numel(unique(lookup(unique([a;b]))));end
function D=nominal_matrix(dij),D=dij.physicalDose;if iscell(D),D=D{1};end,end
function x=pct(a,b),x=100*(double(a)-double(b))/max(abs(double(b)),eps);end
function s=index_list(ix),if isempty(ix),s='';else,s=strjoin(arrayfun(@(x)sprintf('%d',x),reshape(ix,1,[]),'UniformOutput',false),';');end,end

function r=empty_transition()
r=struct('transition',0,'deathBatchRequested',0,'incumbentSpots',0,'incumbentLayers',0, ...
    'proposedDeaths',0,'successfulDeaths',0,'returnedDeaths',0,'rescueBirths',0, ...
    'otherDeaths',0,'netSupportChange',0,'repairSetSize',0,'restrictedSupportSize',0, ...
    'objectiveBefore',NaN,'objectiveAfterDeletion',NaN,'objectiveAfterRepair',NaN, ...
    'predictedDeletionLoss',NaN,'actualDeletionLoss',NaN,'actualRepairedLoss',NaN, ...
    'predictionRatio',NaN,'maxConstraintViolation',NaN,'worstAdversePercent',NaN, ...
    'objectivePass',false,'dosePass',false,'solverAccepted',false,'accepted',false, ...
    'failureClass','','ipoptStatus','','ipoptIterations',NaN,'repairWallSeconds',NaN, ...
    'objectiveEvaluations',0,'gradientEvaluations',0,'fullGradientEvaluations',0, ...
    'bdtBefore',NaN,'bdtAfterRepair',NaN,'ptvD95',NaN,'ptvD98',NaN,'oarDmax',NaN, ...
    'affectedLayers',0,'deathGlobalIx','','successfulDeathGlobalIx','', ...
    'returnedDeathGlobalIx','','birthGlobalIx','','wallSeconds',NaN);
end
function r=empty_component(),r=struct('transition',0,'phase','','component','','incumbentValue',NaN,'trialValue',NaN,'change',NaN);end
function r=empty_mapping(),r=struct('transition',0,'deathGlobalIx','','rescueAllowedIx',0,'rescueGlobalIx',0,'wasDeleted',false,'wasInactiveBefore',false,'weightBefore',NaN,'weightAfter',NaN,'becameActive',false,'defectGradient',NaN,'dominantDefectComponent','','dominantComponentGradient',NaN,'cpIndex',0,'layerId',0,'sameDeathCP',false,'sameDeathLayer',false);end
function add_paths(root,here),f={'config','planning','mapping','objective','delivery','reduction','evaluation','integrated','robustness'};for i=1:numel(f),addpath(fullfile(root,f{i}));end,addpath(here);end
function create_empty(d),if exist(d,'dir')==7,e=dir(d);e=e(~ismember({e.name},{'.','..'}));assert(isempty(e));else,mkdir(d);end,end
function m=merge_structs(a,b),m=a;f=fieldnames(b);for i=1:numel(f),k=f{i};if isfield(m,k)&&isstruct(m.(k))&&isstruct(b.(k)),m.(k)=merge_structs(m.(k),b.(k));else,m.(k)=b.(k);end,end,end
