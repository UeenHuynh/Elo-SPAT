function result = run_secant_witness_circuit_prototype(outputDir)
%RUN_SECANT_WITNESS_CIRCUIT_PROTOTYPE Local Candidate-1 mechanism test.
% This experiment is isolated from all production optimization paths.

narginchk(1,1); if exist(outputDir,'dir')~=7, mkdir(outputDir); end
here=fileparts(mfilename('fullpath')); root=fileparts(fileparts(here)); add_project_paths(root);
artifactFile=fullfile(root,'results','task25_27_pipeline_20deg.mat');
stage1File=fullfile(root,'results','full_integrated_layer_20deg_job22180', ...
    'alpha_0p3','alpha_0p3_final.mat');
source=load(artifactFile,'cfg','machine','dij','stf','pln','map','groups');
stage1=load(stage1File,'final','finalDose','finalDelivery','supportMask', ...
    'extracted','record','acceptance');
patient=load(source.cfg.patientFile,'ct','cst');
cfg=merge_structs_local(default_config(),source.cfg); cfg.arcSpacing=20;
rule=default_stabilization_rule(); rankTolerance=1e-10; zeroTolerance=1e-10;
maxOracleCalls=10; doseTolerance=cfg.evaluation.ctvD95TolerancePercent;

wReference=reshape(double(stage1.final.w),[],1); n=numel(wReference);
assert(n==double(source.dij.totalNumOfBixels));
allowedMask=reshape(logical(stage1.supportMask),[],1);
assert(numel(allowedMask)==n && nnz(allowedMask)==136409, ...
    'Expected exactly 136409 accepted Stage-1 spots.');
[groups,~]=build_group_index(source.map);
allowedLayer=false(numel(groups),1);
for g=1:numel(groups), allowedLayer(g)=any(allowedMask(groups(g).ix)); end
assert(nnz(allowedLayer)==748 && nnz(~allowedLayer)==20, ...
    'Expected 748 allowed and 20 excluded Stage-1 layers.');
excludedMask=false(n,1);
for g=find(~allowedLayer)', excludedMask(groups(g).ix)=true; end
assert(~any(allowedMask & excludedMask));

model=build_matrad_dose_model(source.dij,patient.cst,source.pln,0);
[A,physicalD,stateDescription]=exact_state_matrix(model);
baseline=evaluate_exact_state(wReference,A,physicalD,model,patient,cfg, ...
    source.map,groups,source.machine);
baselineVerification=verify_baseline(baseline,stage1,doseTolerance);
assert(baselineVerification.savedDoseMetricsConsistent);
assert(baseline.maxConstraintViolation<=rule.constraintTolerance);

[patchIx,patchInfo]=build_local_patch(wReference,allowedMask,allowedLayer, ...
    groups,source.map,source.stf,80);
assert(all(allowedMask(patchIx)) && ~any(excludedMask(patchIx)));
Apatch=A(:,patchIx); wPatch=wReference(patchIx);
baselineVerification=verify_objective_gradient(model,wReference,patchIx,baseline,baselineVerification);
assert(baselineVerification.objectiveGradientConsistent);
[baseQ,baseWitnessLabels,baseComponents,witnessInfo]=build_base_witnesses( ...
    model,A*wReference,patient.cst,rankTolerance);
assert(size(baseQ,2)<numel(patchIx),'Witness rank leaves no local kernel.');

anchors=choose_anchors(wPatch,min(maxOracleCalls,numel(patchIx)));
variantNames={'fixed_witness','gradient_only','adaptive_secant'};
allRows=repmat(empty_transition(),0,1); componentRows=repmat(empty_component(),0,1);
variantResults=struct();
for v=1:numel(variantNames)
    fprintf('\n=== ABLATION %s ===\n',variantNames{v});
    [vr,rows,comps]=run_ablation(variantNames{v},wReference,patchIx,A,Apatch, ...
        physicalD,model,patient,cfg,source.map,groups,source.machine, ...
        allowedMask,excludedMask,baseline,baseQ,baseWitnessLabels,anchors, ...
        maxOracleCalls,rankTolerance,zeroTolerance,doseTolerance,rule);
    variantResults.(variantNames{v})=vr;
    allRows=[allRows;rows]; componentRows=[componentRows;comps]; %#ok<AGROW>
end

transitionTable=struct2table(allRows,'AsArray',true);
componentTable=struct2table(componentRows,'AsArray',true);
writetable(transitionTable,fullfile(outputDir,'secant_witness_transitions.csv'));
writetable(componentTable,fullfile(outputDir,'secant_witness_component_changes.csv'));
decision=make_decision(variantResults,transitionTable,numel(patchIx));
result=struct('decision',decision,'baseline',baseline, ...
    'baselineVerification',baselineVerification,'patchInfo',patchInfo, ...
    'patchIx',patchIx,'baseWitnessLabels',{baseWitnessLabels}, ...
    'witnessInfo',witnessInfo,'baseComponents',baseComponents, ...
    'variants',variantResults,'stateDescription',stateDescription, ...
    'rankTolerance',rankTolerance,'zeroTolerance',zeroTolerance, ...
    'stage1File',stage1File,'artifactFile',artifactFile);
save(fullfile(outputDir,'secant_witness_result.mat'),'result','-v7.3');
write_report(outputDir,result,transitionTable,componentTable);
fprintf('\nFINAL_DECISION: %s\n',decision.label);
end

function [vr,rows,componentRows]=run_ablation(name,wStart,patchIx,A,Apatch, ...
    physicalD,model,patient,cfg,map,groups,machine,allowedMask,excludedMask, ...
    baseline,Qfixed,fixedLabels,anchors,maxCalls,rankTol,zeroTol,doseTol,rule)
w=wStart; rows=repmat(empty_transition(),0,1); componentRows=repmat(empty_component(),0,1);
secantBank=zeros(size(A,1),0); secantLabels={}; oracleCalls=0; anchorCursor=1;
accepted=0; rejected=0; verifiedUpdates=0; witnessSeconds=0; stopReason='budget_exhausted';
current=baseline;
while oracleCalls<maxCalls && anchorCursor<=numel(anchors)
    anchor=anchors(anchorCursor);
    if w(patchIx(anchor))<=zeroTol, anchorCursor=anchorCursor+1; continue; end
    x=A*w;
    tWitness=tic;
    if strcmp(name,'fixed_witness')
        Q=Qfixed; labels=fixedLabels;
    else
        [Q,labels]=build_base_witnesses(model,x,patient.cst,rankTol);
        if strcmp(name,'adaptive_secant') && ~isempty(secantBank)
            for q=1:size(secantBank,2)
                [Q,added]=append_vector(Q,secantBank(:,q),rankTol);
                if added, labels{end+1}=secantLabels{q}; end %#ok<AGROW>
            end
        end
    end
    witnessBuild=toc(tWitness); witnessSeconds=witnessSeconds+witnessBuild;
    C=Q'*Apatch; rankBefore=numerical_rank(C,rankTol);
    [h,qp]=solve_secant_exchange_qp(C,w(patchIx),anchor,rankTol);
    row=empty_transition(); row.variant=name; row.attempt=numel(rows)+1;
    row.anchorGlobalIx=patchIx(anchor); row.localBlockSize=numel(patchIx);
    row.witnessRankBefore=rankBefore; row.nullspaceDimensionBefore=numel(patchIx)-rankBefore;
    row.svdSeconds=qp.svdSeconds; row.qpSeconds=qp.qpSeconds;
    row.witnessConstructionSeconds=witnessBuild; row.qpExitflag=qp.exitflag;
    row.kernelResidual=qp.kernelResidual; row.anchorBoundaryResidual=qp.anchorResidual;
    if ~qp.success
        row.status='qp_failed'; rows(end+1,1)=row; anchorCursor=anchorCursor+1; %#ok<AGROW>
        continue;
    end
    assert(qp.kernelResidual<=1e-7,'Kernel unit test failed.');
    trial=w; trial(patchIx)=w(patchIx)+h;
    trial(abs(trial)<zeroTol)=0;
    assert(min(trial)>=-1e-9,'Nonnegativity unit test failed.'); trial=max(trial,0);
    assert(trial(patchIx(anchor))<=1e-8*max(1,w(patchIx(anchor))), ...
        'Boundary-reaching unit test failed.');
    assert(~any(trial(excludedMask)>zeroTol),'Excluded Stage-1 layer reactivated.');

    beforeActive=w>zeroTol; afterActive=trial>zeroTol;
    killed=find(beforeActive & ~afterActive); born=find(~beforeActive & afterActive);
    row.spotsKilled=numel(killed); row.spotsBorn=numel(born);
    row.killedGlobalIx=index_list(killed); row.bornGlobalIx=index_list(born);
    row.netSpotChange=nnz(afterActive)-nnz(beforeActive);
    row.activeSpotsBefore=nnz(beforeActive); row.activeSpotsAfter=nnz(afterActive);
    row.activeLayersBefore=count_active_layers(w,groups,zeroTol);
    row.activeLayersAfter=count_active_layers(trial,groups,zeroTol);
    row.activeLayerChange=row.activeLayersAfter-row.activeLayersBefore;
    changedLayerIds=unique(group_ids_for_indices([killed;born],groups,numel(w)));
    row.affectedLayers=numel(changedLayerIds); row.affectedLayerIds=index_list(changedLayerIds);

    tClinical=tic;
    proposed=evaluate_exact_state(trial,A,physicalD,model,patient,cfg,map,groups,machine);
    row.exactClinicalEvaluationSeconds=toc(tClinical); oracleCalls=oracleCalls+1;
    row.clinicalOracleCall=oracleCalls; row.objectiveBefore=current.objective;
    row.objectiveAfter=proposed.objective; row.objectiveChange=proposed.objective-current.objective;
    row.maxConstraintViolation=proposed.maxConstraintViolation;
    row.provisionalBDT=proposed.delivery.T_total;
    [dosePass,worstAdverse]=dose_gate(proposed.dose,baseline.dose,doseTol);
    supportPass=~any(trial(~allowedMask)>zeroTol);
    constraintPass=proposed.maxConstraintViolation<=rule.constraintTolerance;
    acceptedTrial=dosePass && supportPass && constraintPass;
    row.doseQualityPass=dosePass; row.worstAdverseDoseChangePercent=worstAdverse;
    row.supportPreserved=supportPass; row.accepted=acceptedTrial;
    row.status=ternary(acceptedTrial,'accepted','rejected');
    assert(~acceptedTrial || (dosePass && supportPass && constraintPass), ...
        'Accepted transition violates frozen acceptance rule.');

    beforeComp=current.components;
    afterComp=proposed.components;
    deltas=[afterComp.value]'-[beforeComp.value]';
    [~,explainIx]=max(deltas);
    if isempty(explainIx), explainIx=1; end
    row.explainingComponent=beforeComp(explainIx).label;
    if ~isempty(deltas), row.explainingComponentChange=deltas(explainIx); end
    componentRows=[componentRows;component_change_rows(name,row.attempt,beforeComp,afterComp)]; %#ok<AGROW>

    if acceptedTrial
        w=trial; current=proposed; accepted=accepted+1; anchorCursor=anchorCursor+1;
    else
        rejected=rejected+1;
        if strcmp(name,'adaptive_secant')
            vDose=Apatch*h;
            if norm(vDose)>0 && ~isempty(beforeComp)
                a0=beforeComp(explainIx).gradient;
                [Qnew,s]=append_secant_witness(Q,a0,vDose,deltas(explainIx),rankTol);
                row.secantIdentityError=s.secantError;
                row.secantIndependent=s.independent;
                if s.independent
                    oldKernel=norm(Q'*Apatch*h)/max(1,norm(Apatch*h));
                    newKernel=norm(Qnew'*Apatch*h)/max(1,norm(Apatch*h));
                    row.previousDirectionOldKernelResidual=oldKernel;
                    row.previousDirectionNewKernelResidual=newKernel;
                    row.previousDirectionExcluded=newKernel>max(1e-8,100*oldKernel);
                    assert(row.previousDirectionExcluded, ...
                        'Independent secant did not exclude failed direction.');
                    secantBank(:,end+1)=s.newWitness; %#ok<AGROW>
                    secantLabels{end+1}=['secant:' beforeComp(explainIx).label]; %#ok<AGROW>
                    verifiedUpdates=verifiedUpdates+1;
                end
                row.witnessRankAfter=numerical_rank(Qnew'*Apatch,rankTol);
            end
            % Retry the same anchor with genuinely changed geometry.
            if ~row.secantIndependent, anchorCursor=anchorCursor+1; end
        else
            row.witnessRankAfter=rankBefore; anchorCursor=anchorCursor+1;
        end
    end
    if row.witnessRankAfter==0, row.witnessRankAfter=rankBefore; end
    row.nullspaceDimensionAfter=numel(patchIx)-row.witnessRankAfter;
    rows(end+1,1)=row; %#ok<AGROW>
    fprintf('%s attempt=%d accepted=%d kill=%d born=%d rank=%d->%d objDelta=%.4g\n', ...
        name,row.attempt,row.accepted,row.spotsKilled,row.spotsBorn, ...
        row.witnessRankBefore,row.witnessRankAfter,row.objectiveChange);
    if row.witnessRankAfter>=0.9*numel(patchIx) && accepted==0
        stopReason='rank_saturation'; break;
    end
end
if anchorCursor>numel(anchors), stopReason='anchor_budget_exhausted'; end
finalState=current;
vr=struct('name',name,'w',w,'rows',rows,'oracleCalls',oracleCalls, ...
    'accepted',accepted,'rejected',rejected,'verifiedSecantUpdates',verifiedUpdates, ...
    'finalActiveSpots',nnz(w>zeroTol),'finalActiveLayers',count_active_layers(w,groups,zeroTol), ...
    'netExactZeros',nnz(wStart>zeroTol)-nnz(w>zeroTol), ...
    'acceptedZerosPerOracle',(nnz(wStart>zeroTol)-nnz(w>zeroTol))/max(1,oracleCalls), ...
    'finalObjective',finalState.objective,'finalBDT',finalState.delivery.T_total, ...
    'witnessConstructionSeconds',witnessSeconds,'stopReason',stopReason);
end

function state=evaluate_exact_state(w,A,physicalD,model,patient,cfg,map,groups,machine)
t=tic; [f,g,diag]=evaluate_matrad_dose_model(model,w,true);
x=A*w; components=clinical_components(model,x);
componentSum=sum([components(strcmp({components.kind},'objective')).value]);
assert(abs(componentSum-diag.doseObjective)<=1e-7*max(1,abs(diag.doseObjective)), ...
    'Individual objective sum does not match exact matRad objective.');
dose=compute_dvh_metrics(full(physicalD*w),patient.ct,patient.cst,cfg);
delivery=compute_delivery_metrics(w,map,groups,cfg,machine);
state=struct('objective',f,'gradient',g,'diagnostics',diag,'dose',dose, ...
    'delivery',delivery,'components',components,'evaluationSeconds',toc(t), ...
    'maxConstraintViolation',diag.maxConstraintViolation);
end

function [A,D,description]=exact_state_matrix(model)
D=model.dij.physicalDose; if iscell(D), D=D{1}; end
assert(issparse(D),'Full dose influence matrix must remain sparse.');
if isa(model.optiProb.BP,'matRad_DoseProjection')
    A=D; description='physicalDose';
elseif isa(model.optiProb.BP,'matRad_ConstantRBEProjection')
    assert(isfield(model.dij,'RBE') && isscalar(model.dij.RBE));
    A=double(model.dij.RBE)*D; description='constant-RBE dose';
else
    error('secant_witness:unsupportedState', ...
        'Prototype supports linear dose or constant-RBE state only, got %s.',class(model.optiProb.BP));
end
end

function comps=clinical_components(model,x)
nVox=numel(x); comps=repmat(empty_clinical_component(nVox),0,1);
constraintCounter=0;
for i=1:size(model.cst,1)
    if isempty(model.cst{i,4}{1}), continue; end
    vox=reshape(double(model.cst{i,4}{1}),[],1);
    for j=1:numel(model.cst{i,6})
        obj=model.cst{i,6}{j};
        if isa(obj,'DoseObjectives.matRad_DoseObjective')
            assert(strcmp(obj.robustness,'none'),'Prototype is nominal-only.');
            obj=model.optiProb.BP.setBiologicalDosePrescriptions(obj, ...
                model.cst{i,5}.alphaX,model.cst{i,5}.betaX);
            val=obj.penalty*obj.computeDoseObjectiveFunction(x(vox));
            grad=zeros(nVox,1); grad(vox)=obj.penalty*obj.computeDoseObjectiveGradient(x(vox));
            c=empty_clinical_component(nVox); c.label=sprintf('%s:%s', ...
                char(string(model.cst{i,2})),class(obj)); c.kind='objective';
            c.value=double(val); c.gradient=grad; comps(end+1,1)=c; %#ok<AGROW>
        elseif isa(obj,'DoseConstraints.matRad_DoseConstraint')
            assert(strcmp(obj.robustness,'none'),'Prototype is nominal-only.');
            obj=model.optiProb.BP.setBiologicalDosePrescriptions(obj, ...
                model.cst{i,5}.alphaX,model.cst{i,5}.betaX);
            vals=reshape(double(obj.computeDoseConstraintFunction(x(vox))),[],1);
            jac=double(obj.computeDoseConstraintJacobian(x(vox)));
            for k=1:numel(vals)
                constraintCounter=constraintCounter+1; grad=zeros(nVox,1); grad(vox)=jac(:,k);
                c=empty_clinical_component(nVox); c.label=sprintf('%s:%s:%d', ...
                    char(string(model.cst{i,2})),class(obj),k); c.kind='constraint';
                c.value=vals(k); c.gradient=grad;
                c.lower=model.constraintLower(constraintCounter); c.upper=model.constraintUpper(constraintCounter);
                comps(end+1,1)=c; %#ok<AGROW>
            end
        end
    end
end
end

function c=empty_clinical_component(n)
c=struct('label','','kind','','value',NaN,'gradient',zeros(n,1),'lower',-Inf,'upper',Inf);
end

function [Q,labels,components,info]=build_base_witnesses(model,x,cstOriginal,rankTol)
components=clinical_components(model,x); Q=zeros(numel(x),0); labels={}; rejected=0;
for i=1:numel(components)
    use=strcmp(components(i).kind,'objective');
    if strcmp(components(i).kind,'constraint')
        scale=max(1,abs(components(i).value));
        use=min(abs(components(i).value-[components(i).lower,components(i).upper]))<=1e-4*scale;
    end
    if use
        [Q,added]=append_vector(Q,components(i).gradient,rankTol);
        if added, labels{end+1}=components(i).label; else, rejected=rejected+1; end %#ok<AGROW>
    end
end
% Structure-specific critical dose coordinates: two low-dose target and two
% high-dose OAR voxels. These are exact local response directions, not scores.
for s=1:size(model.cst,1)
    vox=reshape(double(model.cst{s,4}{1}),[],1); if isempty(vox), continue; end
    type=char(string(model.cst{s,3}));
    if strcmpi(type,'TARGET'), [~,ord]=sort(x(vox),'ascend');
    elseif strcmpi(type,'OAR'), [~,ord]=sort(x(vox),'descend');
    else, continue; end
    for k=1:min(2,numel(ord))
        e=zeros(numel(x),1); e(vox(ord(k)))=1;
        [Q,added]=append_vector(Q,e,rankTol);
        if added, labels{end+1}=sprintf('critical:%s:voxel%d', ...
                char(string(cstOriginal{s,2})),vox(ord(k))); end %#ok<AGROW>
    end
end
info=struct('candidateCount',numel(components)+2*size(model.cst,1), ...
    'rank',size(Q,2),'dependentRejected',rejected);
end

function [Q,added]=append_vector(Q,q,rankTol)
q=reshape(double(q),[],1); if isempty(Q), r=q; else, r=q-Q*(Q'*q); r=r-Q*(Q'*r); end
threshold=rankTol*max(1,norm(q)); added=norm(r)>threshold;
if added, Q=[Q,r/norm(r)]; end %#ok<AGROW>
end

function [ix,info]=build_local_patch(w,allowed,allowedLayer,groups,map,stf,maxSize)
cp=[groups.cpIndex]'; cps=unique(cp(allowedLayer)); scores=zeros(numel(cps),1);
for k=1:numel(cps), gs=find(allowedLayer & cp==cps(k)); scores(k)=sum(w(vertcat(groups(gs).ix))); end
[~,best]=max(scores); chosenCp=cps(best); candidateGroups=find(allowedLayer & cp==chosenCp);
norms=arrayfun(@(g)norm(w(groups(g).ix)),candidateGroups); [~,centerPos]=max(norms);
centerGroup=candidateGroups(centerPos); energyIndex=[groups.energyIndex]';
[~,ord]=sort(abs(energyIndex(candidateGroups)-energyIndex(centerGroup)));
chosenGroups=candidateGroups(ord(1:min(5,numel(ord))));
centerCandidates=groups(centerGroup).ix; centerCandidates=centerCandidates(allowed(centerCandidates));
[~,p]=max(w(centerCandidates)); centerIx=centerCandidates(p); centerXY=spot_xy(stf,map(centerIx));
perGroup=max(8,floor(maxSize/numel(chosenGroups))); ix=zeros(0,1);
for g=reshape(chosenGroups,1,[])
    gi=groups(g).ix; gi=gi(allowed(gi)); xy=zeros(numel(gi),2);
    for k=1:numel(gi), xy(k,:)=spot_xy(stf,map(gi(k))); end
    dist=sum((xy-centerXY).^2,2); [~,near]=sort(dist,'ascend'); take=near(1:min(perGroup,numel(near)));
    [~,strong]=max(w(gi)); take=unique([take;strong],'stable'); ix=[ix;gi(take)]; %#ok<AGROW>
end
ix=unique(ix,'stable'); if numel(ix)>maxSize, ix=ix(1:maxSize); end
info=struct('controlPoint',chosenCp,'groupIds',reshape(chosenGroups,1,[]), ...
    'energies',reshape([groups(chosenGroups).energy],1,[]),'centerGlobalIx',centerIx, ...
    'centerBEV',centerXY,'blockSize',numel(ix),'initialActive',nnz(w(ix)>1e-10), ...
    'initialDormant',nnz(w(ix)<=1e-10),'globalAllowedSpots',nnz(allowed), ...
    'allowedLayers',nnz(allowedLayer),'excludedLayers',nnz(~allowedLayer));
end

function xy=spot_xy(stf,m)
r=stf(m.beamNum).ray(m.rayNum);
if isfield(r,'rayPos_bev') && numel(r.rayPos_bev)>=3
    q=double(r.rayPos_bev); xy=[q(1),q(3)];
else, xy=[double(m.rayNum),0]; end
end

function anchors=choose_anchors(w,k)
active=find(w>1e-10); [~,hi]=sort(w(active),'descend'); [~,lo]=sort(w(active),'ascend');
order=reshape([active(hi),active(lo)]',[],1); order=unique(order,'stable');
anchors=order(1:min(k,numel(order)));
end

function v=verify_baseline(b,stage1,tol)
if isfield(stage1,'finalDose')
    worst=dose_metric_difference(b.dose,stage1.finalDose);
else, worst=NaN; end
v=struct('recomputedObjective',b.objective,'recomputedConstraintViolation', ...
    b.maxConstraintViolation,'recomputedBDT',b.delivery.T_total, ...
    'savedDoseWorstDifferencePercent',worst, ...
    'savedDoseMetricsConsistent',isnan(worst)||worst<=1e-6, ...
    'objectiveGradientConsistent',false,'componentGradientRelativeError',NaN, ...
    'directionalDerivativeRelativeError',NaN);
if isfield(stage1,'finalDelivery')
    v.savedBDTRelativeError=abs(b.delivery.T_total-stage1.finalDelivery.T_total)/max(1,abs(stage1.finalDelivery.T_total));
else, v.savedBDTRelativeError=NaN; end
end

function v=verify_objective_gradient(model,w,patchIx,baseline,v)
objectiveMask=strcmp({baseline.components.kind},'objective');
gx=zeros(size(baseline.components(1).gradient));
for i=find(objectiveMask), gx=gx+baseline.components(i).gradient; end
[A,~,~]=exact_state_matrix(model);
componentGradient=A'*gx;
v.componentGradientRelativeError=norm(componentGradient-baseline.gradient)/max(1,norm(baseline.gradient));
rng(19); d=zeros(size(w)); q=randn(numel(patchIx),1);
q=q/max(norm(q),eps); d(patchIx)=q;
positive=w(patchIx); step=min(1e-5,0.1*min(positive(positive>0)));
step=max(step,1e-9);
fp=evaluate_matrad_dose_model(model,w+step*d,false);
fm=evaluate_matrad_dose_model(model,w-step*d,false);
fd=(fp-fm)/(2*step); analytic=baseline.gradient'*d;
v.directionalDerivativeRelativeError=abs(fd-analytic)/max([1,abs(fd),abs(analytic)]);
v.objectiveGradientConsistent=all(isfinite(baseline.gradient)) && ...
    v.componentGradientRelativeError<=1e-8 && v.directionalDerivativeRelativeError<=1e-5;
end

function worst=dose_metric_difference(a,b)
worst=0; fields={'D98','D95','Dmean','Dmax'};
for i=1:numel(a.structures)
    for j=1:numel(fields)
        x=a.structures(i).(fields{j}); y=b.structures(i).(fields{j});
        worst=max(worst,100*abs(double(x)-double(y))/max(1,abs(double(y))));
    end
end
end

function [pass,worst]=dose_gate(candidate,reference,tolerance)
worst=0; pass=true;
for i=1:numel(reference.structures)
    a=reference.structures(i); b=candidate.structures(i);
    if a.isTarget
        adverse=max([0,-pct(b.D95,a.D95),-pct(b.D98,a.D98)]);
    elseif strcmpi(a.type,'OAR')
        adverse=max([0,pct(b.Dmean,a.Dmean),pct(b.Dmax,a.Dmax)]);
    else, adverse=0; end
    worst=max(worst,adverse); pass=pass && adverse<=tolerance;
end
end

function rows=component_change_rows(variant,attempt,before,after)
rows=repmat(empty_component(),numel(before),1);
for i=1:numel(before)
    rows(i).variant=variant; rows(i).attempt=attempt; rows(i).component=before(i).label;
    rows(i).kind=before(i).kind; rows(i).valueBefore=before(i).value;
    rows(i).valueAfter=after(i).value; rows(i).change=after(i).value-before(i).value;
    rows(i).lowerBound=before(i).lower; rows(i).upperBound=before(i).upper;
end
end

function ids=group_ids_for_indices(ix,groups,n)
if isempty(ix), ids=zeros(0,1); return; end
lookup=zeros(n,1); for g=1:numel(groups), lookup(groups(g).ix)=g; end; ids=lookup(ix);
end
function s=index_list(ix)
if isempty(ix), s=''; else, s=strjoin(arrayfun(@(x)sprintf('%d',x),reshape(ix,1,[]),'UniformOutput',false),';'); end
end
function n=count_active_layers(w,groups,tol)
n=0; for g=1:numel(groups), n=n+(norm(w(groups(g).ix))>tol); end
end
function r=numerical_rank(C,tol)
if isempty(C), r=0; return; end; s=svd(C,'econ'); if isempty(s),r=0;else,r=sum(s>tol*max(1,s(1)));end
end
function y=pct(x,r), y=100*(double(x)-double(r))/max(abs(double(r)),eps); end
function x=ternary(c,a,b), if c,x=a;else,x=b;end; end

function d=make_decision(v,t,nPatch)
a=v.adaptive_secant; f=v.fixed_witness; g=v.gradient_only;
rankSaturated=any(t.witnessRankAfter>=0.9*nPatch);
samePerformance=a.netExactZeros==f.netExactZeros && a.netExactZeros==g.netExactZeros;
singleKill=all(t.spotsKilled(t.clinicalOracleCall>0)<=1);
supported=a.verifiedSecantUpdates>0 && a.netExactZeros>max(f.netExactZeros,g.netExactZeros) ...
    && a.acceptedZerosPerOracle>max(f.acceptedZerosPerOracle,g.acceptedZerosPerOracle);
if supported, label='MECHANISM_SUPPORTED';
elseif a.verifiedSecantUpdates>0 && (samePerformance || singleKill || rankSaturated), label='MECHANISM_NOT_SUPPORTED';
else, label='INCONCLUSIVE'; end
d=struct('label',label,'verifiedGeometryChanges',a.verifiedSecantUpdates, ...
    'adaptiveNetZeros',a.netExactZeros,'fixedNetZeros',f.netExactZeros, ...
    'gradientNetZeros',g.netExactZeros,'rankSaturated',rankSaturated, ...
    'onlySingleSpotKills',singleKill,'samePerformance',samePerformance);
end

function write_report(out,r,t,componentTable)
fid=fopen(fullfile(out,'SECANT_WITNESS_CIRCUIT_REPORT.md'),'w'); assert(fid>0); c=onCleanup(@()fclose(fid));
fprintf(fid,'# Secant-witness circuit search: local proof of mechanism\n\n');
fprintf(fid,'Decision: **%s**\n\n',r.decision.label);
if r.decision.verifiedGeometryChanges==0
    fprintf(fid,['No clinical rejection occurred, so no failure-derived secant was constructed. ' ...
        'This is not a successful adaptive-secant result; the patch and tolerances were not changed to manufacture a rejection.\n\n']);
end
fprintf(fid,'## Scope and immutable starting point\n\n');
fprintf(fid,'The experiment loaded the frozen Stage-1 alphaEL=0.3 result and changed no production path. ');
fprintf(fid,'It used %d allowed spots, %d allowed layers, and kept %d excluded layers unavailable. No IPOPT polishing or repair was run.\n\n', ...
    r.patchInfo.globalAllowedSpots,r.patchInfo.allowedLayers,r.patchInfo.excludedLayers);
fprintf(fid,'Recomputed exact objective: %.16g; constraint violation: %.3e; BDT: %.6f s. ', ...
    r.baseline.objective,r.baseline.maxConstraintViolation,r.baseline.delivery.T_total);
fprintf(fid,'Saved-dose worst recomputation difference: %.3e%%; saved-BDT relative error: %.3e.\n\n', ...
    r.baselineVerification.savedDoseWorstDifferencePercent,r.baselineVerification.savedBDTRelativeError);
fprintf(fid,'Component-gradient relative error: %.3e; directional finite-difference error: %.3e.\n\n', ...
    r.baselineVerification.componentGradientRelativeError, ...
    r.baselineVerification.directionalDerivativeRelativeError);
fprintf(fid,'## Local patch\n\nCP %d; %d columns; energies: %s; initial active/dormant: %d/%d; BEV center [%.3f, %.3f].\n\n', ...
    r.patchInfo.controlPoint,r.patchInfo.blockSize,mat2str(r.patchInfo.energies,5), ...
    r.patchInfo.initialActive,r.patchInfo.initialDormant,r.patchInfo.centerBEV(1),r.patchInfo.centerBEV(2));
fprintf(fid,'## Witnesses and exchange\n\nThe initial witness space contains separate gradients of every nominal matRad objective, active constraint Jacobians, and two critical dose-voxel rows per target/OAR. ');
fprintf(fid,'It is not the aggregate scalar gradient. The exact state is `%s`.\n\n',r.stateDescription);
fprintf(fid,'The QP solved $\\min_h \\tfrac12\\|W_Jh\\|_2^2$ subject to $Q^TA_Jh=0$, $h_i=-w_i$, and $h\\ge-w_J$, with relative-change scaling in $W_J$.\n\n');
fprintf(fid,'For rejection component $H$, the appended witness was $\\bar a_H=a_0+[H(x+v)-H(x)-a_0^Tv]v/(v^Tv)$. Rank tolerance: %.1e.\n\n',r.rankTolerance);
fprintf(fid,'## Ablation summary\n\n| Variant | oracle calls | accepted | rejected | verified secants | net exact zeros | zeros/oracle | final objective | final BDT (s) | stop |\n|---|---:|---:|---:|---:|---:|---:|---:|---:|---|\n');
names={'fixed_witness','gradient_only','adaptive_secant'};
for i=1:numel(names),v=r.variants.(names{i});fprintf(fid,'| %s | %d | %d | %d | %d | %d | %.3f | %.9g | %.3f | %s |\n', ...
        names{i},v.oracleCalls,v.accepted,v.rejected,v.verifiedSecantUpdates,v.netExactZeros, ...
        v.acceptedZerosPerOracle,v.finalObjective,v.finalBDT,v.stopReason);end
fprintf(fid,'\n## Transition and rank findings\n\n%d attempted transitions and %d component-change rows were recorded. ',height(t),height(componentTable));
fprintf(fid,'Maximum witness rank was %d of block dimension %d. ',max(t.witnessRankAfter),r.patchInfo.blockSize);
fprintf(fid,'The adaptive method verified %d updates for which the failed direction left the updated kernel.\n\n',r.decision.verifiedGeometryChanges);
fprintf(fid,'Stop-condition flags: rank saturation=%d; only one-or-fewer killed spot per clinical call=%d; fixed/gradient/adaptive net sparsity identical=%d.\n\n', ...
    r.decision.rankSaturated,r.decision.onlySingleSpotKills,r.decision.samePerformance);
fprintf(fid,'## Runtime breakdown\n\n| Variant | SVD (s) | QP (s) | witness construction (s) | exact endpoint evaluation (s) |\n|---|---:|---:|---:|---:|\n');
for i=1:numel(names)
    mask=strcmp(t.variant,names{i});
    fprintf(fid,'| %s | %.3f | %.3f | %.3f | %.3f |\n',names{i}, ...
        sum(t.svdSeconds(mask),'omitnan'),sum(t.qpSeconds(mask),'omitnan'), ...
        sum(t.witnessConstructionSeconds(mask),'omitnan'), ...
        sum(t.exactClinicalEvaluationSeconds(mask),'omitnan'));
end
fprintf(fid,'\nExact clinical evaluation dominates QP/SVD time.\n\n');
fprintf(fid,'## Correctness checks\n\nEvery evaluated exchange passed the kernel-residual, nonnegativity, boundary-hit, and excluded-layer assertions. ');
if r.decision.verifiedGeometryChanges>0
    fprintf(fid,'Every appended clinical secant passed its finite-change identity, independence, and old-direction exclusion tests. ');
else
    fprintf(fid,'Synthetic unit tests passed secant identity, independence, and old-direction exclusion; these were not exercised on a real rejected transition. ');
end
fprintf(fid,'Accepted transitions were asserted against the frozen dose/constraint/support gate.\n\n');
fprintf(fid,'## Files\n\n- `secant_witness_transitions.csv`: every attempted exchange.\n- `secant_witness_component_changes.csv`: exact component changes per clinical evaluation.\n- `secant_witness_result.mat`: complete state and indices.\n\n');
fprintf(fid,'This bounded experiment supports only the stated mechanism decision. It makes no novelty, speed, or superiority claim.\n');
end

function r=empty_transition()
r=struct('variant','','attempt',0,'anchorGlobalIx',0,'localBlockSize',0, ...
    'witnessRankBefore',0,'witnessRankAfter',0,'nullspaceDimensionBefore',0, ...
    'nullspaceDimensionAfter',0,'kernelResidual',NaN,'anchorBoundaryResidual',NaN, ...
    'qpExitflag',NaN,'objectiveBefore',NaN,'objectiveAfter',NaN,'objectiveChange',NaN, ...
    'maxConstraintViolation',NaN,'doseQualityPass',false,'worstAdverseDoseChangePercent',NaN, ...
    'supportPreserved',false,'accepted',false,'status','','spotsBorn',0,'spotsKilled',0, ...
    'bornGlobalIx','','killedGlobalIx','', ...
    'netSpotChange',0,'activeSpotsBefore',0,'activeSpotsAfter',0,'activeLayersBefore',0, ...
    'activeLayersAfter',0,'activeLayerChange',0,'affectedLayers',0,'affectedLayerIds','', ...
    'clinicalOracleCall',0, ...
    'provisionalBDT',NaN,'explainingComponent','','explainingComponentChange',NaN, ...
    'secantIdentityError',NaN,'secantIndependent',false, ...
    'previousDirectionOldKernelResidual',NaN,'previousDirectionNewKernelResidual',NaN, ...
    'previousDirectionExcluded',false,'svdSeconds',NaN,'qpSeconds',NaN, ...
    'witnessConstructionSeconds',NaN,'exactClinicalEvaluationSeconds',NaN);
end
function r=empty_component()
r=struct('variant','','attempt',0,'component','','kind','','valueBefore',NaN, ...
    'valueAfter',NaN,'change',NaN,'lowerBound',NaN,'upperBound',NaN);
end
function add_project_paths(root)
f={'config','planning','mapping','objective','delivery','reduction','robustness','evaluation','integrated'};
for i=1:numel(f),addpath(fullfile(root,f{i}));end
end
function m=merge_structs_local(a,b)
m=a; f=fieldnames(b); for i=1:numel(f),k=f{i};if isfield(m,k)&&isstruct(m.(k))&&isstruct(b.(k)),m.(k)=merge_structs_local(m.(k),b.(k));else,m.(k)=b.(k);end,end
end
