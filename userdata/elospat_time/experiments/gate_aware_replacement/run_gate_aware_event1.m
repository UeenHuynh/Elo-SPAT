function output=run_gate_aware_event1(outputDir)
%RUN_GATE_AWARE_EVENT1 Gate-aware historical-dormant repair of frozen event 1.
% The matRad objective and all acceptance gates remain unchanged. Clinical
% gate information is used only to propose a bounded support for one exact
% support-constrained matRad/IPOPT solve.
narginchk(1,1);here=fileparts(mfilename('fullpath'));root=fileparts(fileparts(here));
add_paths(root);create_empty(outputDir);
artifactFile=fullfile(root,'results','task25_27_pipeline_20deg.mat');
stage1File=fullfile(root,'results','full_integrated_layer_20deg_job22180','alpha_0p3','alpha_0p3_final.mat');
eventDir=fullfile(root,'results','dormant_replacement_job22325');
source=load(artifactFile,'cfg','machine','dij','stf','pln','map');
stage1=load(stage1File,'final','finalDose','supportMask');
sv=load(fullfile(eventDir,'event_01_survivor.mat'),'checkpoint');sv=sv.checkpoint;
sc=load(fullfile(eventDir,'event_01_scalar_gradient.mat'),'rescueCheckpoint');sc=sc.rescueCheckpoint;
patient=load(source.cfg.patientFile,'ct','cst');cfg=merge_structs(default_config(),source.cfg);
model=build_matrad_dose_model(source.dij,patient.cst,source.pln,0);assert(~model.hasConstraints);
[A,D]=state_matrices(model);[groups,~]=build_group_index(source.map);zeroTol=1e-10;
wStage1=reshape(double(stage1.final.w),[],1);allowed=reshape(logical(stage1.supportMask),[],1);
wInc=reshape(double(sv.wInc),[],1);supportInc=reshape(logical(sv.supportInc),[],1);
historical=reshape(logical(sv.historicalDormant),[],1);deaths=reshape(double(sv.deathIx),[],1);
wFailed=reshape(double(sv.survivor.w),[],1);survivorSupport=reshape(logical(sv.survivor.support),[],1);
assert(nnz(allowed)==136409&&nnz(supportInc)==136360&&nnz(historical)==49&&numel(deaths)==677);
assert(isequal(historical,allowed&~supportInc));
assert(~any(survivorSupport(deaths))&&~any(survivorSupport(~allowed)));

reference=evaluate_plan(wStage1,model,D,patient,cfg,source,groups);
incumbent=evaluate_plan(wInc,model,D,patient,cfg,source,groups);
survivor=evaluate_plan(wFailed,model,D,patient,cfg,source,groups);
scalar=evaluate_plan(sc.solution.w,model,D,patient,cfg,source,groups);
settings=struct('objectiveToleranceFraction',0.02,'doseTolerancePercent',2, ...
    'oarAbsoluteToleranceGy',0.2,'zeroTolerance',zeroTol, ...
    'hotspotToleranceGy',0.05,'targetBandToleranceGy',0.10, ...
    'maxSurvivorDecrease',640,'maxSurvivorIncrease',384,'maxDormantRescue',256);
gateSurvivor=acceptance_gate(survivor,reference,settings,true);
gateScalar=acceptance_gate(scalar,reference,settings,true);
assert(~gateSurvivor.pass&&gateSurvivor.objectivePass&&~gateSurvivor.dosePass);
assert(~gateScalar.pass&&gateScalar.objectivePass&&~gateScalar.dosePass);
assert(abs(survivor.objective-2873.55556690567)<1e-6&&abs(scalar.objective-2861.04035816224)<1e-6);

% Identify the authoritative violated metric and its near-tied hotspot set.
[violations,worst]=gate_violations(survivor.dose,reference.dose,settings);
assert(strcmp(worst.structure,'Lt femoral head')&&strcmp(worst.metric,'Dmax'));
sIx=find(strcmp(patient.cst(:,2),'Lt femoral head'),1);assert(~isempty(sIx));
vox=reshape(double(patient.cst{sIx,4}{1}),[],1);doseFailed=full(D*wFailed);
gValue=max(doseFailed(vox));gAllowed=worst.allowedValue;
hotLocal=find(doseFailed(vox)>=gValue-settings.hotspotToleranceGy);hotVox=vox(hotLocal);
gateRows=D(hotVox,:);gateRhs=gAllowed-doseFailed(hotVox);
assert(all(gateRhs<0));

% Local target response rows around the exact D95/D98 quantile boundaries.
[targetRows,targetRhs,targetInfo]=target_response_rows(D,doseFailed,patient.cst, ...
    survivor.dose,reference.dose,settings);

% Keep objective, gate and target information as separate arrays.
[~,gObjective]=evaluate_matrad_dose_model(model,wFailed,true);
gateSensitivity=max(gateRows,[],1)';targetSensitivity=sqrt(sum(targetRows.^2,1))';

% Bounded survivor working set. Decrease candidates efficiently remove gate
% dose per unit local target response; increase candidates efficiently supply
% target response per unit gate dose. Objective-gradient candidates are kept
% as a separate tie-breaking block, never blended into gate feasibility.
active=survivorSupport&(wFailed>zeroTol);epsScale=1e-12;
decreaseScore=zeros(size(wFailed));decreaseScore(active)= ...
    gateSensitivity(active).*wFailed(active)./(targetSensitivity(active).*wFailed(active)+epsScale);
increaseScore=zeros(size(wFailed));increaseScore(active)= ...
    targetSensitivity(active)./(gateSensitivity(active)+epsScale);
dec=top_indices(decreaseScore,active,settings.maxSurvivorDecrease);
inc=top_indices(increaseScore,active,settings.maxSurvivorIncrease);
objMask=active&gObjective<0;objScore=zeros(size(wFailed));objScore(objMask)=-gObjective(objMask);
obj=top_indices(objScore,objMask,128);
survivorLocal=unique([dec;inc;obj],'stable');

% The complete legal dormant universe before this event consists of the 49
% historical spots plus the 677 current deaths. Excluded layers never enter.
dormant=allowed&~survivorSupport;dormantIx=find(dormant);assert(numel(dormantIx)==726);
J=unique([survivorLocal;dormantIx],'stable');[isDormant,~]=ismember(J,dormantIx);
[Aineq,bineq,rowScale]=local_constraints(gateRows,gateRhs,targetRows,targetRhs,J);
[lb,ub,varScale]=local_bounds(J,isDormant,wFailed,wStage1);

% Phase I: minimize the common normalized feasibility violation rho. This is
% a pure feasibility problem and contains no matRad-objective coefficient.
nJ=numel(J);Aphase=[Aineq,-rowScale];fphase=[zeros(nJ,1);1];
optsLP=optimoptions('linprog','Algorithm','dual-simplex','Display','none');tLP=tic;
[xPhase,phaseValue,phaseExit,phaseOut,lambda]=linprog(fphase,Aphase,bineq,[],[], ...
    [lb;0],[ub;10],optsLP);phaseSeconds=toc(tLP);assert(phaseExit>0);
rho=xPhase(end);deltaPhase=xPhase(1:nJ);

% Dormant reduced costs are induced by the dual clinical-feasibility prices.
% Negative values mean that activating a zero variable relaxes the priced
% hotspot/target defect. No objective gradient enters this price.
dualIneq=lambda.ineqlin;reducedCost=Aineq'*dualIneq;
dormantLocal=find(isDormant);useful=dormantLocal(reducedCost(dormantLocal)<0);
[~,ord]=sort(reducedCost(useful),'ascend');useful=useful(ord);
useful=useful(1:min(settings.maxDormantRescue,numel(useful)));
rescueProposed=J(useful);

% If even the full bounded local universe cannot satisfy the linearized
% feasibility model, a hard-feasible Phase-II QP is impossible.  Stop here
% rather than spending the single authorized exact IPOPT solve on a support
% that the proposal model itself certifies as infeasible.
phaseFeasTol=1e-7;
if rho>phaseFeasTol
    proposal=struct('J',J,'survivorLocal',survivorLocal,'dormantIx',dormantIx, ...
        'rescueProposedByDual',rescueProposed,'rescueIx',zeros(0,1), ...
        'historicalProposed',zeros(0,1),'sameEventProposed',zeros(0,1), ...
        'rho',rho,'phaseFeasibilityTolerance',phaseFeasTol, ...
        'phaseValue',phaseValue,'phaseExit',phaseExit,'phaseOutput',phaseOut, ...
        'phaseSeconds',phaseSeconds,'deltaPhase',deltaPhase, ...
        'reducedCost',reducedCost,'dualIneq',dualIneq,'gateAllowed',gAllowed, ...
        'hotVox',hotVox,'targetInfo',targetInfo,'settings',settings); %#ok<NASGU>
    save(fullfile(outputDir,'gate_aware_proposal.mat'),'proposal','-v7.3');
    write_hotspots(outputDir,hotVox,doseFailed,gAllowed);
    write_pricing(outputDir,dormantIx,J,reducedCost,gateSensitivity,targetSensitivity, ...
        gObjective,historical,deaths,rescueProposed,zeros(0,1),source.map);
    output=finish_without_exact('DORMANT_REPLACEMENT_NOT_SUPPORTED', ...
        'bounded_gate_aware_linear_feasibility_model_infeasible',reference,incumbent, ...
        survivor,scalar,violations,worst,proposal,eventDir,stage1File,artifactFile);
    save(fullfile(outputDir,'gate_aware_complete.mat'),'output','-v7.3');
    write_report(outputDir,output);fprintf('FINAL_STATUS: %s\n',output.status);return;
end

% Phase II is lexicographic: retain the original hard clinical linearization,
% then minimize the frozen objective's first-order model with a small relative
% trust regularizer. Keep the complete bounded dormant pool in this QP.
% Phase-I duals are degenerate when rho=0 and must not pre-eliminate candidate
% births before the Phase-II constrained reduced-cost test.
J2=J;loc2=(1:numel(J))';
A2=Aineq(:,loc2);lb2=lb(loc2);ub2=ub(loc2);scale2=varScale(loc2);
Hbase=spdiags(1./max(scale2,1e-8).^2,0,numel(J2),numel(J2));
eta=max(1e-8,1e-3*max(abs(gObjective(J2).*scale2)));H=eta*Hbase;
optsQP=optimoptions('quadprog','Algorithm','interior-point-convex','Display','none', ...
    'ConstraintTolerance',1e-8,'OptimalityTolerance',1e-8);tQP=tic;
[delta2,qpValue,qpExit,qpOut,qpLambda]=quadprog(H,gObjective(J2),A2,bineq,[],[],lb2,ub2,[],optsQP); %#ok<ASGLU>
qpSeconds=toc(tQP);assert(qpExit>0);
predictedGate=max(doseFailed(hotVox)+gateRows(:,J2)*delta2);
predictedTarget=targetInfo.currentValue+targetRows(:,J2)*delta2;
assert(predictedGate<=gAllowed+1e-6&&all(predictedTarget>=targetInfo.allowedValue-1e-6));

localDormant=~survivorSupport(J2);positiveDormant=localDormant&delta2>max(1e-8,1e-8*max(abs(delta2)));
rescueIx=J2(positiveDormant);
if numel(rescueIx)>settings.maxDormantRescue
    [~,birthOrder]=sort(delta2(positiveDormant),'descend');
    rescueIx=rescueIx(birthOrder(1:settings.maxDormantRescue));
    % Enforce the bounded rescue batch in the proposal model and re-solve.
    % Nonselected dormant variables are fixed at zero; surviving variables
    % remain free to repair the gate and preserve local target response.
    keepDormant=ismember(J2,rescueIx);lbBound=lb2;ubBound=ub2;
    lbBound(localDormant&~keepDormant)=0;ubBound(localDormant&~keepDormant)=0;
    tQPBound=tic;
    [deltaBound,qpValueBound,qpExitBound,qpOutBound,qpLambdaBound]= ...
        quadprog(H,gObjective(J2),A2,bineq,[],[],lbBound,ubBound,[],optsQP); %#ok<ASGLU>
    qpSeconds=qpSeconds+toc(tQPBound);assert(qpExitBound>0);
    delta2=deltaBound;qpValue=qpValueBound;qpExit=qpExitBound;qpOut=qpOutBound;qpLambda=qpLambdaBound;
    predictedGate=max(doseFailed(hotVox)+gateRows(:,J2)*delta2);
    predictedTarget=targetInfo.currentValue+targetRows(:,J2)*delta2;
    assert(predictedGate<=gAllowed+1e-6&&all(predictedTarget>=targetInfo.allowedValue-1e-6));
    positiveDormant=localDormant&delta2>max(1e-8,1e-8*max(abs(delta2)));
    rescueIx=J2(positiveDormant);
    assert(numel(rescueIx)<=settings.maxDormantRescue);
end
historicalProposed=rescueIx(historical(rescueIx));sameEventProposed=rescueIx(ismember(rescueIx,deaths));

% Save the complete cheap proposal before the one authorized exact solve.
proposal=struct('J',J,'J2',J2,'survivorLocal',survivorLocal,'dormantIx',dormantIx, ...
    'rescueProposedByDual',rescueProposed,'rescueIx',rescueIx, ...
    'historicalProposed',historicalProposed,'sameEventProposed',sameEventProposed, ...
    'delta',delta2,'rho',rho,'phaseFeasibilityTolerance',phaseFeasTol, ...
    'phaseValue',phaseValue,'phaseExit',phaseExit, ...
    'phaseOutput',phaseOut,'phaseSeconds',phaseSeconds,'reducedCost',reducedCost, ...
    'dualIneq',dualIneq,'qpValue',qpValue,'qpExit',qpExit,'qpOutput',qpOut, ...
    'qpLambda',qpLambda,'qpSeconds',qpSeconds,'predictedGate',predictedGate, ...
    'predictedTarget',predictedTarget,'gateAllowed',gAllowed,'hotVox',hotVox, ...
    'targetInfo',targetInfo,'settings',settings); %#ok<NASGU>
save(fullfile(outputDir,'gate_aware_proposal.mat'),'proposal','-v7.3');
write_hotspots(outputDir,hotVox,doseFailed,gAllowed);
write_pricing(outputDir,dormantIx,J,reducedCost,gateSensitivity,targetSensitivity, ...
    gObjective,historical,deaths,rescueProposed,rescueIx,source.map);

if isempty(rescueIx)||isempty(historicalProposed)
    output=finish_without_exact('DORMANT_REPLACEMENT_NOT_SUPPORTED', ...
        'gate_aware_proposal_selected_no_historical_dormant_birth',reference,incumbent, ...
        survivor,scalar,violations,worst,proposal,eventDir,stage1File,artifactFile);
    save(fullfile(outputDir,'gate_aware_complete.mat'),'output','-v7.3');write_report(outputDir,output);return;
end

% ONE expensive exact solve. The gate model only proposes support/warm start;
% the unchanged matRad objective and unchanged exact clinical gates decide.
trialSupport=survivorSupport;trialSupport(rescueIx)=true;
wWarm=wFailed;wWarm(J2)=max(0,wFailed(J2)+delta2);wWarm(~trialSupport)=0;
tExact=tic;exact=run_final_reoptimization(source.dij,patient.cst,source.pln,wWarm,cfg,trialSupport,[]);exactSeconds=toc(tExact);
metrics=compute_solver_metrics(exact.resultGUI,struct(),exact.optimizer);
solverAccepted=metrics.converged||metrics.acceptableButNotConverged||contains(lower(metrics.statusMessage),'acceptable');
plan=evaluate_plan(exact.w,model,D,patient,cfg,source,groups);gate=acceptance_gate(plan,reference,settings,solverAccepted);
positive=exact.w>zeroTol;sameReturn=deaths(positive(deaths));historicalBirth=find(historical&positive);
permanent=setdiff(deaths,sameReturn);net=nnz(positive)-nnz(wInc>zeroTol);
supported=gate.pass&&solverAccepted&&~isempty(historicalBirth)&&net<0;
if supported,status='DORMANT_REPLACEMENT_SUPPORTED';else,status='DORMANT_REPLACEMENT_NOT_SUPPORTED';end

write_transition(outputDir,reference,incumbent,survivor,scalar,plan,gate,metrics, ...
    deaths,permanent,sameReturn,historicalBirth,net,exactSeconds,rho,predictedGate,gAllowed);
write_plan_metrics(outputDir,stage1.finalDose,plans_cell(survivor,scalar,plan));
output=struct('status',status,'reference',reference,'incumbent',incumbent, ...
    'survivor',survivor,'scalarComparator',scalar,'violations',violations,'worstGate',worst, ...
    'proposal',proposal,'exact',exact,'exactPlan',plan,'exactGate',gate,'solverMetrics',metrics, ...
    'sameEventReturns',sameReturn,'historicalDormantBirths',historicalBirth, ...
    'permanentDeaths',permanent,'netSupportChange',net,'exactSeconds',exactSeconds, ...
    'stage1File',stage1File,'eventDir',eventDir,'artifactFile',artifactFile);
save(fullfile(outputDir,'gate_aware_complete.mat'),'output','-v7.3');write_report(outputDir,output);
fprintf('FINAL_STATUS: %s\n',status);
end

function [Aineq,bineq,scale]=local_constraints(G,gRhs,T,tRhs,J)
Aineq=[G(:,J);-T(:,J)];bineq=[gRhs;-tRhs];
scale=max(0.05,abs(bineq));scale=reshape(scale,[],1);
end
function [lb,ub,s]=local_bounds(J,isDormant,w,wRef)
s=max(abs(wRef(J)),1);lb=-w(J);lb(isDormant)=0;
ub=max(abs(w(J)),s);ub(isDormant)=max(2*abs(wRef(J(isDormant))),1);
end
function ix=top_indices(score,mask,k)
q=find(mask&isfinite(score)&score>0);[~,o]=sort(score(q),'descend');ix=q(o(1:min(k,numel(o))));
end

function [rows,rhs,info]=target_response_rows(D,dose,cst,current,reference,s)
rows=sparse(0,size(D,2));rhs=zeros(0,1);labels={};cur=[];allowed=[];
for i=1:numel(current.structures)
    if ~current.structures(i).isTarget,continue;end
    cIx=find(strcmp(cst(:,2),current.structures(i).name),1);vox=reshape(double(cst{cIx,4}{1}),[],1);
    for metric={'D95','D98'}
        key=metric{1};q=current.structures(i).(key);band=vox(abs(dose(vox)-q)<=s.targetBandToleranceGy);
        if isempty(band),[~,o]=sort(abs(dose(vox)-q),'ascend');band=vox(o(1:min(64,numel(o))));end
        rows(end+1,:)=mean(D(band,:),1); %#ok<AGROW>
        a=(1-s.doseTolerancePercent/100)*reference.structures(i).(key);
        rhs(end+1,1)=a-q;labels{end+1,1}=sprintf('%s:%s',current.structures(i).name,key); %#ok<AGROW>
        cur(end+1,1)=q;allowed(end+1,1)=a; %#ok<AGROW>
    end
end
info=struct('labels',{labels},'currentValue',cur,'allowedValue',allowed);
end

function [v,worst]=gate_violations(c,r,s)
v=repmat(struct('structure','','metric','','value',NaN,'referenceValue',NaN, ...
    'allowedValue',NaN,'adversePercent',NaN,'violated',false),0,1);
for i=1:numel(c.structures),a=r.structures(i);b=c.structures(i);
    if a.isTarget,tests={{'D95','lower'},{'D98','lower'}};
    elseif strcmpi(a.type,'OAR')&&~strcmpi(a.name,'BODY')&&~strcmpi(a.name,'EXTERNAL'),tests={{'Dmean','upper'},{'Dmax','upper'}};
    else,continue;end
    for j=1:numel(tests),metric=tests{j}{1};sense=tests{j}{2};x=b.(metric);x0=a.(metric);
        if strcmp(sense,'lower'),limit=(1-s.doseTolerancePercent/100)*x0;adv=max(0,100*(x0-x)/max(abs(x0),eps));bad=x<limit;
        else,tol=max(s.oarAbsoluteToleranceGy,s.doseTolerancePercent/100*max(abs(x0),1));limit=x0+tol;adv=max(0,100*(x-x0)/max(abs(x0),eps));bad=x>limit;end
        v(end+1,1)=struct('structure',a.name,'metric',metric,'value',x,'referenceValue',x0,'allowedValue',limit,'adversePercent',adv,'violated',bad); %#ok<AGROW>
    end
end
[~,q]=max([v.adversePercent]);worst=v(q);
end

function p=evaluate_plan(w,model,D,patient,cfg,source,groups)
[f,~,d]=evaluate_matrad_dose_model(model,w,false);p=struct('objective',f, ...
 'maxConstraintViolation',d.maxConstraintViolation,'dose',compute_dvh_metrics(full(D*w),patient.ct,patient.cst,cfg), ...
 'delivery',compute_delivery_metrics(w,source.map,groups,cfg,source.machine),'activeSpots',nnz(w>1e-10),'activeLayers',count_layers(w,groups,1e-10));
end
function gate=acceptance_gate(p,r,s,solverAccepted)
[dosePass,worst]=dose_quality(p.dose,r.dose,s.doseTolerancePercent,s.oarAbsoluteToleranceGy);objectivePass=p.objective<=(1+s.objectiveToleranceFraction)*r.objective;constraintPass=p.maxConstraintViolation<=1e-6;
gate=struct('pass',solverAccepted&&dosePass&&objectivePass&&constraintPass,'dosePass',dosePass,'objectivePass',objectivePass,'constraintPass',constraintPass,'worstAdverse',worst);
end
function [pass,worst]=dose_quality(c,r,pctTol,absTol)
pass=true;worst=0;for i=1:numel(c.structures),a=r.structures(i);b=c.structures(i);if a.isTarget,adverse=max([0,-pct(b.D95,a.D95),-pct(b.D98,a.D98)]);ok=adverse<=pctTol;
 elseif strcmpi(a.type,'OAR')&&~strcmpi(a.name,'BODY')&&~strcmpi(a.name,'EXTERNAL'),dm=b.Dmean-a.Dmean;dx=b.Dmax-a.Dmax;ok=dm<=max(absTol,pctTol/100*max(abs(a.Dmean),1))&&dx<=max(absTol,pctTol/100*max(abs(a.Dmax),1));adverse=max([0,pct(b.Dmean,a.Dmean),pct(b.Dmax,a.Dmax)]);
 else,adverse=0;ok=true;end;pass=pass&&ok;worst=max(worst,adverse);end
end

function write_hotspots(out,vox,dose,allowed)
T=table(vox,dose(vox),repmat(allowed,numel(vox),1),'VariableNames',{'voxelIx','failedDose','allowedDmax'});writetable(T,fullfile(out,'gate_hotspots.csv'));
end
function write_pricing(out,dormant,J,rc,gSens,tSens,gObj,historical,deaths,dualSet,finalSet,map)
[tf,loc]=ismember(dormant,J);assert(all(tf));n=numel(dormant);source=cell(n,1);for i=1:n,source{i}=ternary(historical(dormant(i)),'HISTORICAL_DORMANT','SAME_EVENT_DEATH');end
T=table([map(dormant).globalIx]',source,gObj(dormant),gSens(dormant),tSens(dormant),rc(loc), ...
    ismember(dormant,dualSet),ismember(dormant,finalSet),ismember(dormant,deaths), ...
    'VariableNames',{'globalIx','sourceClass','objectiveGradient','gateSensitivity', ...
    'targetSensitivity','phaseIFeasibilityReducedCost','selectedByPhaseIDual','selectedFinalQP','sameEventDeath'});
writetable(T,fullfile(out,'gate_aware_candidate_pricing.csv'));
end
function write_transition(out,ref,inc,surv,scalar,plan,gate,metrics,deaths,permanent,returns,births,net,seconds,rho,pred,allowed)
method={'survivor_job22325';'scalar_job22325';'gate_aware'};objective=[surv.objective;scalar.objective;plan.objective];dosePass=[0;0;gate.dosePass];objectivePass=[1;1;gate.objectivePass];planPass=[0;0;gate.pass];sameReturns=[0;56;numel(returns)];histBirths=[0;3;numel(births)];permanentDeaths=[numel(deaths);621;numel(permanent)];netChange=[-677;-618;net];bdt=[surv.delivery.T_total;scalar.delivery.T_total;plan.delivery.T_total];wall=[15486.843622;4901.365124;seconds];
T=table(method,objective,dosePass,objectivePass,planPass,sameReturns,histBirths,permanentDeaths,netChange,bdt,wall,repmat(rho,3,1),repmat(pred,3,1),repmat(allowed,3,1));
T.Properties.VariableNames={'method','objective','dosePass','objectivePass','planAccepted','sameEventReturns','historicalDormantBirths','permanentDeaths','netSupportChange','bdt','solverWallSeconds','phaseIRho','predictedGate','allowedGate'};writetable(T,fullfile(out,'gate_aware_transition.csv'));
S=table({metrics.statusMessage},solver_iterations(metrics),ref.objective,inc.objective,'VariableNames',{'solverStatus','iterations','stage1Objective','incumbentObjective'});writetable(S,fullfile(out,'gate_aware_solver.csv'));
end
function write_plan_metrics(out,ref,plans)
names={'survivor_job22325','scalar_job22325','gate_aware'};rows={};for p=1:3,for i=1:numel(ref.structures),a=ref.structures(i);b=plans{p}.dose.structures(i);rows(end+1,:)={names{p},a.name,a.type,a.D95,b.D95,a.D98,b.D98,a.Dmean,b.Dmean,a.Dmax,b.Dmax};end,end %#ok<AGROW>
T=cell2table(rows,'VariableNames',{'method','structure','type','stage1D95','trialD95','stage1D98','trialD98','stage1Dmean','trialDmean','stage1Dmax','trialDmax'});writetable(T,fullfile(out,'gate_aware_clinical_metrics.csv'));
end
function x=plans_cell(a,b,c),x={a,b,c};end
function o=finish_without_exact(status,reason,ref,inc,surv,scalar,v,worst,proposal,eventDir,stage1File,artifactFile)
o=struct('status',status,'reason',reason,'reference',ref,'incumbent',inc,'survivor',surv,'scalarComparator',scalar,'violations',v,'worstGate',worst,'proposal',proposal,'eventDir',eventDir,'stage1File',stage1File,'artifactFile',artifactFile);
end
function write_report(out,o)
fid=fopen(fullfile(out,'GATE_AWARE_REPLACEMENT_REPORT.md'),'w');assert(fid>0);c=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'# Gate-aware dormant replacement, frozen Event 1\n\nStatus: **%s**\n\n',o.status);fprintf(fid,'The matRad objective and clinical gates were not modified. Gate information was used only for support proposal.\n');if isfield(o,'reason'),fprintf(fid,'\nStop reason: `%s`.\n',o.reason);else,fprintf(fid,'\nHistorical births: %d; same-event returns: %d; net support change: %d; exact gate pass: %d.\n',numel(o.historicalDormantBirths),numel(o.sameEventReturns),o.netSupportChange,o.exactGate.pass);end
end
function [A,D]=state_matrices(model),D=model.dij.physicalDose;if iscell(D),D=D{1};end;if isa(model.optiProb.BP,'matRad_DoseProjection'),A=D;elseif isa(model.optiProb.BP,'matRad_ConstantRBEProjection'),A=double(model.dij.RBE)*D;else,error('Linear state required.');end,end
function n=count_layers(w,g,t),n=0;for i=1:numel(g),n=n+any(w(g(i).ix)>t);end,end
function x=pct(a,b),x=100*(double(a)-double(b))/max(abs(double(b)),eps);end
function x=ternary(c,a,b),if c,x=a;else,x=b;end,end
function x=solver_iterations(m),x=NaN;if isfield(m,'iterations'),x=m.iterations;elseif isfield(m,'iterationCount'),x=m.iterationCount;end,end
function add_paths(root),f={'config','planning','mapping','objective','delivery','reduction','evaluation','integrated','robustness'};for i=1:numel(f),addpath(fullfile(root,f{i}));end,end
function create_empty(d),if exist(d,'dir')==7,e=dir(d);e=e(~ismember({e.name},{'.','..'}));assert(isempty(e));else,mkdir(d);end,end
function m=merge_structs(a,b),m=a;f=fieldnames(b);for i=1:numel(f),k=f{i};if isfield(m,k)&&isstruct(m.(k))&&isstruct(b.(k)),m.(k)=merge_structs(m.(k),b.(k));else,m.(k)=b.(k);end,end,end
