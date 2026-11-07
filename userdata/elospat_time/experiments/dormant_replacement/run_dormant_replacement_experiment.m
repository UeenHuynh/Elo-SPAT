function output=run_dormant_replacement_experiment(outputDir)
%RUN_DORMANT_REPLACEMENT_EXPERIMENT Full-universe Stage-2 mechanism test.
% Isolated research path: coherent deaths, survivor control, then matched
% scalar/component-aware dormant add-back.  It does not modify Stage 1.
narginchk(1,1);here=fileparts(mfilename('fullpath'));root=fileparts(fileparts(here));
add_paths(root);create_empty(outputDir);
artifactFile=fullfile(root,'results','task25_27_pipeline_20deg.mat');
stage1File=fullfile(root,'results','full_integrated_layer_20deg_job22180','alpha_0p3','alpha_0p3_final.mat');
continuationFile=fullfile(root,'results','bidirectional_support_job22312','bidirectional_support_complete.mat');
source=load(artifactFile,'cfg','machine','dij','stf','pln','map');
stage1=load(stage1File,'final','finalDose','finalDelivery','supportMask');
prior=load(continuationFile,'output');
patient=load(source.cfg.patientFile,'ct','cst');cfg=merge_structs(default_config(),source.cfg);
model=build_matrad_dose_model(source.dij,patient.cst,source.pln,0);assert(~model.hasConstraints);
[A,D]=state_matrices(model);[groups,~]=build_group_index(source.map);n=numel(stage1.final.w);
wReference=reshape(double(stage1.final.w),[],1);allowed=reshape(logical(stage1.supportMask),[],1);
allowedLayer=false(numel(groups),1);for g=1:numel(groups),allowedLayer(g)=any(allowed(groups(g).ix));end
assert(nnz(allowed)==136409&&nnz(allowedLayer)==748&&nnz(~allowedLayer)==20);
excluded=~allowed;zeroTol=1e-10;reference=evaluate_plan(wReference,model,D,patient,cfg,source,groups);
assert(abs(reference.objective-3139.658625654965)<=1e-8*reference.objective);

settings=struct('zeroTolerance',zeroTol,'objectiveToleranceFraction',0.02, ...
    'doseTolerancePercent',2,'oarAbsoluteToleranceGy',0.2, ...
    'eventModes',{{'layer_block_4','layer_block_12','full_cp_1','full_cp_2','full_cp_3','energy_band_4cp'}}, ...
    'rescueFraction',0.20,'minimumRescue',64,'maximumRescue',512, ...
    'minimumTrueReplacementEvents',1,'maximumEvents',6);

wInc=reshape(double(prior.output.final.w),[],1);supportInc=allowed&(wInc>zeroTol);
assert(strcmp(prior.output.status,'BIDIRECTIONAL_SUPPORT_MECHANISM_SUPPORTED'));
assert(nnz(supportInc)==136360&&nnz(allowed&~supportInc)==49&&~any(supportInc(excluded)));
continuationPlan=evaluate_plan(wInc,model,D,patient,cfg,source,groups);
continuationGate=acceptance_gate(continuationPlan,reference,settings,true);assert(continuationGate.pass);
historicalDormant=allowed&~supportInc;
rows=repmat(empty_transition(),0,1);components=repmat(empty_component(),0,1);
mapping=repmat(empty_mapping(),0,1);trajectory=repmat(empty_trajectory(),0,1);
comparisons=repmat(empty_comparison(),0,1);
trajectory(end+1,1)=trajectory_row(0,'accepted_job22312_continuation',supportInc,wInc,groups,source,reference,allowed,cfg); %#ok<AGROW>
totalTimer=tic;trueReplacementEvents=0;clinicalFailures=0;numericalFailures=0;
for event=1:settings.maximumEvents
    incumbent=evaluate_plan(wInc,model,D,patient,cfg,source,groups);
    [gInc,compInc]=objective_information(model,A,wInc);
    [deathIx,deathInfo]=coherent_death_set(settings.eventModes{event},wInc,supportInc, ...
        historicalDormant,compInc,source.map,source.stf,groups,event);
    if isempty(deathIx),break;end
    survivorSupport=supportInc;survivorSupport(deathIx)=false;
    wDeleted=wInc;wDeleted(deathIx)=0;
    raw=evaluate_plan(wDeleted,model,D,patient,cfg,source,groups);

    survivor=solve_support(source,patient,cfg,wDeleted,survivorSupport,model,D,groups,reference,settings);
    if survivor.solverAccepted&&survivor.gate.pass
        survivorClass='SURVIVOR_REPAIR_SUCCESS';failureType='none';
    elseif survivor.solverAccepted
        survivorClass='SURVIVOR_REPAIR_FAILURE';failureType='clinical_support_failure';clinicalFailures=clinicalFailures+1;
    else
        survivorClass='SURVIVOR_REPAIR_FAILURE';failureType='numerical_solver_failure';numericalFailures=numericalFailures+1;
    end
    rows(end+1,1)=transition_row(event,settings.eventModes{event},'survivor', ...
        incumbent,raw,survivor,numel(deathIx),numel(deathIx),0,0,0, ...
        survivor.plan.activeSpots-incumbent.activeSpots,nnz(historicalDormant), ...
        nnz(allowed&~(survivor.w>zeroTol)),survivorClass,failureType,deathInfo,toc(totalTimer)); %#ok<AGROW>
    components=[components;component_rows(event,'raw_deletion','none',compInc,raw.components); ...
        component_rows(event,'survivor_repair','survivor',compInc,survivor.plan.components)]; %#ok<AGROW>

    checkpoint=struct('event',event,'mode',settings.eventModes{event},'wInc',wInc, ...
        'supportInc',supportInc,'historicalDormant',historicalDormant,'deathIx',deathIx, ...
        'deathInfo',deathInfo,'raw',raw,'survivor',survivor); %#ok<NASGU>
    save(fullfile(outputDir,sprintf('event_%02d_survivor.mat',event)),'checkpoint','-v7.3');
    flush_tables(outputDir,rows,components,mapping,trajectory);
    fprintf('DORMANT event=%d mode=%s deaths=%d survivor=%s obj=%.6f worst=%.4f\n', ...
        event,settings.eventModes{event},numel(deathIx),survivorClass, ...
        survivor.plan.objective,survivor.gate.worstAdverse);

    % Survivor success is a control only: do not turn this branch into another
    % survivor-only compression loop. Numerical failure is not evidence that
    % dormant columns are needed and therefore cannot trigger a valid rescue.
    if ~strcmp(failureType,'clinical_support_failure')
        trajectory(end+1,1)=trajectory_row(event,'no_incumbent_change',supportInc,wInc,groups,source,reference,allowed,cfg); %#ok<AGROW>
        continue;
    end

    failed=survivor.plan;wFailed=survivor.w;[gFailed,compFailed]=objective_information(model,A,wFailed);
    components=[components;component_rows(event,'failed_defect','survivor',compInc,compFailed)]; %#ok<AGROW>
    dormant=allowed&~survivorSupport; % Current deaths plus prior accepted dormant support.
    rescueN=min(settings.maximumRescue,max(settings.minimumRescue,ceil(settings.rescueFraction*numel(deathIx))));
    [scalarIx,scalarScore]=rescue_scalar(gFailed,A,dormant,rescueN);
    [componentIx,componentScore,damageWeights]=rescue_component(compInc,compFailed,A,dormant,rescueN);
    variants={'scalar_gradient','component_aware'};ixSets={scalarIx,componentIx};scoreSets={scalarScore,componentScore};
    candidates=repmat(struct('name','','solve',[],'ix',[],'score',[],'returns',[], ...
        'trueBirths',[],'accepted',false,'net',0),2,1);
    deletedDose=A(:,deathIx)*wInc(deathIx);
    for v=1:2
        rix=ixSets{v};trialSupport=survivorSupport;trialSupport(rix)=true;
        wStart=wFailed;seed=max(wReference(rix),1e-8);wStart(rix)=seed;
        sol=solve_support(source,patient,cfg,wStart,trialSupport,model,D,groups,reference,settings);
        positive=sol.w>zeroTol;currentReturn=rix(ismember(rix,deathIx)&positive(rix));
        trueBirth=rix(historicalDormant(rix)&positive(rix));permanentDeath=setdiff(deathIx,currentReturn);
        net=sol.plan.activeSpots-incumbent.activeSpots;accepted=sol.solverAccepted&&sol.gate.pass&&net<0;
        candidates(v)=struct('name',variants{v},'solve',sol,'ix',rix,'score',scoreSets{v}, ...
            'returns',currentReturn,'trueBirths',trueBirth,'accepted',accepted,'net',net);
        rows(end+1,1)=transition_row(event,settings.eventModes{event},variants{v}, ...
            incumbent,raw,sol,numel(deathIx),numel(permanentDeath),numel(rix), ...
            numel(currentReturn),numel(trueBirth),net,nnz(historicalDormant), ...
            nnz(allowed&~positive),survivorClass,ternary(accepted,'none','rescue_failed'), ...
            deathInfo,toc(totalTimer)); %#ok<AGROW>
        components=[components;component_rows(event,'joint_repair',variants{v},compInc,sol.plan.components)]; %#ok<AGROW>
        mapping=[mapping;mapping_rows(event,variants{v},deathIx,rix,scoreSets{v}, ...
            historicalDormant,sol.w,deletedDose,A,source.map,source.stf,groups)]; %#ok<AGROW>
        rescueCheckpoint=struct('event',event,'variant',variants{v},'deathIx',deathIx, ...
            'rescueIx',rix,'scores',scoreSets{v},'damageWeights',damageWeights, ...
            'historicalDormantBefore',historicalDormant,'solution',sol); %#ok<NASGU>
        save(fullfile(outputDir,sprintf('event_%02d_%s.mat',event,variants{v})), ...
            'rescueCheckpoint','-v7.3');
    end
    comparisons(end+1,1)=comparison_row(event,numel(deathIx),scalarIx,componentIx,candidates); %#ok<AGROW>
    writetable(struct2table(comparisons,'AsArray',true),fullfile(outputDir,'rescue_method_comparison.csv'));

    chosen=choose_candidate(candidates);
    if chosen>0
        oldSupport=supportInc;wInc=candidates(chosen).solve.w;supportInc=allowed&(wInc>zeroTol);
        historicalDormant=allowed&~supportInc;
        if ~isempty(candidates(chosen).trueBirths),trueReplacementEvents=trueReplacementEvents+1;end
        action=['accepted_' candidates(chosen).name];
        assert(nnz(supportInc)<nnz(oldSupport)&&~any(supportInc(excluded)));
    else
        action='no_incumbent_change';
    end
    trajectory(end+1,1)=trajectory_row(event,action,supportInc,wInc,groups,source,reference,allowed,cfg); %#ok<AGROW>
    flush_tables(outputDir,rows,components,mapping,trajectory);
end

% One common final polish of the selected incumbent support.
final=solve_support(source,patient,cfg,wInc,supportInc,model,D,groups,reference,settings);
finalSupportPreserved=~any(abs(final.w(~supportInc))>1e-8*max(1,max(abs(final.w))));
assert(finalSupportPreserved&&~any(supportInc(excluded)));
acceptedRescues=sum([trajectory.event]>0 & contains({trajectory.action},'accepted_'));
if trueReplacementEvents>=settings.minimumTrueReplacementEvents&&final.gate.pass&&final.solverAccepted
    status='DORMANT_REPLACEMENT_SUPPORTED';
elseif clinicalFailures>=2&&trueReplacementEvents==0
    status='DORMANT_REPLACEMENT_NOT_SUPPORTED';
elseif clinicalFailures>0&&acceptedRescues==0
    status='DORMANT_REPLACEMENT_NOT_SUPPORTED';
else
    status='INCONCLUSIVE';
end
output=struct('status',status,'settings',settings,'reference',reference,'rows',rows, ...
    'components',components,'mapping',mapping,'trajectory',trajectory,'comparisons',comparisons,'final',final, ...
    'supportFinal',supportInc,'historicalDormantFinal',historicalDormant, ...
    'trueReplacementEvents',trueReplacementEvents,'clinicalFailures',clinicalFailures, ...
    'numericalFailures',numericalFailures,'wallSeconds',toc(totalTimer), ...
    'stage1File',stage1File,'continuationFile',continuationFile,'continuationPlan',continuationPlan, ...
    'continuationGate',continuationGate,'artifactFile',artifactFile);
save(fullfile(outputDir,'dormant_replacement_complete.mat'),'output','-v7.3');
flush_tables(outputDir,rows,components,mapping,trajectory);write_report(outputDir,output);
fprintf('FINAL_STATUS: %s\n',status);
end

function sol=solve_support(source,patient,cfg,wStart,support,model,D,groups,reference,settings)
t=tic;wStart(~support)=0;result=run_final_reoptimization(source.dij,patient.cst,source.pln,wStart,cfg,support,[]);
wall=toc(t);metrics=compute_solver_metrics(result.resultGUI,struct(),result.optimizer);
solverAccepted=metrics.converged||metrics.acceptableButNotConverged||contains(lower(metrics.statusMessage),'acceptable');
plan=evaluate_plan(result.w,model,D,patient,cfg,source,groups);gate=acceptance_gate(plan,reference,settings,solverAccepted);
sol=struct('w',result.w,'result',result,'metrics',metrics,'solverAccepted',solverAccepted, ...
    'status',metrics.statusMessage,'wallSeconds',wall,'plan',plan,'gate',gate,'support',support);
end

function [death,info]=coherent_death_set(mode,w,support,historical,comps,map,stf,groups,event)
n=numel(w);risk=zeros(n,1);for c=1:numel(comps),risk=risk+(w.*comps(c).weightGradient/max(1,abs(comps(c).value))).^2;end
risk=sqrt(risk);groupRisk=inf(numel(groups),1);groupActive=zeros(numel(groups),1);cp=[groups.cpIndex]';ei=[groups.energyIndex]';
for g=1:numel(groups),ix=groups(g).ix;ix=ix(support(ix)&~historical(ix));groupActive(g)=numel(ix);if ~isempty(ix),groupRisk(g)=sum(risk(ix))/sqrt(numel(ix));end,end
if startsWith(mode,'layer_block_')
    L=str2double(extractAfter(mode,'layer_block_'));[block,score]=best_layer_block(cp,ei,groupRisk,groupActive,L,event);
elseif startsWith(mode,'full_cp_')
    C=str2double(extractAfter(mode,'full_cp_'));[block,score]=best_cp_block(cp,groupRisk,groupActive,C,event);
else
    [block,score]=best_energy_band(cp,ei,groupRisk,groupActive,4,event);
end
death=vertcat(groups(block).ix);death=death(support(death));death=unique(death,'stable');
assert(~isempty(death)&&~any(historical(death)));
pos=zeros(numel(death),2);for i=1:numel(death),pos(i,:)=spot_xy(stf,map(death(i)));end
info=struct('blockScore',score,'groupIds',index_list(block),'controlPoints',index_list(unique([map(death).cpIndex])), ...
    'energies',index_list(unique([map(death).energy])),'meanX',mean(pos(:,1)),'meanY',mean(pos(:,2)));
end

function [block,best]=best_layer_block(cp,ei,risk,count,L,event)
cands={};scores=[];for c=unique(cp)'
    g=find(cp==c&count>0);[~,o]=sort(ei(g));g=g(o);if numel(g)<L,continue;end
    for j=1:numel(g)-L+1,b=g(j:j+L-1);cands{end+1}=b;scores(end+1)=sum(risk(b))*sqrt(sum(count(b)));end %#ok<AGROW>
end
[block,best]=pick_structural(cands,scores,event);
end
function [block,best]=best_cp_block(cp,risk,count,C,event)
u=unique(cp);s=zeros(numel(u),1);for i=1:numel(u),g=find(cp==u(i)&count>0);s(i)=sum(risk(g))*sqrt(sum(count(g)));end
[~,o]=sort(s,'descend');chosen=u(o(1:min(C,numel(o))));block=find(ismember(cp,chosen)&count>0);best=sum(s(o(1:min(C,numel(o)))));
if event<3,[~,q]=sort(s,'ascend');chosen=u(q(1:min(C,numel(q))));block=find(ismember(cp,chosen)&count>0);best=sum(s(q(1:min(C,numel(q)))));end
end
function [block,best]=best_energy_band(cp,ei,risk,count,C,event)
u=unique(cp);cpScore=zeros(numel(u),1);for i=1:numel(u),g=find(cp==u(i)&count>0);cpScore(i)=sum(risk(g));end
[~,o]=sort(cpScore,'descend');chosen=u(o(1:min(C,numel(o))));block=zeros(0,1);
for i=1:numel(chosen),g=find(cp==chosen(i)&count>0);[~,q]=sort(risk(g),'descend');block=[block;g(q(1:min(6,numel(q))))];end %#ok<AGROW>
best=sum(risk(block))*sqrt(sum(count(block)));block=unique(block);
end
function [block,best]=pick_structural(cands,scores,event)
assert(~isempty(cands));[scores,o]=sort(scores,'ascend');q=min(numel(o),max(1,round((0.35+0.10*event)*numel(o))));block=cands{o(q)};best=scores(q);
end

function [ix,score]=rescue_scalar(g,A,dormant,k)
norm2=max(reshape(full(sum(A.^2,1)),[],1),eps);score=zeros(size(g));neg=dormant&g<0;score(neg)=g(neg).^2./(2*norm2(neg));[~,o]=sort(score,'descend');o=o(score(o)>0);ix=o(1:min(k,numel(o)));score=score(ix);
end
function [ix,score,damage]=rescue_component(before,failed,A,dormant,k)
m=numel(before);damage=zeros(m,1);merit=zeros(numel(dormant),1);
for c=1:m
    j=find(strcmp({failed.label},before(c).label),1);assert(~isempty(j));
    damage(c)=max(0,(failed(j).value-before(c).value)/max(1,abs(before(c).value)));
    merit=merit+damage(c)*failed(j).weightGradient/max(1,abs(before(c).value));
end
norm2=max(reshape(full(sum(A.^2,1)),[],1),eps);fullScore=zeros(size(merit));neg=dormant&merit<0;fullScore(neg)=merit(neg).^2./(2*norm2(neg));
[~,o]=sort(fullScore,'descend');o=o(fullScore(o)>0);ix=o(1:min(k,numel(o)));score=fullScore(ix);
end

function chosen=choose_candidate(c)
valid=find([c.accepted]);if isempty(valid),chosen=0;return;end
net=[c(valid).net];best=min(net);valid=valid(net==best);if numel(valid)>1,[~,q]=min(arrayfun(@(x)x.solve.plan.objective,c(valid)));chosen=valid(q);else,chosen=valid;end
end

function r=comparison_row(event,deaths,a,b,c)
r=empty_comparison();r.event=event;r.proposedDeaths=deaths;r.scalarCandidates=numel(a);r.componentCandidates=numel(b);r.overlapCount=numel(intersect(a,b));r.jaccard=r.overlapCount/max(1,numel(union(a,b)));r.scalarAccepted=c(1).accepted;r.componentAccepted=c(2).accepted;r.scalarSameEventReturns=numel(c(1).returns);r.componentSameEventReturns=numel(c(2).returns);r.scalarHistoricalBirths=numel(c(1).trueBirths);r.componentHistoricalBirths=numel(c(2).trueBirths);r.scalarObjective=c(1).solve.plan.objective;r.componentObjective=c(2).solve.plan.objective;r.scalarNetChange=c(1).net;r.componentNetChange=c(2).net;
end

function rows=mapping_rows(event,method,deaths,rescue,scores,historical,wAfter,deletedDose,A,map,stf,groups)
groupOf=zeros(numel(map),1);for g=1:numel(groups),groupOf(groups(g).ix)=g;end
deathCp=unique([map(deaths).cpIndex]);deathLayer=unique(groupOf(deaths));dn=norm(deletedDose);rows=repmat(empty_mapping(),numel(rescue),1);
for q=1:numel(rescue),i=rescue(q);col=A(:,i);rows(q).event=event;rows(q).method=method;rows(q).deathGlobalIx=index_list([map(deaths).globalIx]);
    rows(q).birthGlobalIx=map(i).globalIx;rows(q).sourceClass=ternary(historical(i),'historical_dormant','current_death_return');
    rows(q).score=scores(q);rows(q).finalWeight=wAfter(i);rows(q).becamePositive=wAfter(i)>1e-10;rows(q).cpIndex=map(i).cpIndex;rows(q).energy=map(i).energy;
    rows(q).sameDeathCP=ismember(map(i).cpIndex,deathCp);rows(q).sameDeathLayer=ismember(groupOf(i),deathLayer);
    rows(q).doseCosine=full(col'*deletedDose)/max(eps,norm(col)*dn);[rows(q).spatialDistance,rows(q).energyDistance]=nearest_death(i,deaths,map,stf);end
end
function [sd,ed]=nearest_death(i,deaths,map,stf)
same=find([map(deaths).cpIndex]==map(i).cpIndex);if isempty(same),sd=Inf;ed=Inf;return;end
p=spot_xy(stf,map(i));sd=Inf;ed=Inf;for q=reshape(same,1,[]),j=deaths(q);d=norm(p-spot_xy(stf,map(j)));if d<sd,sd=d;end;ed=min(ed,abs(map(i).energy-map(j).energy));end
end
function p=spot_xy(stf,m),r=stf(m.beamNum).ray(m.rayNum);x=double(r.rayPos_bev);p=[x(1),x(3)];end

function [g,comps]=objective_information(model,A,w)
[~,g]=evaluate_matrad_dose_model(model,w,true);x=A*w;nVox=numel(x);comps=repmat(struct('label','','value',NaN,'weightGradient',zeros(numel(w),1)),0,1);
for i=1:size(model.cst,1),vox=reshape(double(model.cst{i,4}{1}),[],1);if isempty(vox),continue;end
 for j=1:numel(model.cst{i,6}),o=model.cst{i,6}{j};if ~isa(o,'DoseObjectives.matRad_DoseObjective'),continue;end
  o=model.optiProb.BP.setBiologicalDosePrescriptions(o,model.cst{i,5}.alphaX,model.cst{i,5}.betaX);dg=zeros(nVox,1);dg(vox)=o.penalty*o.computeDoseObjectiveGradient(x(vox));
  comps(end+1,1)=struct('label',sprintf('%s:%s',char(string(model.cst{i,2})),class(o)),'value',double(o.penalty*o.computeDoseObjectiveFunction(x(vox))),'weightGradient',A'*dg); %#ok<AGROW>
 end
end
assert(norm(sum(reshape([comps.weightGradient],numel(w),[]),2)-g)<=1e-8*max(1,norm(g)));
end
function vals=component_values(model,x)
vals=repmat(struct('label','','value',NaN),0,1);for i=1:size(model.cst,1),vox=reshape(double(model.cst{i,4}{1}),[],1);if isempty(vox),continue;end
 for j=1:numel(model.cst{i,6}),o=model.cst{i,6}{j};if ~isa(o,'DoseObjectives.matRad_DoseObjective'),continue;end;o=model.optiProb.BP.setBiologicalDosePrescriptions(o,model.cst{i,5}.alphaX,model.cst{i,5}.betaX);vals(end+1,1)=struct('label',sprintf('%s:%s',char(string(model.cst{i,2})),class(o)),'value',double(o.penalty*o.computeDoseObjectiveFunction(x(vox))));end
end,end

function p=evaluate_plan(w,model,D,patient,cfg,source,groups)
[f,~,d]=evaluate_matrad_dose_model(model,w,false);p=struct('objective',f,'components',component_values(model,state_matrix(model)*w), ...
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

function r=transition_row(event,mode,method,inc,raw,sol,deaths,permanent,rescue,returns,births,net,dormantBefore,dormantAfter,survivorClass,failure,info,wall)
r=empty_transition();r.event=event;r.mode=mode;r.method=method;r.incumbentSpots=inc.activeSpots;r.incumbentLayers=inc.activeLayers;r.proposedDeaths=deaths;r.permanentDeaths=permanent;r.rescueCandidates=rescue;r.sameEventReturns=returns;r.historicalDormantBirths=births;r.dormantPoolBefore=dormantBefore;r.dormantPoolAfter=dormantAfter;r.netSupportChange=net;r.objectiveBefore=inc.objective;r.objectiveRawDeletion=raw.objective;r.objectiveAfterRepair=sol.plan.objective;r.worstAdversePercent=sol.gate.worstAdverse;r.objectivePass=sol.gate.objectivePass;r.dosePass=sol.gate.dosePass;r.constraintPass=sol.gate.constraintPass;r.solverAccepted=sol.solverAccepted;r.planAccepted=sol.gate.pass;r.survivorClassification=survivorClass;r.failureClass=failure;r.solverStatus=sol.status;r.solverIterations=solver_iterations(sol.metrics);r.solverWallSeconds=sol.wallSeconds;r.activeSpotsAfter=sol.plan.activeSpots;r.activeLayersAfter=sol.plan.activeLayers;r.bdtAfter=sol.plan.delivery.T_total;r.deathGroupIds=info.groupIds;r.deathControlPoints=info.controlPoints;r.deathEnergies=info.energies;r.wallSeconds=wall;
end
function rows=component_rows(event,phase,method,before,after)
rows=repmat(empty_component(),numel(before),1);for i=1:numel(before),j=find(strcmp({after.label},before(i).label),1);rows(i)=struct('event',event,'phase',phase,'method',method,'component',before(i).label,'incumbentValue',before(i).value,'trialValue',after(j).value,'change',after(j).value-before(i).value);end
end
function r=trajectory_row(event,action,support,w,groups,source,reference,allowed,cfg)
r=empty_trajectory();r.event=event;r.action=action;r.supportSpots=nnz(support);r.activeSpots=nnz(w>1e-10);r.activeLayers=count_layers(w,groups,1e-10);r.dormantSpots=nnz(allowed&~support);r.bdt=compute_delivery_metrics(w,source.map,groups,cfg,source.machine).T_total;r.objectiveReference=reference.objective;
end
function flush_tables(out,rows,comps,map,traj)
if ~isempty(rows),writetable(struct2table(rows,'AsArray',true),fullfile(out,'dormant_transitions.csv'));end
if ~isempty(comps),writetable(struct2table(comps,'AsArray',true),fullfile(out,'clinical_defect_components.csv'));end
if ~isempty(map),writetable(struct2table(map,'AsArray',true),fullfile(out,'death_birth_mapping.csv'));end
if ~isempty(traj),writetable(struct2table(traj,'AsArray',true),fullfile(out,'support_trajectory.csv'));end
end
function write_report(out,o)
fid=fopen(fullfile(out,'DORMANT_REPLACEMENT_REPORT.md'),'w');assert(fid>0);c=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'# Dormant replacement mechanism experiment\n\nStatus: **%s**\n\n',o.status);fprintf(fid,'True replacement events: %d; clinical survivor failures: %d; numerical survivor failures: %d.\n\n',o.trueReplacementEvents,o.clinicalFailures,o.numericalFailures);fprintf(fid,'Final spots %d, layers %d, objective %.9g, BDT %.3f s, plan accepted %d.\n\n',o.final.plan.activeSpots,o.final.plan.activeLayers,o.final.plan.objective,o.final.plan.delivery.T_total,o.final.gate.pass);fprintf(fid,'See `dormant_transitions.csv`, `clinical_defect_components.csv`, `death_birth_mapping.csv`, and `support_trajectory.csv`.\n');
end

function [A,D]=state_matrices(model),D=model.dij.physicalDose;if iscell(D),D=D{1};end;A=state_matrix(model);end
function A=state_matrix(model),D=model.dij.physicalDose;if iscell(D),D=D{1};end;if isa(model.optiProb.BP,'matRad_DoseProjection'),A=D;elseif isa(model.optiProb.BP,'matRad_ConstantRBEProjection'),A=double(model.dij.RBE)*D;else,error('Linear exact state required.');end,end
function n=count_layers(w,groups,t),n=0;for g=1:numel(groups),n=n+any(w(groups(g).ix)>t);end,end
function x=pct(a,b),x=100*(double(a)-double(b))/max(abs(double(b)),eps);end
function x=ternary(c,a,b),if c,x=a;else,x=b;end,end
function x=solver_iterations(m),x=NaN;if isfield(m,'iterations'),x=m.iterations;elseif isfield(m,'iterationCount'),x=m.iterationCount;end,end
function s=index_list(ix),if isempty(ix),s='';else,s=strjoin(arrayfun(@(x)sprintf('%g',x),reshape(ix,1,[]),'UniformOutput',false),';');end,end
function r=empty_transition(),r=struct('event',0,'mode','','method','','incumbentSpots',0,'incumbentLayers',0,'proposedDeaths',0,'permanentDeaths',0,'rescueCandidates',0,'sameEventReturns',0,'historicalDormantBirths',0,'dormantPoolBefore',0,'dormantPoolAfter',0,'netSupportChange',0,'objectiveBefore',NaN,'objectiveRawDeletion',NaN,'objectiveAfterRepair',NaN,'worstAdversePercent',NaN,'objectivePass',false,'dosePass',false,'constraintPass',false,'solverAccepted',false,'planAccepted',false,'survivorClassification','','failureClass','','solverStatus','','solverIterations',NaN,'solverWallSeconds',NaN,'activeSpotsAfter',0,'activeLayersAfter',0,'bdtAfter',NaN,'deathGroupIds','','deathControlPoints','','deathEnergies','','wallSeconds',NaN);end
function r=empty_component(),r=struct('event',0,'phase','','method','','component','','incumbentValue',NaN,'trialValue',NaN,'change',NaN);end
function r=empty_mapping(),r=struct('event',0,'method','','deathGlobalIx','','birthGlobalIx',0,'sourceClass','','score',NaN,'finalWeight',NaN,'becamePositive',false,'cpIndex',0,'energy',NaN,'sameDeathCP',false,'sameDeathLayer',false,'doseCosine',NaN,'spatialDistance',NaN,'energyDistance',NaN);end
function r=empty_trajectory(),r=struct('event',0,'action','','supportSpots',0,'activeSpots',0,'activeLayers',0,'dormantSpots',0,'bdt',NaN,'objectiveReference',NaN);end
function r=empty_comparison(),r=struct('event',0,'proposedDeaths',0,'scalarCandidates',0,'componentCandidates',0,'overlapCount',0,'jaccard',NaN,'scalarAccepted',false,'componentAccepted',false,'scalarSameEventReturns',0,'componentSameEventReturns',0,'scalarHistoricalBirths',0,'componentHistoricalBirths',0,'scalarObjective',NaN,'componentObjective',NaN,'scalarNetChange',0,'componentNetChange',0);end
function add_paths(root),f={'config','planning','mapping','objective','delivery','reduction','evaluation','integrated','robustness'};for i=1:numel(f),addpath(fullfile(root,f{i}));end,end
function create_empty(d),if exist(d,'dir')==7,e=dir(d);e=e(~ismember({e.name},{'.','..'}));assert(isempty(e));else,mkdir(d);end,end
function m=merge_structs(a,b),m=a;f=fieldnames(b);for i=1:numel(f),k=f{i};if isfield(m,k)&&isstruct(m.(k))&&isstruct(b.(k)),m.(k)=merge_structs(m.(k),b.(k));else,m.(k)=b.(k);end,end,end
