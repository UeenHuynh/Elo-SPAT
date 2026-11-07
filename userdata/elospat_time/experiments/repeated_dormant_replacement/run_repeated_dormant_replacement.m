function output=run_repeated_dormant_replacement(outputDir)
% Repeated clinically feasible death/birth support transitions after job 22558.
narginchk(1,1);here=fileparts(mfilename('fullpath'));root=fileparts(fileparts(here));add_paths(root,here);create_empty(outputDir);
artifactFile=fullfile(root,'results','task25_27_pipeline_20deg.mat');stage1File=fullfile(root,'results','full_integrated_layer_20deg_job22180','alpha_0p3','alpha_0p3_final.mat');baseFile=fullfile(root,'results','gate_constrained_fixed_support_job22558','gate_constrained_complete.mat');
source=load(artifactFile,'cfg','machine','dij','stf','pln','map');stage1=load(stage1File,'final','finalDose','supportMask');base=load(baseFile,'output');base=base.output;
patient=load(source.cfg.patientFile,'ct','cst');cfg=merge_structs(default_config(),source.cfg);allowedFull=reshape(logical(stage1.supportMask),[],1);assert(nnz(allowedFull)==136409);
[dij,allowedToFull,~,map]=reduce_dij_to_support(source.dij,allowedFull,source.map);for i=1:numel(map),map(i).originalGlobalIx=map(i).globalIx;map(i).globalIx=i;end
[groups,~]=build_group_index(map);[groupsFull,~]=build_group_index(source.map);D=nominal(dij);model=build_matrad_dose_model(dij,patient.cst,source.pln,0);A=state_matrix(model);n=numel(allowedToFull);zeroTol=1e-10;
wRef=reshape(double(stage1.final.w(allowedToFull)),[],1);reference=evaluate_plan(wRef,model,D,patient,cfg,source,groupsFull,allowedToFull,stage1.final.w);
spec=build_clinical_gate_spec(D,patient.cst,reference.dose,2,0.2);w=reshape(double(base.w(allowedToFull)),[],1);inc=evaluate_plan(w,model,D,patient,cfg,source,groupsFull,allowedToFull,stage1.final.w);assert(clinical_accept(inc,reference,[],true).pass&&nnz(w>zeroTol)==135939);verify_gate_callbacks(D,w,spec,inc.dose);
settings=struct('maxTransitions',3,'targetAcceptedTransitions',2,'initialDeathBatch',128,'minimumDeathBatch',64,'maximumDeathBatch',512, ...
 'deathShortlist',192,'maxPerLayerFraction',0.20,'maxRescue',128,'maxLocalSurvivors',1024,'zeroTolerance',zeroTol);
rows=repmat(empty_transition(),0,1);maps=repmat(empty_mapping(),0,1);gates=repmat(empty_gate(),0,1);supportRows=repmat(empty_support(),0,1);runtime=repmat(empty_runtime(),0,1);
supportRows(end+1)=support_row(0,'accepted_job22558',w,groups,allowedFull,allowedToFull,inc,reference,8);batch=settings.initialDeathBatch;accepted=0;newBirthTotal=0;numericalFailures=0;total=tic;
for k=1:settings.maxTransitions
 wBefore=w;incumbent=evaluate_plan(w,model,D,patient,cfg,source,groupsFull,allowedToFull,stage1.final.w);supportInc=w>zeroTol;dormant=~supportInc;
 [g,comps]=objective_information(model,A,w);[gateVals,gateJ,gateState]=evaluate_clinical_gates(D,w,spec,true);
 [deathIx,deathTable,finiteCalls]=propose_deaths(w,g,comps,gateJ,spec,model,map,groups,batch,settings);if isempty(deathIx),break;end
 writetable(deathTable,fullfile(outputDir,sprintf('transition_%02d_death_proposal.csv',k)));
 survivorSupport=supportInc;survivorSupport(deathIx)=false;wDeleted=w;wDeleted(deathIx)=0;
 raw=evaluate_plan(wDeleted,model,D,patient,cfg,source,groupsFull,allowedToFull,stage1.final.w);
 [survivor,survRun]=solve_support(dij,patient.cst,source.pln,wDeleted,survivorSupport,false,reference,patient,cfg,source,groupsFull,allowedToFull,stage1.final.w,[]);
 rawGate=clinical_accept(raw,reference,[],true);survGate=clinical_accept(survivor.plan,reference,survivor.solver,survivor.solverUsable);gates=[gates;gate_rows(k,'raw_deletion',raw,reference,rawGate);gate_rows(k,'survivor',survivor.plan,reference,survGate)];runtime(end+1)=runtime_row(k,'survivor_repair',survRun,finiteCalls); %#ok<AGROW>
 save(fullfile(outputDir,sprintf('transition_%02d_survivor.mat',k)),'survivor','deathIx','survivorSupport','-v7.3');
 if survGate.pass&&survivor.plan.activeSpots<incumbent.activeSpots
  trial=survivor;method='SURVIVOR_REPAIR_SUCCESS';rescueIx=zeros(0,1);historicalBirth=zeros(0,1);
 else
  method='SURVIVOR_REPAIR_FAILURE';[rescueIx,wWarm,proposal]=gate_aware_rescue(w,survivor.w,survivorSupport,dormant,g,gateVals,gateJ,spec,D,wRef,settings);
  save(fullfile(outputDir,sprintf('transition_%02d_rescue_proposal.mat',k)),'proposal','rescueIx','-v7.3');
  trialSupport=survivorSupport;trialSupport(rescueIx)=true;[trial,gateRun]=solve_support(dij,patient.cst,source.pln,wWarm,trialSupport,true,reference,patient,cfg,source,groupsFull,allowedToFull,stage1.final.w,spec);
  trialGate=clinical_accept(trial.plan,reference,trial.solver,trial.solverUsable);gates=[gates;gate_rows(k,'gate_constrained_repair',trial.plan,reference,trialGate)]; %#ok<AGROW>
  historicalBirth=rescueIx(trial.w(rescueIx)>zeroTol);runtime(end+1)=runtime_row(k,'gate_constrained_repair',gateRun,finiteCalls); %#ok<AGROW>
 end
 if strcmp(method,'SURVIVOR_REPAIR_SUCCESS'),trialGate=survGate;end
 positive=trial.w>zeroTol;permanent=deathIx(~positive(deathIx));sameReturn=deathIx(positive(deathIx));net=nnz(positive)-nnz(supportInc);
 clinicallyAccepted=trialGate.pass&&net<0;numericallyStationary=trial.solver.numericallyStationary;
 if clinicallyAccepted
  w=trial.w;accepted=accepted+1;newBirthTotal=newBirthTotal+numel(historicalBirth);batch=min(settings.maximumDeathBatch,2*batch);action='ACCEPTED';
 else
  batch=max(settings.minimumDeathBatch,floor(batch/2));action='REJECTED';if ~trial.solverUsable,numericalFailures=numericalFailures+1;end
 end
 rows(end+1)=transition_row(k,method,action,incumbent,raw,trial.plan,deathIx,permanent,sameReturn,historicalBirth,nnz(dormant),net,trialGate,trial.solver,toc(total)); %#ok<AGROW>
 maps=[maps;mapping_rows(k,deathIx,historicalBirth,map,groups,wBefore,trial.w)]; %#ok<AGROW>
 current=evaluate_plan(w,model,D,patient,cfg,source,groupsFull,allowedToFull,stage1.final.w);supportRows(end+1)=support_row(k,action,w,groups,allowedFull,allowedToFull,current,reference,numel(historicalBirth)); %#ok<AGROW>
 flush(outputDir,rows,maps,gates,supportRows,runtime);checkpoint=struct('transition',k,'w',w,'batchNext',batch,'accepted',accepted,'newBirthTotal',newBirthTotal,'row',rows(end));save(fullfile(outputDir,sprintf('transition_%02d_checkpoint.mat',k)),'checkpoint','-v7.3');
 fprintf('REPEATED k=%d method=%s action=%s deaths=%d permanent=%d births=%d net=%d obj=%.6f gates=%d stationary=%d\n',k,method,action,numel(deathIx),numel(permanent),numel(historicalBirth),net,trial.plan.objective,trialGate.dosePass,numericallyStationary);
 if accepted>=settings.targetAcceptedTransitions,break;end
end
final=evaluate_plan(w,model,D,patient,cfg,source,groupsFull,allowedToFull,stage1.final.w);finalGate=clinical_accept(final,reference,[],true);cumulative=135939-final.activeSpots;
if accepted>=2&&cumulative>0&&newBirthTotal>0&&finalGate.pass,status='REPEATED_DORMANT_REPLACEMENT_SUPPORTED';elseif numericalFailures>0&&accepted<2,status='INCONCLUSIVE';else,status='REPEATED_DORMANT_REPLACEMENT_NOT_SUPPORTED';end
output=struct('status',status,'settings',settings,'w',w,'reference',reference,'final',final,'finalGate',finalGate,'rows',rows,'mapping',maps,'gates',gates,'supportTrajectory',supportRows,'runtime',runtime,'acceptedTransitions',accepted,'newHistoricalBirths',newBirthTotal,'cumulativeReductionFromJob22558',cumulative,'wallSeconds',toc(total),'baseFile',baseFile);
save(fullfile(outputDir,'repeated_dormant_replacement_complete.mat'),'output','-v7.3');flush(outputDir,rows,maps,gates,supportRows,runtime);write_report(outputDir,output);fprintf('FINAL_STATUS: %s\n',status);
end

function [sol,run]=solve_support(dij,cst,pln,wStart,support,useGates,reference,patient,cfg,source,groupsFull,toFull,wTemplate,spec)
[dijR,r2f]=reduce_dij_to_support(dij,support);m=build_matrad_dose_model(dijR,cst,pln,0);w0=wStart(r2f);D=nominal(dijR);
set=struct('ipopt',struct('max_iter',260,'max_cpu_time',120000,'print_level',3,'acceptable_iter',8,'acceptable_obj_change_tol',1e-5));
if useGates,sp=build_clinical_gate_spec(D,cst,reference.dose,2,0.2);run=solve_clinical_gate_ipopt(m,w0,D,sp,set);else,run=solve_weighted_matrad_ipopt(m,w0,zeros(size(w0)),set);run.solverUsable=run.postKKT.primalInfeasibility<=1e-6&&any(run.status==[0 1 -1 -4]);run.numericallyStationary=run.desiredConverged;end
w=zeros(double(dij.totalNumOfBixels),1);w(r2f)=run.w;plan=evaluate_plan(w,build_matrad_dose_model(dij,cst,pln,0),nominal(dij),patient,cfg,source,groupsFull,toFull,wTemplate);sol=struct('w',w,'plan',plan,'solver',run,'solverUsable',run.solverUsable,'support',support);
end

function [death,tbl,calls]=propose_deaths(w,g,comps,gateJ,spec,model,map,groups,batch,s)
active=w>s.zeroTolerance;n=numel(w);obj=abs(g.*w);comp=zeros(n,1);for i=1:numel(comps),comp=comp+(w.*comps(i).gradient/max(1,abs(comps(i).value))).^2;end;comp=sqrt(comp);
target=zeros(n,1);relief=zeros(n,1);for i=1:numel(spec),q=abs(gateJ(i,:)'.*w);if strcmp(spec(i).sense,'lower'),target=max(target,q);else,relief=max(relief,q);end,end
base=rank01(obj)+rank01(comp)+rank01(target)-0.25*rank01(relief);base(~active)=Inf;[~,o]=sort(base,'ascend');short=o(1:min(s.deathShortlist,nnz(active)));f0=evaluate_matrad_dose_model(model,w,false);finite=inf(n,1);calls=1;
for q=1:numel(short),z=w;z(short(q))=0;finite(short(q))=evaluate_matrad_dose_model(model,z,false)-f0;calls=calls+1;end
score=rank_subset(finite,short)+rank_subset(comp,short)+rank_subset(target,short)-0.25*rank_subset(relief,short);[~,q]=sort(score(short),'ascend');order=short(q);
groupOf=zeros(n,1);for i=1:numel(groups),groupOf(groups(i).ix)=i;end;cp=[map.cpIndex]';layerCap=max(2,ceil(s.maxPerLayerFraction*batch));cpCap=max(4,ceil(batch/4));lc=zeros(numel(groups),1);cc=zeros(max(cp),1);death=[];
for i=reshape(order,1,[]),gi=groupOf(i);ci=cp(i);if lc(gi)>=layerCap||cc(ci)>=cpCap,continue;end;death(end+1,1)=i;lc(gi)=lc(gi)+1;cc(ci)=cc(ci)+1;if numel(death)>=batch,break;end,end
tbl=table(death,reshape([map(death).originalGlobalIx],[],1),w(death),finite(death),obj(death),comp(death),target(death),relief(death),groupOf(death),cp(death),'VariableNames',{'allowedIx','globalIx','weight','finiteRemovalObjectiveChange','objectiveSensitivity','componentSensitivity','targetContribution','gateRelief','layerId','controlPoint'});
end

function [rescue,wWarm,p]=gate_aware_rescue(wInc,wFail,support,dormant,gInc,~,~,spec,D,wRef,s)
[vals,J,rows]=evaluate_clinical_gates(D,wFail,spec,true);m=numel(vals);A=zeros(m,size(D,2));b=zeros(m,1);for i=1:m,if strcmp(spec(i).sense,'lower'),A(i,:)=-J(i,:);b(i)=vals(i)-spec(i).limit;else,A(i,:)=J(i,:);b(i)=spec(i).limit-vals(i);end,end;A=sparse(A);
importance=max(abs(A),[],1)';active=support&(wFail>s.zeroTolerance);surv=top(importance.*max(wFail,1e-8),active,s.maxLocalSurvivors);dorm=find(dormant);Jix=unique([surv;dorm],'stable');isDorm=ismember(Jix,dorm);
Aj=A(:,Jix);scale=max(0.05,abs(b));lb=-wFail(Jix);lb(isDorm)=0;ub=max([abs(wFail(Jix)),abs(wRef(Jix)),ones(numel(Jix),1)],[],2);ub(isDorm)=max(2*abs(wRef(Jix(isDorm))),1);
lp=optimoptions('linprog','Algorithm','dual-simplex','Display','none');[x,~,ef]=linprog([zeros(numel(Jix),1);1],[Aj,-scale],b,[],[],[lb;0],[ub;10],lp);rho=Inf;if ef>0,rho=x(end);end
rescue=zeros(0,1);wWarm=wFail;if rho<=1e-6
 H=spdiags(1./max(ub,1e-8).^2,0,numel(Jix),numel(Jix));eta=max(1e-8,1e-3*max(abs(gInc(Jix).*ub)));qp=optimoptions('quadprog','Algorithm','interior-point-convex','Display','none');[d,~,qe]=quadprog(eta*H,gInc(Jix),Aj,b,[],[],lb,ub,[],qp);
 if qe>0,cand=find(isDorm&d>1e-8);[~,o]=sort(d(cand),'descend');cand=cand(o(1:min(s.maxRescue,numel(o))));rescue=Jix(cand);keep=ismember(Jix,rescue);lb2=lb;ub2=ub;lb2(isDorm&~keep)=0;ub2(isDorm&~keep)=0;[d,~,qe]=quadprog(eta*H,gInc(Jix),Aj,b,[],[],lb2,ub2,[],qp);if qe>0,wWarm(Jix)=max(0,wFail(Jix)+d);else,rescue=zeros(0,1);end,end
end
p=struct('rho',rho,'localIx',Jix,'rescueIx',rescue,'gateRows',rows,'gateValues',vals,'A',Aj,'b',b);
end

function [g,c]=objective_information(model,A,w)
[~,g]=evaluate_matrad_dose_model(model,w,true);x=A*w;c=repmat(struct('label','','value',NaN,'gradient',zeros(numel(w),1)),0,1);
for i=1:size(model.cst,1),vox=reshape(double(model.cst{i,4}{1}),[],1);for j=1:numel(model.cst{i,6}),o=model.cst{i,6}{j};if ~isa(o,'DoseObjectives.matRad_DoseObjective'),continue;end;o=model.optiProb.BP.setBiologicalDosePrescriptions(o,model.cst{i,5}.alphaX,model.cst{i,5}.betaX);dg=zeros(size(x));dg(vox)=o.penalty*o.computeDoseObjectiveGradient(x(vox));c(end+1)=struct('label',sprintf('%s:%s',model.cst{i,2},class(o)),'value',double(o.penalty*o.computeDoseObjectiveFunction(x(vox))),'gradient',A'*dg);end,end
end
function p=evaluate_plan(w,model,D,patient,cfg,source,groupsFull,toFull,wTemplate),[f,~,d]=evaluate_matrad_dose_model(model,w,false);z=zeros(size(wTemplate));z(toFull)=w;p=struct('objective',f,'maxConstraintViolation',d.maxConstraintViolation,'dose',compute_dvh_metrics(full(D*w),patient.ct,patient.cst,cfg),'delivery',compute_delivery_metrics(z,source.map,groupsFull,cfg,source.machine),'activeSpots',nnz(w>1e-10),'activeLayers',count_layers(w,build_groups(source.map,toFull),1e-10));end
function gate=clinical_accept(p,r,solver,usable),[pass,worst,detail]=dose_quality(p.dose,r.dose);op=p.objective<=1.02*r.objective;cp=p.maxConstraintViolation<=1e-6;gate=struct('pass',usable&&pass&&op&&cp,'dosePass',pass,'objectivePass',op,'constraintPass',cp,'worstAdverse',worst,'detail',{detail});if isempty(solver),gate.numericallyStationary=false;else,gate.numericallyStationary=solver.numericallyStationary;end,end
function [pass,worst,rows]=dose_quality(c,r),pass=true;worst=0;rows={};for i=1:numel(c.structures),a=r.structures(i);b=c.structures(i);if a.isTarget,keys={'D95','D98'};sense='lower';elseif strcmpi(a.type,'OAR')&&~ismember(upper(a.name),{'BODY','EXTERNAL'}),keys={'Dmean','Dmax'};sense='upper';else,continue;end;for k=1:numel(keys),m=keys{k};if strcmp(sense,'lower'),lim=.98*a.(m);ok=b.(m)>=lim;adv=max(0,100*(a.(m)-b.(m))/max(abs(a.(m)),eps));else,lim=a.(m)+max(.2,.02*max(abs(a.(m)),1));ok=b.(m)<=lim;adv=max(0,100*(b.(m)-a.(m))/max(abs(a.(m)),eps));end;rows(end+1,:)={a.name,m,b.(m),lim,ok};pass=pass&&ok;worst=max(worst,adv);end,end,end
function r=transition_row(k,method,action,inc,raw,trial,deaths,permanent,returns,births,dormant,net,gate,solver,wall),r=empty_transition();r.iteration=k;r.method=method;r.action=action;r.incumbentSpots=inc.activeSpots;r.incumbentLayers=inc.activeLayers;r.proposedDeaths=numel(deaths);r.permanentDeaths=numel(permanent);r.sameEventReturns=numel(returns);r.historicalDormantBirths=numel(births);r.dormantPoolBefore=dormant;r.netSupportChange=net;r.objectiveBefore=inc.objective;r.objectiveRaw=raw.objective;r.objectiveAfter=trial.objective;r.activeSpotsAfter=trial.activeSpots;r.activeLayersAfter=trial.activeLayers;r.bdtBefore=inc.delivery.T_total;r.bdtAfter=trial.delivery.T_total;r.dosePass=gate.dosePass;r.objectivePass=gate.objectivePass;r.clinicallyAcceptedSupport=gate.pass;r.numericallyStationary=gate.numericallyStationary;r.solverStatus=solver.statusMessage;r.primalInfeasibility=readfield(solver,'primalInfeasibility',readnested(solver,'postKKT','primalInfeasibility'));r.dualInfeasibility=readfield(solver,'dualInfeasibility',readnested(solver,'postKKT','dualInfeasibility'));r.wallSeconds=wall;end
function r=empty_transition(),r=struct('iteration',0,'method','','action','','incumbentSpots',0,'incumbentLayers',0,'proposedDeaths',0,'permanentDeaths',0,'sameEventReturns',0,'historicalDormantBirths',0,'dormantPoolBefore',0,'netSupportChange',0,'objectiveBefore',NaN,'objectiveRaw',NaN,'objectiveAfter',NaN,'activeSpotsAfter',0,'activeLayersAfter',0,'bdtBefore',NaN,'bdtAfter',NaN,'dosePass',false,'objectivePass',false,'clinicallyAcceptedSupport',false,'numericallyStationary',false,'solverStatus','','primalInfeasibility',NaN,'dualInfeasibility',NaN,'wallSeconds',NaN);end
function r=empty_mapping(),r=struct('iteration',0,'deathGlobalIx','','birthGlobalIx',0,'birthLayer',0,'birthControlPoint',0,'weightBefore',0,'weightAfter',0);end
function r=mapping_rows(k,d,b,map,groups,w0,w1),r=repmat(empty_mapping(),numel(b),1);go=zeros(numel(map),1);for i=1:numel(groups),go(groups(i).ix)=i;end;ds=index_list([map(d).originalGlobalIx]);for i=1:numel(b),q=b(i);r(i)=struct('iteration',k,'deathGlobalIx',ds,'birthGlobalIx',map(q).originalGlobalIx,'birthLayer',go(q),'birthControlPoint',map(q).cpIndex,'weightBefore',w0(q),'weightAfter',w1(q));end,end
function r=empty_gate(),r=struct('iteration',0,'phase','','structure','','metric','','value',NaN,'limit',NaN,'pass',false);end
function r=gate_rows(k,phase,p,ref,gate),d=gate.detail;r=repmat(empty_gate(),size(d,1),1);for i=1:size(d,1),r(i)=struct('iteration',k,'phase',phase,'structure',d{i,1},'metric',d{i,2},'value',d{i,3},'limit',d{i,4},'pass',d{i,5});end,end
function r=empty_support(),r=struct('iteration',0,'action','','activeSpots',0,'activeLayers',0,'dormantPool',0,'objective',NaN,'bdt',NaN,'objectiveChangePercent',NaN,'bdtChangePercent',NaN,'historicalBirthsThisIteration',0);end
function r=support_row(k,a,w,g,allowed,toFull,p,ref,b),r=empty_support();r.iteration=k;r.action=a;r.activeSpots=nnz(w>1e-10);r.activeLayers=count_layers(w,g,1e-10);r.dormantPool=nnz(allowed)-r.activeSpots;r.objective=p.objective;r.bdt=p.delivery.T_total;r.objectiveChangePercent=100*(p.objective-ref.objective)/ref.objective;r.bdtChangePercent=100*(p.delivery.T_total-ref.delivery.T_total)/ref.delivery.T_total;r.historicalBirthsThisIteration=b;end
function r=empty_runtime(),r=struct('iteration',0,'phase','','solverWallSeconds',NaN,'iterations',NaN,'objectiveCalls',0,'gradientCalls',0,'constraintCalls',0,'jacobianCalls',0,'finiteRemovalObjectiveCalls',0);end
function r=runtime_row(k,p,s,c),r=empty_runtime();r.iteration=k;r.phase=p;r.solverWallSeconds=s.wallSeconds;r.iterations=s.iterations;r.objectiveCalls=s.callbackCounts.objective;r.gradientCalls=s.callbackCounts.gradient;r.constraintCalls=readfield(s.callbackCounts,'constraints',0);r.jacobianCalls=readfield(s.callbackCounts,'jacobian',0);r.finiteRemovalObjectiveCalls=c;end
function flush(o,a,b,c,d,e),if ~isempty(a),writetable(struct2table(a,'AsArray',true),fullfile(o,'iteration_transitions.csv'));end;if ~isempty(b),writetable(struct2table(b,'AsArray',true),fullfile(o,'death_birth_mapping.csv'));end;if ~isempty(c),writetable(struct2table(c,'AsArray',true),fullfile(o,'clinical_gate_trajectory.csv'));end;if ~isempty(d),writetable(struct2table(d,'AsArray',true),fullfile(o,'support_bdt_trajectory.csv'));end;if ~isempty(e),writetable(struct2table(e,'AsArray',true),fullfile(o,'runtime_evaluations.csv'));end;end
function write_report(o,x),f=fopen(fullfile(o,'REPEATED_DORMANT_REPLACEMENT_REPORT.md'),'w');c=onCleanup(@()fclose(f));fprintf(f,'# Repeated dormant replacement\n\nStatus: **%s**\n\nAccepted transitions: %d; new historical births: %d; cumulative reduction from job 22558: %d spots.\n\nFinal objective %.9g; BDT %.3f s; all gates pass=%d.\n',x.status,x.acceptedTransitions,x.newHistoricalBirths,x.cumulativeReductionFromJob22558,x.final.objective,x.final.delivery.T_total,x.finalGate.pass);end
function verify_gate_callbacks(D,w,spec,dose)
[v,J]=evaluate_clinical_gates(D,w,spec,true);assert(all(isfinite(v))&&all(size(J)==[numel(spec),numel(w)]));
for i=1:numel(spec),q=find(strcmp({dose.structures.name},spec(i).structure),1);assert(~isempty(q));assert(abs(v(i)-dose.structures(q).(spec(i).metric))<=1e-9*max(1,abs(v(i))));end
dir=sin((1:numel(w))')/sqrt(numel(w));h=1e-7;vp=evaluate_clinical_gates(D,w+h*dir,spec,false);vm=evaluate_clinical_gates(D,w-h*dir,spec,false);fd=(vp-vm)/(2*h);rel=abs(fd-J*dir)./max(1,abs(fd));assert(max(rel)<1e-4);
end
function x=rank01(v),x=zeros(size(v));q=isfinite(v);if any(q),a=min(v(q));b=max(v(q));x(q)=(v(q)-a)/max(b-a,eps);end;end
function x=rank_subset(v,ix),x=zeros(size(v));x(ix)=rank01(v(ix));end
function ix=top(score,mask,k),q=find(mask&isfinite(score));[~,o]=sort(score(q),'descend');ix=q(o(1:min(k,numel(o))));end
function D=nominal(d),D=d.physicalDose;if iscell(D),D=D{1};end,end
function A=state_matrix(m),A=nominal(m.dij);if isa(m.optiProb.BP,'matRad_ConstantRBEProjection'),A=double(m.dij.RBE)*A;end,end
function g=build_groups(map,toFull),m=map(toFull);for i=1:numel(m),m(i).globalIx=i;end;g=build_group_index(m);end
function n=count_layers(w,g,t),n=0;for i=1:numel(g),n=n+any(w(g(i).ix)>t);end,end
function s=index_list(ix),s=strjoin(arrayfun(@num2str,reshape(ix,1,[]),'UniformOutput',false),';');end
function v=readfield(s,n,d),if isstruct(s)&&isfield(s,n),v=s.(n);else,v=d;end,end
function v=readnested(s,a,b),v=NaN;if isstruct(s)&&isfield(s,a)&&isfield(s.(a),b),v=s.(a).(b);end,end
function add_paths(root,here),f={'config','planning','mapping','objective','delivery','reduction','evaluation','integrated','robustness'};for i=1:numel(f),addpath(fullfile(root,f{i}));end;addpath(here);end
function create_empty(d),if exist(d,'dir')==7,e=dir(d);e=e(~ismember({e.name},{'.','..'}));assert(isempty(e));else,mkdir(d);end,end
function m=merge_structs(a,b),m=a;f=fieldnames(b);for i=1:numel(f),k=f{i};if isfield(m,k)&&isstruct(m.(k))&&isstruct(b.(k)),m.(k)=merge_structs(m.(k),b.(k));else,m.(k)=b.(k);end,end,end
