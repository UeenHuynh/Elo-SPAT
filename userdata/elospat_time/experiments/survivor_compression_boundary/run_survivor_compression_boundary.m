function output=run_survivor_compression_boundary(outputDir)
% Bounded search for the survivor-compression boundary, then one dormant rescue.
narginchk(1,1);here=fileparts(mfilename('fullpath'));root=fileparts(fileparts(here));add_paths(root,here);create_empty(outputDir);
artifactFile=fullfile(root,'results','task25_27_pipeline_20deg.mat');stage1File=fullfile(root,'results','full_integrated_layer_20deg_job22180','alpha_0p3','alpha_0p3_final.mat');baseFile=fullfile(root,'results','repeated_dormant_replacement_job22761','repeated_dormant_replacement_complete.mat');
source=load(artifactFile,'cfg','machine','dij','stf','pln','map');stage1=load(stage1File,'final','finalDose','supportMask');base=load(baseFile,'output');base=base.output;
patient=load(source.cfg.patientFile,'ct','cst');cfg=merge_structs(default_config(),source.cfg);allowedFull=reshape(logical(stage1.supportMask),[],1);assert(nnz(allowedFull)==136409);
[dij,toFull,~,map]=reduce_dij_to_support(source.dij,allowedFull,source.map);for i=1:numel(map),map(i).originalGlobalIx=map(i).globalIx;map(i).globalIx=i;end
[groups,~]=build_group_index(map);[groupsFull,~]=build_group_index(source.map);D=nominal(dij);model=build_matrad_dose_model(dij,patient.cst,source.pln,0);A=state_matrix(model);zeroTol=1e-10;
wRef=reshape(double(stage1.final.w(toFull)),[],1);reference=evaluate_plan(wRef,model,D,patient,cfg,source,groupsFull,toFull,stage1.final.w);
w=reshape(double(base.w),[],1);assert(numel(w)==numel(toFull));inc=evaluate_plan(w,model,D,patient,cfg,source,groupsFull,toFull,stage1.final.w);assert(inc.activeSpots==135619&&clinical_accept(inc,reference,[],true).pass);
spec=build_clinical_gate_spec(D,patient.cst,reference.dose,2,0.2);
settings=struct('schedule',[256 384 512 677],'candidateBlocks',10,'confirmationBlocks',1,'maxRescue',256,'maxLocalSurvivors',1400,'zeroTolerance',zeroTol);
rows=repmat(empty_attempt(),0,1);gates=repmat(empty_gate(),0,1);runtime=repmat(empty_runtime(),0,1);maps=repmat(empty_mapping(),0,1);boundaryFound=false;boundaryReproduced=false;rescueStatus='NOT_RUN';stationarity=false;total=tic;
for level=1:numel(settings.schedule)
 batch=settings.schedule(level);wBefore=w;incumbent=evaluate_plan(w,model,D,patient,cfg,source,groupsFull,toFull,stage1.final.w);supportInc=w>zeroTol;historicalDormant=~supportInc;
 [g,comps]=objective_information(model,A,w);[~,gateJ]=evaluate_clinical_gates(D,w,spec,true);
 [deathIx,proposal,proposalPool]=propose_coherent_deaths(w,g,comps,gateJ,spec,model,map,groups,batch,settings,1);
 writetable(proposal,fullfile(outputDir,sprintf('level_%02d_death_proposal.csv',level)));
 [raw,survivor,survRun,survGate]=run_survivor(deathIx,w,supportInc,dij,patient,cfg,source,groupsFull,toFull,stage1.final.w,reference);
 gates=[gates;gate_rows(level,1,'raw_deletion',raw,reference);gate_rows(level,1,'survivor',survivor.plan,reference)];runtime(end+1)=runtime_row(level,1,'survivor',survRun,proposalPool.finiteCalls); %#ok<AGROW>
 survivorFailure=~survGate.pass;rows(end+1)=attempt_row(level,1,batch,'PRIMARY',incumbent,raw,survivor.plan,deathIx,survivorFailure,survGate,survRun,toc(total)); %#ok<AGROW>
 save(fullfile(outputDir,sprintf('level_%02d_primary.mat',level)),'deathIx','raw','survivor','survGate','proposalPool','-v7.3');flush(outputDir,rows,gates,runtime,maps);
 fprintf('BOUNDARY level=%d batch=%d survivorPass=%d obj=%.6f BDT=%.3f stationary=%d\n',level,batch,survGate.pass,survivor.plan.objective,survivor.plan.delivery.T_total,survRun.numericallyStationary);
 if ~survivorFailure
  w=survivor.w;stationarity=survRun.numericallyStationary;continue;
 end
 % One matched, independent coherent block at the same pressure establishes reproducibility.
 [death2,proposal2,pool2]=propose_coherent_deaths(w,g,comps,gateJ,spec,model,map,groups,batch,settings,2);
 writetable(proposal2,fullfile(outputDir,sprintf('level_%02d_confirmation_death_proposal.csv',level)));
 [raw2,survivor2,run2,gate2]=run_survivor(death2,w,supportInc,dij,patient,cfg,source,groupsFull,toFull,stage1.final.w,reference);
 gates=[gates;gate_rows(level,2,'raw_deletion',raw2,reference);gate_rows(level,2,'survivor',survivor2.plan,reference)];runtime(end+1)=runtime_row(level,2,'survivor_confirmation',run2,pool2.finiteCalls); %#ok<AGROW>
 rows(end+1)=attempt_row(level,2,batch,'CONFIRMATION',incumbent,raw2,survivor2.plan,death2,~gate2.pass,gate2,run2,toc(total)); %#ok<AGROW>
 boundaryFound=true;boundaryReproduced=~gate2.pass;
 % Rescue the primary failed event only. Historical dormant and same-event returns stay distinct.
 [rescueIx,wWarm,rescueProposal]=gate_aware_rescue(w,survivor.w,survivor.support,historicalDormant,deathIx,g,spec,D,wRef,settings);
 save(fullfile(outputDir,'boundary_rescue_proposal.mat'),'rescueProposal','rescueIx','-v7.3');
 trialSupport=survivor.support;trialSupport(rescueIx)=true;
 [trial,gateRun]=solve_support(dij,patient.cst,source.pln,wWarm,trialSupport,true,reference,patient,cfg,source,groupsFull,toFull,stage1.final.w);
 trialGate=clinical_accept(trial.plan,reference,gateRun,gateRun.solverUsable);historicalBirth=rescueIx(historicalDormant(rescueIx)&trial.w(rescueIx)>zeroTol);sameReturn=deathIx(trial.w(deathIx)>zeroTol);permanent=deathIx(trial.w(deathIx)<=zeroTol);net=trial.plan.activeSpots-incumbent.activeSpots;
 excludedPreserved=all(trial.w>=-zeroTol);accepted=trialGate.pass&&~isempty(historicalBirth)&&net<0&&excludedPreserved&&numel(sameReturn)<0.9*numel(deathIx);
 if accepted,rescueStatus='SUPPORTED';w=trial.w;else,rescueStatus='NOT_SUPPORTED';end
 stationarity=gateRun.numericallyStationary;gates=[gates;gate_rows(level,3,'gate_constrained_rescue',trial.plan,reference)];runtime(end+1)=runtime_row(level,3,'gate_constrained_rescue',gateRun,0); %#ok<AGROW>
 rows(end+1)=rescue_row(level,batch,incumbent,raw,trial.plan,deathIx,permanent,sameReturn,historicalBirth,trialGate,gateRun,toc(total)); %#ok<AGROW>
 maps=mapping_rows(level,deathIx,historicalBirth,map,groups,wBefore,trial.w); %#ok<NASGU>
 save(fullfile(outputDir,'boundary_rescue_result.mat'),'trial','trialGate','deathIx','historicalBirth','sameReturn','permanent','accepted','-v7.3');flush(outputDir,rows,gates,runtime,maps);break;
end
final=evaluate_plan(w,model,D,patient,cfg,source,groupsFull,toFull,stage1.final.w);finalGate=clinical_accept(final,reference,[],true);
if boundaryFound,boundaryLabel='SURVIVOR_COMPRESSION_BOUNDARY_FOUND';else,boundaryLabel='SURVIVOR_COMPRESSION_BOUNDARY_NOT_FOUND';end
if strcmp(rescueStatus,'SUPPORTED'),dormantLabel='REPEATED_DORMANT_REPLACEMENT_SUPPORTED';elseif boundaryFound&&strcmp(rescueStatus,'NOT_SUPPORTED'),dormantLabel='REPEATED_DORMANT_REPLACEMENT_NOT_SUPPORTED';else,dormantLabel='INCONCLUSIVE';end
if stationarity,stationarityLabel='NUMERICALLY_STATIONARY';else,stationarityLabel='NOT_ESTABLISHED';end
output=struct('boundaryLabel',boundaryLabel,'boundaryReproduced',boundaryReproduced,'dormantLabel',dormantLabel,'stationarityLabel',stationarityLabel,'settings',settings,'w',w,'final',final,'finalGate',finalGate,'attempts',rows,'gates',gates,'runtime',runtime,'mapping',maps,'wallSeconds',toc(total));
save(fullfile(outputDir,'survivor_compression_boundary_complete.mat'),'output','-v7.3');flush(outputDir,rows,gates,runtime,maps);write_report(outputDir,output);fprintf('FINAL_LABELS: %s | %s | %s\n',boundaryLabel,dormantLabel,stationarityLabel);
end

function [raw,sol,run,gate]=run_survivor(deaths,w,support,dij,patient,cfg,source,groupsFull,toFull,wTemplate,reference)
s=support;s(deaths)=false;z=w;z(deaths)=0;model=build_matrad_dose_model(dij,patient.cst,source.pln,0);raw=evaluate_plan(z,model,nominal(dij),patient,cfg,source,groupsFull,toFull,wTemplate);
[sol,run]=solve_support(dij,patient.cst,source.pln,z,s,false,reference,patient,cfg,source,groupsFull,toFull,wTemplate);gate=clinical_accept(sol.plan,reference,run,run.solverUsable);
end

function [death,tbl,pool]=propose_coherent_deaths(w,g,comps,gateJ,spec,model,map,groups,batch,s,rankChoice)
n=numel(w);active=w>s.zeroTolerance;comp=zeros(n,1);for i=1:numel(comps),comp=comp+(w.*comps(i).gradient/max(1,abs(comps(i).value))).^2;end;comp=sqrt(comp);target=zeros(n,1);relief=zeros(n,1);
for i=1:numel(spec),q=abs(gateJ(i,:)'.*w);if strcmp(spec(i).sense,'lower'),target=max(target,q);else,relief=max(relief,q);end,end
spotScore=rank01(abs(g.*w))+rank01(comp)+rank01(target)-0.25*rank01(relief);groupScore=inf(numel(groups),1);cp=zeros(numel(groups),1);en=zeros(numel(groups),1);
for j=1:numel(groups),ix=groups(j).ix;ix=ix(active(ix));if ~isempty(ix),groupScore(j)=median(spotScore(ix));cp(j)=map(ix(1)).cpIndex;en(j)=map(ix(1)).energyIndex;end,end
[~,order]=sort(groupScore,'ascend');candidates=cell(0,1);seed=zeros(0,1);for q=reshape(order,1,[]),if ~isfinite(groupScore(q)),continue;end;same=find(cp==cp(q)&isfinite(groupScore));[~,o]=sortrows([abs(en(same)-en(q)),groupScore(same)],[1 2]);ix=[];for h=reshape(same(o),1,[]),local=groups(h).ix;local=local(active(local));[~,so]=sort(spotScore(local),'ascend');ix=[ix;reshape(local(so),[],1)];if numel(ix)>=batch,break;end,end;if numel(ix)>=batch,candidates{end+1,1}=ix(1:batch);seed(end+1,1)=q;end;if numel(candidates)>=s.candidateBlocks,break;end,end
assert(numel(candidates)>=rankChoice,'Not enough coherent deletion blocks');f0=evaluate_matrad_dose_model(model,w,false);quality=zeros(numel(candidates),3);rawObj=zeros(numel(candidates),1);
for q=1:numel(candidates),z=w;z(candidates{q})=0;rawObj(q)=evaluate_matrad_dose_model(model,z,false);quality(q,:)=[max(0,rawObj(q)-f0),mean(target(candidates{q})),mean(comp(candidates{q}))];end
[~,o]=sortrows(quality,[1 2 3]);pick=o(rankChoice);death=candidates{pick};pool=struct('finiteCalls',numel(candidates)+1,'seedGroup',seed(pick),'quality',quality,'rankChoice',rankChoice);
tbl=table(death,reshape([map(death).originalGlobalIx],[],1),w(death),spotScore(death),abs(g(death).*w(death)),comp(death),target(death),relief(death),repmat(seed(pick),numel(death),1),'VariableNames',{'allowedIx','globalIx','weight','combinedScore','objectiveSensitivity','componentSensitivity','targetContribution','gateRelief','seedLayer'});
end

function [rescue,wWarm,p]=gate_aware_rescue(wInc,wFail,support,historical,deaths,gInc,spec,D,wRef,s)
[vals,J,rows]=evaluate_clinical_gates(D,wFail,spec,true);m=numel(vals);A=zeros(m,size(D,2));b=zeros(m,1);for i=1:m,if strcmp(spec(i).sense,'lower'),A(i,:)=-J(i,:);b(i)=vals(i)-spec(i).limit;else,A(i,:)=J(i,:);b(i)=spec(i).limit-vals(i);end,end;A=sparse(A);
importance=max(abs(A),[],1)';surv=top(importance.*max(wFail,1e-8),support&(wFail>s.zeroTolerance),s.maxLocalSurvivors);eligible=find(historical);eligible=unique([eligible;deaths(:)],'stable');Jix=unique([surv;eligible],'stable');isEligible=ismember(Jix,eligible);
Aj=A(:,Jix);scale=max(0.05,abs(b));lb=-wFail(Jix);lb(isEligible)=0;ub=max([abs(wFail(Jix)),abs(wRef(Jix)),ones(numel(Jix),1)],[],2);ub(isEligible)=max(2*abs(wRef(Jix(isEligible))),1);
lp=optimoptions('linprog','Algorithm','dual-simplex','Display','none');[x,~,ef]=linprog([zeros(numel(Jix),1);1],[Aj,-scale],b,[],[],[lb;0],[ub;10],lp);rho=Inf;if ef>0,rho=x(end);end
rescue=zeros(0,1);wWarm=wFail;if rho<=1e-6
 H=spdiags(1./max(ub,1e-8).^2,0,numel(Jix),numel(Jix));eta=max(1e-8,1e-3*max(abs(gInc(Jix).*ub)));qp=optimoptions('quadprog','Algorithm','interior-point-convex','Display','none');[d,~,qe]=quadprog(eta*H,gInc(Jix),Aj,b,[],[],lb,ub,[],qp);
 if qe>0,cand=find(isEligible&d>1e-8);[~,o]=sort(d(cand),'descend');cand=cand(o(1:min(s.maxRescue,numel(o))));rescue=Jix(cand);keep=ismember(Jix,rescue);lb2=lb;ub2=ub;lb2(isEligible&~keep)=0;ub2(isEligible&~keep)=0;[d,~,qe]=quadprog(eta*H,gInc(Jix),Aj,b,[],[],lb2,ub2,[],qp);if qe>0,wWarm(Jix)=max(0,wFail(Jix)+d);else,rescue=zeros(0,1);end,end
end
p=struct('rho',rho,'localIx',Jix,'eligibleIx',eligible,'rescueIx',rescue,'gateRows',rows,'gateValues',vals,'A',Aj,'b',b);
end

function [sol,run]=solve_support(dij,cst,pln,wStart,support,useGates,reference,patient,cfg,source,groupsFull,toFull,wTemplate)
[dijR,r2f]=reduce_dij_to_support(dij,support);m=build_matrad_dose_model(dijR,cst,pln,0);w0=wStart(r2f);D=nominal(dijR);set=struct('ipopt',struct('max_iter',260,'max_cpu_time',120000,'print_level',3,'acceptable_iter',8,'acceptable_obj_change_tol',1e-5));
if useGates,sp=build_clinical_gate_spec(D,cst,reference.dose,2,0.2);run=solve_clinical_gate_ipopt(m,w0,D,sp,set);else,run=solve_weighted_matrad_ipopt(m,w0,zeros(size(w0)),set);run.solverUsable=run.postKKT.primalInfeasibility<=1e-6&&any(run.status==[0 1 -1 -4]);run.numericallyStationary=run.desiredConverged;end
w=zeros(double(dij.totalNumOfBixels),1);w(r2f)=run.w;plan=evaluate_plan(w,build_matrad_dose_model(dij,cst,pln,0),nominal(dij),patient,cfg,source,groupsFull,toFull,wTemplate);sol=struct('w',w,'plan',plan,'support',support);
end

function p=evaluate_plan(w,model,D,patient,cfg,source,groupsFull,toFull,wTemplate),[f,~,d]=evaluate_matrad_dose_model(model,w,false);z=zeros(size(wTemplate));z(toFull)=w;p=struct('objective',f,'maxConstraintViolation',d.maxConstraintViolation,'dose',compute_dvh_metrics(full(D*w),patient.ct,patient.cst,cfg),'delivery',compute_delivery_metrics(z,source.map,groupsFull,cfg,source.machine),'activeSpots',nnz(w>1e-10),'activeLayers',count_layers(w,build_groups(source.map,toFull),1e-10));end
function gate=clinical_accept(p,r,solver,usable),[pass,worst,detail]=dose_quality(p.dose,r.dose);op=p.objective<=1.02*r.objective;cp=p.maxConstraintViolation<=1e-6;gate=struct('pass',usable&&pass&&op&&cp,'dosePass',pass,'objectivePass',op,'constraintPass',cp,'worstAdverse',worst,'detail',{detail});if isempty(solver),gate.numericallyStationary=false;else,gate.numericallyStationary=solver.numericallyStationary;end,end
function [pass,worst,rows]=dose_quality(c,r),pass=true;worst=0;rows={};for i=1:numel(c.structures),a=r.structures(i);b=c.structures(i);if a.isTarget,keys={'D95','D98'};sense='lower';elseif strcmpi(a.type,'OAR')&&~ismember(upper(a.name),{'BODY','EXTERNAL'}),keys={'Dmean','Dmax'};sense='upper';else,continue;end;for k=1:numel(keys),m=keys{k};if strcmp(sense,'lower'),lim=.98*a.(m);ok=b.(m)>=lim;adv=max(0,100*(a.(m)-b.(m))/max(abs(a.(m)),eps));else,lim=a.(m)+max(.2,.02*max(abs(a.(m)),1));ok=b.(m)<=lim;adv=max(0,100*(b.(m)-a.(m))/max(abs(a.(m)),eps));end;rows(end+1,:)={a.name,m,b.(m),lim,ok};pass=pass&&ok;worst=max(worst,adv);end,end,end
function [g,c]=objective_information(model,A,w),[~,g]=evaluate_matrad_dose_model(model,w,true);x=A*w;c=repmat(struct('value',NaN,'gradient',zeros(numel(w),1)),0,1);for i=1:size(model.cst,1),vox=reshape(double(model.cst{i,4}{1}),[],1);for j=1:numel(model.cst{i,6}),o=model.cst{i,6}{j};if ~isa(o,'DoseObjectives.matRad_DoseObjective'),continue;end;o=model.optiProb.BP.setBiologicalDosePrescriptions(o,model.cst{i,5}.alphaX,model.cst{i,5}.betaX);dg=zeros(size(x));dg(vox)=o.penalty*o.computeDoseObjectiveGradient(x(vox));c(end+1)=struct('value',double(o.penalty*o.computeDoseObjectiveFunction(x(vox))),'gradient',A'*dg);end,end,end
function r=attempt_row(level,rep,batch,kind,inc,raw,plan,deaths,failure,gate,solver,wall),r=empty_attempt();r.level=level;r.replicate=rep;r.kind=kind;r.deletionSize=batch;r.proposedDeaths=numel(deaths);r.survivorFailure=failure;r.objectiveBefore=inc.objective;r.objectiveRaw=raw.objective;r.objectiveAfter=plan.objective;r.activeSpotsAfter=plan.activeSpots;r.activeLayersAfter=plan.activeLayers;r.bdtAfter=plan.delivery.T_total;r.dosePass=gate.dosePass;r.objectivePass=gate.objectivePass;r.planAccepted=gate.pass;r.solverStatus=solver.statusMessage;r.primalInfeasibility=readfield(solver,'primalInfeasibility',readnested(solver,'postKKT','primalInfeasibility'));r.dualInfeasibility=readfield(solver,'dualInfeasibility',readnested(solver,'postKKT','dualInfeasibility'));r.numericallyStationary=solver.numericallyStationary;r.wallSeconds=wall;end
function r=rescue_row(level,batch,inc,raw,plan,deaths,permanent,returns,births,gate,solver,wall),r=attempt_row(level,3,batch,'DORMANT_RESCUE',inc,raw,plan,deaths,~gate.pass,gate,solver,wall);r.permanentDeaths=numel(permanent);r.sameEventReturns=numel(returns);r.historicalDormantBirths=numel(births);r.netSupportChange=plan.activeSpots-inc.activeSpots;end
function r=empty_attempt(),r=struct('level',0,'replicate',0,'kind','','deletionSize',0,'proposedDeaths',0,'permanentDeaths',0,'sameEventReturns',0,'historicalDormantBirths',0,'netSupportChange',0,'survivorFailure',false,'objectiveBefore',NaN,'objectiveRaw',NaN,'objectiveAfter',NaN,'activeSpotsAfter',0,'activeLayersAfter',0,'bdtAfter',NaN,'dosePass',false,'objectivePass',false,'planAccepted',false,'solverStatus','','primalInfeasibility',NaN,'dualInfeasibility',NaN,'numericallyStationary',false,'wallSeconds',NaN);end
function r=empty_gate(),r=struct('level',0,'replicate',0,'phase','','structure','','metric','','value',NaN,'limit',NaN,'pass',false);end
function r=gate_rows(level,rep,phase,p,ref),[~,~,d]=dose_quality(p.dose,ref.dose);r=repmat(empty_gate(),size(d,1),1);for i=1:size(d,1),r(i)=struct('level',level,'replicate',rep,'phase',phase,'structure',d{i,1},'metric',d{i,2},'value',d{i,3},'limit',d{i,4},'pass',d{i,5});end,end
function r=empty_runtime(),r=struct('level',0,'replicate',0,'phase','','solverWallSeconds',NaN,'iterations',NaN,'objectiveCalls',0,'gradientCalls',0,'constraintCalls',0,'jacobianCalls',0,'finiteRemovalObjectiveCalls',0);end
function r=runtime_row(k,q,p,s,c),r=empty_runtime();r.level=k;r.replicate=q;r.phase=p;r.solverWallSeconds=s.wallSeconds;r.iterations=s.iterations;r.objectiveCalls=s.callbackCounts.objective;r.gradientCalls=s.callbackCounts.gradient;r.constraintCalls=readfield(s.callbackCounts,'constraints',0);r.jacobianCalls=readfield(s.callbackCounts,'jacobian',0);r.finiteRemovalObjectiveCalls=c;end
function r=empty_mapping(),r=struct('level',0,'deathGlobalIx','','birthGlobalIx',0,'birthLayer',0,'birthControlPoint',0,'weightBefore',0,'weightAfter',0);end
function r=mapping_rows(k,d,b,map,groups,w0,w1),r=repmat(empty_mapping(),numel(b),1);go=zeros(numel(map),1);for i=1:numel(groups),go(groups(i).ix)=i;end;ds=index_list([map(d).originalGlobalIx]);for i=1:numel(b),q=b(i);r(i)=struct('level',k,'deathGlobalIx',ds,'birthGlobalIx',map(q).originalGlobalIx,'birthLayer',go(q),'birthControlPoint',map(q).cpIndex,'weightBefore',w0(q),'weightAfter',w1(q));end,end
function flush(o,a,b,c,d),if ~isempty(a),writetable(struct2table(a,'AsArray',true),fullfile(o,'boundary_attempts.csv'));end;if ~isempty(b),writetable(struct2table(b,'AsArray',true),fullfile(o,'clinical_gate_trajectory.csv'));end;if ~isempty(c),writetable(struct2table(c,'AsArray',true),fullfile(o,'runtime_evaluations.csv'));end;if ~isempty(d),writetable(struct2table(d,'AsArray',true),fullfile(o,'death_birth_mapping.csv'));end,end
function write_report(o,x),f=fopen(fullfile(o,'SURVIVOR_COMPRESSION_BOUNDARY_REPORT.md'),'w');c=onCleanup(@()fclose(f));fprintf(f,'# Survivor-compression boundary\n\n- `%s`\n- repeated dormant replacement: `%s`\n- numerical stationarity: `%s`\n- boundary reproduced by matched block: %d\n- final active spots/layers: %d/%d\n- final objective: %.9g\n- final BDT: %.3f s\n- final clinical gates pass: %d\n',x.boundaryLabel,x.dormantLabel,x.stationarityLabel,x.boundaryReproduced,x.final.activeSpots,x.final.activeLayers,x.final.objective,x.final.delivery.T_total,x.finalGate.pass);end
function x=rank01(v),x=zeros(size(v));q=isfinite(v);if any(q),a=min(v(q));b=max(v(q));x(q)=(v(q)-a)/max(b-a,eps);end,end
function ix=top(score,mask,k),q=find(mask&isfinite(score));[~,o]=sort(score(q),'descend');ix=q(o(1:min(k,numel(o))));end
function D=nominal(d),D=d.physicalDose;if iscell(D),D=D{1};end,end
function A=state_matrix(m),A=nominal(m.dij);if isa(m.optiProb.BP,'matRad_ConstantRBEProjection'),A=double(m.dij.RBE)*A;end,end
function g=build_groups(map,toFull),m=map(toFull);for i=1:numel(m),m(i).globalIx=i;end;g=build_group_index(m);end
function n=count_layers(w,g,t),n=0;for i=1:numel(g),n=n+any(w(g(i).ix)>t);end,end
function s=index_list(ix),s=strjoin(arrayfun(@num2str,reshape(ix,1,[]),'UniformOutput',false),';');end
function v=readfield(s,n,d),if isstruct(s)&&isfield(s,n),v=s.(n);else,v=d;end,end
function v=readnested(s,a,b),v=NaN;if isstruct(s)&&isfield(s,a)&&isfield(s.(a),b),v=s.(a).(b);end,end
function add_paths(root,here),f={'config','planning','mapping','objective','delivery','reduction','evaluation','integrated','robustness'};for i=1:numel(f),addpath(fullfile(root,f{i}));end;addpath(fullfile(root,'experiments','repeated_dormant_replacement'));addpath(here);end
function create_empty(d),if exist(d,'dir')==7,e=dir(d);e=e(~ismember({e.name},{'.','..'}));assert(isempty(e));else,mkdir(d);end,end
function m=merge_structs(a,b),m=a;f=fieldnames(b);for i=1:numel(f),k=f{i};if isfield(m,k)&&isstruct(m.(k))&&isstruct(b.(k)),m.(k)=merge_structs(m.(k),b.(k));else,m.(k)=b.(k);end,end,end
