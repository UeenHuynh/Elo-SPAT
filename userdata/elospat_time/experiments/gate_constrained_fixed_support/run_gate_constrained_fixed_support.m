function output=run_gate_constrained_fixed_support(outputDir)
% Isolate whether job-22353 support can satisfy the frozen LFH Dmax gate.
narginchk(1,1);here=fileparts(mfilename('fullpath'));root=fileparts(fileparts(here));add_paths(root);create_empty(outputDir);
artifactFile=fullfile(root,'results','task25_27_pipeline_20deg.mat');stage1File=fullfile(root,'results','full_integrated_layer_20deg_job22180','alpha_0p3','alpha_0p3_final.mat');
jobDir=fullfile(root,'results','gate_aware_event1_job22353');eventDir=fullfile(root,'results','dormant_replacement_job22325');
source=load(artifactFile,'cfg','machine','dij','stf','pln','map');stage1=load(stage1File,'final','finalDose','supportMask');
j=load(fullfile(jobDir,'gate_aware_complete.mat'),'output');j=j.output;sv=load(fullfile(eventDir,'event_01_survivor.mat'),'checkpoint');sv=sv.checkpoint;
patient=load(source.cfg.patientFile,'ct','cst');cfg=merge_structs(default_config(),source.cfg);zeroTol=1e-10;
D=source.dij.physicalDose;if iscell(D),D=D{1};end;[groups,~]=build_group_index(source.map);
wRef=reshape(double(stage1.final.w),[],1);wOrd=reshape(double(j.exact.w),[],1);support=reshape(logical(j.exact.supportMask),[],1);
historical=reshape(logical(sv.historicalDormant),[],1);deaths=reshape(double(sv.deathIx),[],1);rescue=reshape(double(j.proposal.rescueIx),[],1);
assert(nnz(support)==135939&&numel(deaths)==677&&numel(rescue)==256&&all(support(rescue))&&~any(support(~stage1.supportMask)));

% The frozen evaluator uses nominal physical dose and a literal voxel max.
sIx=find(strcmp(patient.cst(:,2),'Lt femoral head'),1);assert(~isempty(sIx));vox=reshape(double(patient.cst{sIx,4}{1}),[],1);
refDose=compute_dvh_metrics(full(D*wRef),patient.ct,patient.cst,cfg);ordDose=compute_dvh_metrics(full(D*wOrd),patient.ct,patient.cst,cfg);
refLFH=metric(refDose,'Lt femoral head','Dmax');ordLFH=metric(ordDose,'Lt femoral head','Dmax');
assert(abs(refLFH-max(D(vox,:)*wRef))<1e-10&&abs(ordLFH-max(D(vox,:)*wOrd))<1e-10);
gateUpper=refLFH+max(0.2,0.02*max(abs(refLFH),1));

[dijR,r2f,f2r]=reduce_dij_to_support(source.dij,support);model=build_matrad_dose_model(dijR,patient.cst,source.pln,0);
Dr=model.dij.physicalDose;if iscell(Dr),Dr=Dr{1};end;gateA=Dr(vox,:);assert(size(gateA,2)==nnz(support));
wFailed=reshape(double(sv.survivor.w),[],1);wProposal=wFailed;J2=j.proposal.J2;wProposal(J2)=max(0,wFailed(J2)+j.proposal.delta);wProposal(~support)=0;
w0=wProposal(r2f);warmMax=max(gateA*w0);scale=min(1,gateUpper/warmMax);w0=w0*scale;assert(max(gateA*w0)<=gateUpper+1e-10);
% IPOPT reports aggregate CPU time across the eight allocated cores.  Use a
% budget consistent with the 14-hour Slurm wall-time instead of 12 CPU-hours.
settings=struct('ipopt',struct('max_cpu_time',300000));solver=solve_gate_constrained_matrad_ipopt(model,w0,gateA,gateUpper,settings);
w=zeros(double(source.dij.totalNumOfBixels),1);w(r2f)=solver.w;assert(~any(abs(w(~support))>1e-12));
reference=evaluate_plan(wRef,D,patient,cfg,source,groups);ordinary=evaluate_plan(wOrd,D,patient,cfg,source,groups);constrained=evaluate_plan(w,D,patient,cfg,source,groups);
gate=acceptance_gate(constrained,reference,solver.solverAccepted);positive=w>zeroTol;sameReturn=deaths(positive(deaths));births=find(historical&positive);permanent=setdiff(deaths,sameReturn);net=nnz(positive)-nnz(sv.supportInc);
supported=gate.pass&&~isempty(births)&&net<0;
if supported,status='DORMANT_REPLACEMENT_SUPPORTED';elseif solver.solverAccepted,status='INCONCLUSIVE';else,status='INCONCLUSIVE';end
write_comparison(outputDir,ordinary,constrained,j,solver,gate,births,sameReturn,permanent,net,support,groups);
write_clinical(outputDir,reference,ordinary,constrained);
semantics=struct('quantity','nominal physical dose','operation','literal maximum over cst voxel indices', ...
 'structure','Lt femoral head','voxelCount',numel(vox),'gateUpperGy',gateUpper, ...
 'referenceDmaxGy',refLFH,'ordinaryDmaxGy',ordLFH,'equivalenceErrorReference',abs(refLFH-max(D(vox,:)*wRef)), ...
 'equivalenceErrorOrdinary',abs(ordLFH-max(D(vox,:)*wOrd)));
output=struct('status',status,'semantics',semantics,'supportMask',support,'rescueIx',rescue, ...
 'solver',solver,'w',w,'reference',reference,'ordinary',ordinary,'constrained',constrained,'gate',gate, ...
 'historicalDormantBirths',births,'sameEventReturns',sameReturn,'permanentDeaths',permanent, ...
 'netSupportChange',net,'warmStartScale',scale,'reducedToFullIdx',r2f,'fullToReducedIdx',f2r);
save(fullfile(outputDir,'gate_constrained_complete.mat'),'output','-v7.3');write_report(outputDir,output);fprintf('FINAL_STATUS: %s\n',status);
end
function p=evaluate_plan(w,D,patient,cfg,source,groups)
model=build_matrad_dose_model(source.dij,patient.cst,source.pln,0);[f,~,d]=evaluate_matrad_dose_model(model,w,false);
p=struct('objective',f,'maxConstraintViolation',d.maxConstraintViolation,'dose',compute_dvh_metrics(full(D*w),patient.ct,patient.cst,cfg), ...
 'delivery',compute_delivery_metrics(w,source.map,groups,cfg,source.machine),'activeSpots',nnz(w>1e-10),'activeLayers',count_layers(w,groups,1e-10));
end
function g=acceptance_gate(p,r,solverOK)
[dosePass,worst,detail]=dose_quality(p.dose,r.dose,2,0.2);g=struct('pass',solverOK&&dosePass&&p.objective<=1.02*r.objective&&p.maxConstraintViolation<=1e-6, ...
 'dosePass',dosePass,'objectivePass',p.objective<=1.02*r.objective,'constraintPass',p.maxConstraintViolation<=1e-6,'worstAdverse',worst,'detail',{detail});
end
function [pass,worst,rows]=dose_quality(c,r,pctTol,absTol)
pass=true;worst=0;rows={};
for i=1:numel(c.structures)
 a=r.structures(i);b=c.structures(i);
 if a.isTarget
  keys={'D95','D98'};
  for k=1:numel(keys),m=keys{k};lim=(1-pctTol/100)*a.(m);ok=b.(m)>=lim;adv=max(0,100*(a.(m)-b.(m))/max(abs(a.(m)),eps));rows(end+1,:)={a.name,m,b.(m),lim,ok};pass=pass&&ok;worst=max(worst,adv);end %#ok<AGROW>
 elseif strcmpi(a.type,'OAR')&&~strcmpi(a.name,'BODY')&&~strcmpi(a.name,'EXTERNAL')
  keys={'Dmean','Dmax'};
  for k=1:numel(keys),m=keys{k};lim=a.(m)+max(absTol,pctTol/100*max(abs(a.(m)),1));ok=b.(m)<=lim;adv=max(0,100*(b.(m)-a.(m))/max(abs(a.(m)),eps));rows(end+1,:)={a.name,m,b.(m),lim,ok};pass=pass&&ok;worst=max(worst,adv);end %#ok<AGROW>
 end
end
end
function write_comparison(out,a,b,j,s,g,births,returns,permanent,net,support,groups)
method={'job22353_objective_only';'same_support_gate_constrained'};objective=[a.objective;b.objective];lfh=[metric(a.dose,'Lt femoral head','Dmax');metric(b.dose,'Lt femoral head','Dmax')];
activeSpots=[a.activeSpots;b.activeSpots];activeLayers=[a.activeLayers;b.activeLayers];histBirths=[numel(j.historicalDormantBirths);numel(births)];sameReturns=[numel(j.sameEventReturns);numel(returns)];permanentDeaths=[numel(j.permanentDeaths);numel(permanent)];netChange=[j.netSupportChange;net];bdt=[a.delivery.T_total;b.delivery.T_total];
iterations=[j.solverMetrics.iterations;s.iterations];wall=[j.exactSeconds;s.wallSeconds];solverStatus={j.solverMetrics.statusMessage;s.statusMessage};dosePass=[j.exactGate.dosePass;g.dosePass];planPass=[j.exactGate.pass;g.pass];
objectiveCalls=[180;s.callbackCounts.objective];gradientCalls=[180;s.callbackCounts.gradient];constraintCalls=[0;s.callbackCounts.constraints];jacobianCalls=[0;s.callbackCounts.jacobian];
T=table(method,objective,lfh,dosePass,planPass,activeSpots,activeLayers,histBirths,sameReturns,permanentDeaths,netChange,bdt,iterations,objectiveCalls,gradientCalls,constraintCalls,jacobianCalls,wall,solverStatus);writetable(T,fullfile(out,'fixed_support_comparison.csv'));
S=table(nnz(support),count_layers(double(support),groups,0.5),'VariableNames',{'proposedSupportSpots','proposedSupportLayers'});writetable(S,fullfile(out,'support_identity.csv'));
end
function write_clinical(out,r,a,b)
rows={};names={'job22353_objective_only','same_support_gate_constrained'};p={a,b};for q=1:2,for i=1:numel(r.dose.structures),x=r.dose.structures(i);y=p{q}.dose.structures(i);rows(end+1,:)={names{q},x.name,x.type,x.D95,y.D95,x.D98,y.D98,x.Dmean,y.Dmean,x.Dmax,y.Dmax};end,end %#ok<AGROW>
T=cell2table(rows,'VariableNames',{'method','structure','type','referenceD95','trialD95','referenceD98','trialD98','referenceDmean','trialDmean','referenceDmax','trialDmax'});writetable(T,fullfile(out,'fixed_support_clinical_metrics.csv'));
end
function write_report(out,o),fid=fopen(fullfile(out,'GATE_CONSTRAINED_FIXED_SUPPORT_REPORT.md'),'w');assert(fid>0);c=onCleanup(@()fclose(fid));fprintf(fid,'# Gate-constrained solve on frozen job-22353 support\n\nStatus: **%s**\n\n',o.status);fprintf(fid,'Dmax semantics: nominal physical-dose maximum over %d LFH voxels. Gate: %.12g Gy.\n\n',o.semantics.voxelCount,o.semantics.gateUpperGy);fprintf(fid,'Objective %.12g; LFH Dmax %.12g Gy; dose gate pass %d; historical births %d; net support change %d.\n',o.constrained.objective,metric(o.constrained.dose,'Lt femoral head','Dmax'),o.gate.dosePass,numel(o.historicalDormantBirths),o.netSupportChange);end
function v=metric(d,n,m),i=find(strcmp({d.structures.name},n),1);assert(~isempty(i));v=d.structures(i).(m);end
function n=count_layers(w,g,t),n=0;for i=1:numel(g),n=n+any(w(g(i).ix)>t);end,end
function add_paths(root),f={'config','planning','mapping','objective','delivery','reduction','evaluation','integrated','robustness'};for i=1:numel(f),addpath(fullfile(root,f{i}));end,end
function create_empty(d),if exist(d,'dir')==7,e=dir(d);e=e(~ismember({e.name},{'.','..'}));assert(isempty(e));else,mkdir(d);end,end
function m=merge_structs(a,b),m=a;f=fieldnames(b);for i=1:numel(f),k=f{i};if isfield(m,k)&&isstruct(m.(k))&&isstruct(b.(k)),m.(k)=merge_structs(m.(k),b.(k));else,m.(k)=b.(k);end,end,end
