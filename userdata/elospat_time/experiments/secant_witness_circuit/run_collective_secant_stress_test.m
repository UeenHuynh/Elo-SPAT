function result=run_collective_secant_stress_test(outputDir)
%RUN_COLLECTIVE_SECANT_STRESS_TEST Final bounded Candidate-1 activation test.
narginchk(1,1); if exist(outputDir,'dir')~=7,mkdir(outputDir);end
here=fileparts(mfilename('fullpath'));root=fileparts(fileparts(here));add_paths(root);
artifactFile=fullfile(root,'results','task25_27_pipeline_20deg.mat');
stage1File=fullfile(root,'results','full_integrated_layer_20deg_job22180','alpha_0p3','alpha_0p3_final.mat');
source=load(artifactFile,'cfg','machine','dij','stf','pln','map');
stage1=load(stage1File,'final','finalDose','finalDelivery','supportMask');
patient=load(source.cfg.patientFile,'ct','cst'); cfg=merge_structs(default_config(),source.cfg);
rankTol=1e-10;zeroTol=1e-10;maxAttemptsPerSet=3;deathSizes=[2;4;8];
rule=default_stabilization_rule();
w0=reshape(double(stage1.final.w),[],1);allowed=reshape(logical(stage1.supportMask),[],1);
assert(nnz(allowed)==136409);[groups,~]=build_group_index(source.map);
allowedLayer=arrayfun(@(g)any(allowed(groups(g).ix)),(1:numel(groups))');
assert(nnz(allowedLayer)==748&&nnz(~allowedLayer)==20);
excluded=false(size(w0));for g=find(~allowedLayer)',excluded(groups(g).ix)=true;end
model=build_matrad_dose_model(source.dij,patient.cst,source.pln,0);
[A,D,rbe,stateName]=state_matrix(model);x0=A*w0;
baseline=evaluate_state(w0,A,D,model,patient,cfg,source.map,groups,source.machine);
assert(abs(baseline.objective-3139.658625654965)<=1e-8*max(1,baseline.objective));
assert(baseline.maxConstraintViolation<=rule.constraintTolerance);

[criticalVox,criticalInfo]=target_critical_voxels(D*w0,patient.cst);
[deathGlobal,deathScores]=select_prescribed_deaths(D,w0,allowed,criticalVox,source.map,groups,8);
[patchIx,patchInfo]=collective_patch(deathGlobal,allowed,source.map,source.stf,groups,240);
assert(all(ismember(deathGlobal,patchIx))&&all(allowed(patchIx))&&~any(excluded(patchIx)));
Apatch=A(:,patchIx); [Q0,witnessLabels]=base_witnesses(model,x0,patient.cst,rankTol);
assert(size(Q0,2)<numel(patchIx));
deathSets=cell(numel(deathSizes),1);
for k=1:numel(deathSizes),deathSets{k}=deathGlobal(1:deathSizes(k));end

names={'fixed_witness','gradient_only','adaptive_secant'};
rows=repmat(empty_row(),0,1);compRows=repmat(empty_comp(),0,1);variants=struct();
for i=1:numel(names)
    [vr,rr,cc]=run_variant(names{i},w0,patchIx,deathSets,A,Apatch,D,rbe,model, ...
        patient,cfg,source.map,groups,source.machine,allowed,excluded,baseline,Q0, ...
        rankTol,zeroTol,maxAttemptsPerSet,rule.constraintTolerance);
    variants.(names{i})=vr;rows=[rows;rr];compRows=[compRows;cc]; %#ok<AGROW>
end
T=struct2table(rows,'AsArray',true);C=struct2table(compRows,'AsArray',true);
writetable(T,fullfile(outputDir,'collective_secant_transitions.csv'));
writetable(C,fullfile(outputDir,'collective_secant_component_changes.csv'));
decision=decide(variants,T);
result=struct('decision',decision,'baseline',baseline,'patchIx',patchIx, ...
    'patchInfo',patchInfo,'criticalInfo',criticalInfo,'criticalVox',criticalVox, ...
    'deathGlobal',deathGlobal,'deathScores',deathScores,'deathSets',{deathSets}, ...
    'baseWitnessRank',size(Q0,2),'witnessLabels',{witnessLabels},'variants',variants, ...
    'stage1File',stage1File,'artifactFile',artifactFile,'stateName',stateName);
save(fullfile(outputDir,'collective_secant_result.mat'),'result','-v7.3');
write_report(outputDir,result,T,C);fprintf('FINAL_DECISION: %s\n',decision.label);
end

function [vr,rows,compRows]=run_variant(name,w0,patchIx,deathSets,A,Apatch,D,rbe,model, ...
    patient,cfg,map,groups,machine,allowed,excluded,baseline,Qfixed,rankTol,zeroTol,maxAttempts,constraintTolerance)
rows=repmat(empty_row(),0,1);compRows=repmat(empty_comp(),0,1);
oracle=0;acceptedSets=0;rejected=0;verified=0;geometryChanges=0;acceptedZeros=0;
for s=1:numel(deathSets)
    deathGlobal=deathSets{s};[tf,deathLocal]=ismember(deathGlobal,patchIx);assert(all(tf));
    secants=zeros(size(A,1),0);previousRejectedH=[];scenarioDone=false;
    for attempt=1:maxAttempts
        if strcmp(name,'fixed_witness'),Q=Qfixed;
        else,Q=base_witnesses(model,A*w0,patient.cst,rankTol);end
        if strcmp(name,'adaptive_secant')
            for q=1:size(secants,2),[Q,~]=append_vector(Q,secants(:,q),rankTol);end
        end
        Cker=Q'*Apatch;rankBefore=numrank(Cker,rankTol);
        [h,qp]=solve_collective_exchange_qp(Cker,w0(patchIx),deathLocal,rankTol);
        row=empty_row();row.variant=name;row.deathSetSize=numel(deathGlobal);row.attempt=attempt;
        row.prescribedDeathGlobalIx=index_list(deathGlobal);row.patchSize=numel(patchIx);
        row.rankBefore=rankBefore;row.kernelDimensionBefore=numel(patchIx)-rankBefore;
        row.svdSeconds=qp.svdSeconds;row.qpSeconds=qp.qpSeconds;row.qpExitflag=qp.exitflag;
        if ~qp.success,row.status='qp_failed';rows(end+1,1)=row;break;end %#ok<AGROW>
        assert(qp.kernelResidual<=1e-7&&qp.deathResidual<=1e-8);
        if ~isempty(previousRejectedH)
            row.exchangeDirectionRelativeChange=norm(h-previousRejectedH)/max(1,norm(previousRejectedH));
            row.geometryChanged=row.exchangeDirectionRelativeChange>1e-6;geometryChanges=geometryChanges+row.geometryChanged;
        end
        trial=w0;trial(patchIx)=w0(patchIx)+h;trial(abs(trial)<zeroTol)=0;trial=max(trial,0);
        assert(min(trial)>=0&&all(trial(deathGlobal)<=zeroTol));
        assert(~any(trial(excluded)>zeroTol));
        killed=find(w0>zeroTol&trial<=zeroTol);born=find(w0<=zeroTol&trial>zeroTol);
        row.actualDeaths=numel(killed);row.actualBirths=numel(born);row.netExactZeros=numel(killed)-numel(born);
        row.killedGlobalIx=index_list(killed);row.bornGlobalIx=index_list(born);
        row.activeLayers=count_layers(trial,groups,zeroTol);
        t=tic;proposal=evaluate_state(trial,A,D,model,patient,cfg,map,groups,machine);
        row.clinicalEvaluationSeconds=toc(t);oracle=oracle+1;row.oracleCall=oracle;
        [dosePass,worst,failure]=dose_gate_detail(proposal.dose,baseline.dose, ...
            cfg.evaluation.ctvD95TolerancePercent,patient.cst,D*w0,rbe);
        supportPass=~any(trial(~allowed)>zeroTol);
        constraintPass=proposal.maxConstraintViolation<=constraintTolerance;
        accept=dosePass&&supportPass&&constraintPass;
        row.accepted=accept;row.status=choose(accept,'accepted','rejected');row.worstAdversePercent=worst;
        row.objectiveBefore=baseline.objective;row.objectiveAfter=proposal.objective;
        row.objectiveChange=proposal.objective-baseline.objective;row.maxConstraintViolation=proposal.maxConstraintViolation;
        row.provisionalBDT=proposal.delivery.T_total;row.supportPreserved=supportPass;
        row.failureComponent=failure.label;row.failureComponentChange=failure.deltaH;
        compRows=[compRows;component_rows(name,numel(deathGlobal),attempt,baseline.components,proposal.components)]; %#ok<AGROW>
        if accept
            acceptedSets=acceptedSets+1;acceptedZeros=acceptedZeros+row.netExactZeros;scenarioDone=true;
            row.rankAfter=rankBefore;row.kernelDimensionAfter=numel(patchIx)-rankBefore;
        else
            rejected=rejected+1;
            if strcmp(name,'adaptive_secant')
                v=Apatch*h;
                [Qnew,sec]=append_secant_witness(Q,failure.a0,v,failure.deltaH,rankTol);
                row.secantIdentityError=sec.secantError;row.secantIndependent=sec.independent;
                if sec.independent
                    oldResidual=norm(Q'*v)/max(1,norm(v));newResidual=norm(Qnew'*v)/max(1,norm(v));
                    row.oldDirectionOldKernelResidual=oldResidual;row.oldDirectionNewKernelResidual=newResidual;
                    row.oldDirectionExcluded=newResidual>max(1e-8,100*oldResidual);
                    assert(row.oldDirectionExcluded);secants(:,end+1)=sec.newWitness;verified=verified+1; %#ok<AGROW>
                end
                row.rankAfter=numrank(Qnew'*Apatch,rankTol);row.kernelDimensionAfter=numel(patchIx)-row.rankAfter;
                previousRejectedH=h;
                if ~sec.independent,scenarioDone=true;end
            else
                row.rankAfter=rankBefore;row.kernelDimensionAfter=numel(patchIx)-rankBefore;scenarioDone=true;
            end
        end
        rows(end+1,1)=row; %#ok<AGROW>
        fprintf('%s deaths=%d attempt=%d accepted=%d actual=%d rank=%d->%d adverse=%.3f%%\n', ...
            name,numel(deathGlobal),attempt,accept,row.actualDeaths,row.rankBefore,row.rankAfter,worst);
        if scenarioDone,break;end
    end
end
vr=struct('oracleCalls',oracle,'acceptedSets',acceptedSets,'rejectedTransitions',rejected, ...
    'verifiedSecants',verified,'geometryChanges',geometryChanges,'acceptedZeros',acceptedZeros, ...
    'acceptedZerosPerOracle',acceptedZeros/max(1,oracle),'rows',rows);
end

function s=evaluate_state(w,A,D,model,patient,cfg,map,groups,machine)
[f,~,d]=evaluate_matrad_dose_model(model,w,false);x=A*w;components=clinical_components(model,x);
assert(abs(sum([components.value])-d.doseObjective)<=1e-7*max(1,abs(d.doseObjective)));
dose=compute_dvh_metrics(full(D*w),patient.ct,patient.cst,cfg);
delivery=compute_delivery_metrics(w,map,groups,cfg,machine);
s=struct('objective',f,'components',components,'dose',dose,'delivery',delivery, ...
    'maxConstraintViolation',d.maxConstraintViolation);
end

function [A,D,rbe,name]=state_matrix(model)
D=model.dij.physicalDose;if iscell(D),D=D{1};end;assert(issparse(D));
if isa(model.optiProb.BP,'matRad_DoseProjection')
    rbe=1;A=D;name='physicalDose';
elseif isa(model.optiProb.BP,'matRad_ConstantRBEProjection')
    rbe=double(model.dij.RBE);A=rbe*D;name='constant-RBE dose';
else
    error('Unsupported nonlinear state projection.');
end
end

function comps=clinical_components(model,x)
n=numel(x);comps=repmat(struct('label','','value',NaN,'gradient',zeros(n,1)),0,1);
for i=1:size(model.cst,1)
    vox=reshape(double(model.cst{i,4}{1}),[],1);if isempty(vox),continue;end
    for j=1:numel(model.cst{i,6})
        o=model.cst{i,6}{j};if ~isa(o,'DoseObjectives.matRad_DoseObjective'),continue;end
        assert(strcmp(o.robustness,'none'));
        o=model.optiProb.BP.setBiologicalDosePrescriptions(o,model.cst{i,5}.alphaX,model.cst{i,5}.betaX);
        c=struct('label',sprintf('%s:%s',char(string(model.cst{i,2})),class(o)), ...
            'value',double(o.penalty*o.computeDoseObjectiveFunction(x(vox))),'gradient',zeros(n,1));
        c.gradient(vox)=o.penalty*o.computeDoseObjectiveGradient(x(vox));comps(end+1,1)=c; %#ok<AGROW>
    end
end
end

function [Q,labels]=base_witnesses(model,x,cst,rankTol)
comps=clinical_components(model,x);Q=zeros(numel(x),0);labels={};
for i=1:numel(comps)
    [Q,a]=append_vector(Q,comps(i).gradient,rankTol);
    if a,labels{end+1}=comps(i).label;end %#ok<AGROW>
end
for s=1:size(model.cst,1)
    vox=reshape(double(model.cst{s,4}{1}),[],1);if isempty(vox),continue;end
    type=char(string(model.cst{s,3}));
    if strcmpi(type,'TARGET'),[~,o]=sort(x(vox),'ascend');
    elseif strcmpi(type,'OAR'),[~,o]=sort(x(vox),'descend');
    else,continue;end
    for k=1:min(2,numel(o))
        e=zeros(numel(x),1);e(vox(o(k)))=1;[Q,a]=append_vector(Q,e,rankTol);
        if a,labels{end+1}=sprintf('critical:%s:%d',char(string(cst{s,2})),vox(o(k)));end %#ok<AGROW>
    end
end
end

function [Q,a]=append_vector(Q,q,tol)
q=reshape(double(q),[],1);
if isempty(Q),r=q;else,r=q-Q*(Q'*q);r=r-Q*(Q'*r);end
a=norm(r)>tol*max(1,norm(q));
if a,Q=[Q,r/norm(r)];end %#ok<AGROW>
end

function [vox,info]=target_critical_voxels(dose,cst)
vox=zeros(0,1);names={};
for s=1:size(cst,1)
    if ~strcmpi(cst{s,3},'TARGET'),continue;end
    ix=reshape(double(cst{s,4}{1}),[],1);[~,o]=sort(dose(ix),'ascend');n=numel(ix);
    pos=unique(max(1,min(n,round(1+(n-1)*[0.015 0.02 0.03 0.04 0.05 0.06]))));
    vox=[vox;ix(o(pos))];names{end+1}=char(string(cst{s,2})); %#ok<AGROW>
end
vox=unique(vox,'stable');info=struct('nVoxels',numel(vox),'targetNames',{names});
end

function [chosen,scoresChosen]=select_prescribed_deaths(D,w,allowed,criticalVox,map,groups,k)
score=w.*reshape(full(sum(D(criticalVox,:),1)),[],1)/max(1,numel(criticalVox));score(~allowed)=-Inf;
[~,ord]=sort(score,'descend');chosen=zeros(0,1);cpCount=zeros(max([map.cpIndex]),1);
layerCount=zeros(numel(groups),1);groupOf=zeros(numel(w),1);
for g=1:numel(groups),groupOf(groups(g).ix)=g;end
for q=1:numel(ord)
    i=ord(q);cp=map(i).cpIndex;g=groupOf(i);
    if cpCount(cp)>=2||layerCount(g)>=2,continue;end
    chosen(end+1,1)=i;cpCount(cp)=cpCount(cp)+1;layerCount(g)=layerCount(g)+1; %#ok<AGROW>
    if numel(chosen)==k,break;end
end
assert(numel(chosen)==k);scoresChosen=score(chosen);
end

function [patch,info]=collective_patch(deaths,allowed,map,stf,groups,maxSize)
groupOf=zeros(numel(map),1);for g=1:numel(groups),groupOf(groups(g).ix)=g;end
patch=deaths(:);neighborhoods=repmat(struct('death',0,'cp',0,'energy',0,'bev',[0 0]),numel(deaths),1);
for q=1:numel(deaths)
    i=deaths(q);g0=groupOf(i);cp=map(i).cpIndex;e0=groups(g0).energyIndex;center=xy(stf,map(i));
    candidateGroups=find([groups.cpIndex]'==cp&abs([groups.energyIndex]'-e0)<=2);
    cand=vertcat(groups(candidateGroups).ix);cand=cand(allowed(cand));pos=zeros(numel(cand),2);
    for j=1:numel(cand),pos(j,:)=xy(stf,map(cand(j)));end
    [~,o]=sort(sum((pos-center).^2,2));patch=[patch;cand(o(1:min(28,numel(o))))]; %#ok<AGROW>
    neighborhoods(q)=struct('death',i,'cp',cp,'energy',map(i).energy,'bev',center);
end
patch=unique(patch,'stable');
if numel(patch)>maxSize
    keep=ismember(patch,deaths);patch=[patch(keep);patch(~keep)];patch=patch(1:maxSize);
end
info=struct('size',numel(patch),'nNeighborhoods',numel(deaths),'neighborhoods',neighborhoods, ...
    'controlPoints',unique([map(patch).cpIndex]),'energies',unique([map(patch).energy]));
end
function p=xy(stf,m),r=stf(m.beamNum).ray(m.rayNum);q=double(r.rayPos_bev);p=[q(1),q(3)];end

function [pass,worst,f]=dose_gate_detail(candidate,reference,tol,cst,physicalDose,rbe)
worst=0;f=struct('label','','deltaH',0,'a0',zeros(numel(physicalDose),1));
for s=1:numel(reference.structures)
    a=reference.structures(s);b=candidate.structures(s);tests={};
    if a.isTarget,tests={{'D95',-1},{'D98',-1}};
    elseif strcmpi(a.type,'OAR'),tests={{'Dmean',1},{'Dmax',1}};end
    for j=1:numel(tests)
        metric=tests{j}{1};sgn=tests{j}{2};adverse=max(0,sgn*pct(b.(metric),a.(metric)));
        if adverse>worst
            worst=adverse;f.label=sprintf('%s:%s',a.name,metric);
            f.deltaH=sgn*(b.(metric)-a.(metric));
            f.a0=sgn*metric_gradient(physicalDose,cst{s,4}{1},metric,rbe);
        end
    end
end
pass=worst<=tol;
end

function g=metric_gradient(dose,vox,metric,rbe)
g=zeros(numel(dose),1);vox=reshape(double(vox),[],1);d=dose(vox);n=numel(vox);
if strcmp(metric,'Dmean')
    g(vox)=1/(n*rbe);
elseif strcmp(metric,'Dmax')
    [~,i]=max(d);g(vox(i))=1/rbe;
else
    if strcmp(metric,'D95'),p=5;else,p=2;end
    [~,o]=sort(d,'ascend');pos=1+(n-1)*p/100;lo=floor(pos);hi=ceil(pos);a=pos-lo;
    g(vox(o(lo)))=g(vox(o(lo)))+(1-a)/rbe;
    g(vox(o(hi)))=g(vox(o(hi)))+a/rbe;
end
end

function rows=component_rows(v,k,a,before,after)
rows=repmat(empty_comp(),numel(before),1);
for i=1:numel(before)
    rows(i).variant=v;rows(i).deathSetSize=k;rows(i).attempt=a;rows(i).component=before(i).label;
    rows(i).before=before(i).value;rows(i).after=after(i).value;rows(i).change=after(i).value-before(i).value;
end
end

function d=decide(v,T)
a=v.adaptive_secant;g=v.gradient_only;anyReject=any(~T.accepted&T.oracleCall>0);
equiv=a.acceptedZeros==g.acceptedZeros&&a.acceptedSets==g.acceptedSets;
rankSat=any(T.rankAfter>=0.9*T.patchSize);lowYield=a.acceptedZerosPerOracle<=1.25;
if equiv||lowYield||rankSat
    label='MECHANISM_NOT_SUPPORTED';
elseif ~anyReject
    label='INCONCLUSIVE';
elseif a.verifiedSecants==0||a.geometryChanges==0
    label='MECHANISM_NOT_SUPPORTED';
else
    label='MECHANISM_SUPPORTED';
end
d=struct('label',label,'anyClinicalRejection',anyReject,'verifiedSecants',a.verifiedSecants, ...
    'geometryChanges',a.geometryChanges,'adaptiveGradientEquivalent',equiv, ...
    'lowYield',lowYield,'rankSaturation',rankSat);
end

function write_report(out,r,T,C)
fid=fopen(fullfile(out,'COLLECTIVE_SECANT_STRESS_REPORT.md'),'w');assert(fid>0);cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'# Candidate 1 collective-boundary stress test\n\nDecision: **%s**\n\n',r.decision.label);
fprintf(fid,'Frozen Stage-1: objective %.12g, constraint violation %.3e, BDT %.6f s. No production path or tolerance changed.\n\n', ...
    r.baseline.objective,r.baseline.maxConstraintViolation,r.baseline.delivery.T_total);
fprintf(fid,'Patch: %d columns across %d CPs and %d energies; base witness rank %d. Prescribed death IDs: `%s`.\n\n', ...
    r.patchInfo.size,numel(r.patchInfo.controlPoints),numel(r.patchInfo.energies), ...
    r.baseWitnessRank,index_list(r.deathGlobal));
fprintf(fid,'| Variant | oracle calls | accepted sets | rejected transitions | verified secants | changed geometries | accepted zeros | zeros/oracle |\n');
fprintf(fid,'|---|---:|---:|---:|---:|---:|---:|---:|\n');
names={'fixed_witness','gradient_only','adaptive_secant'};
for i=1:numel(names)
    z=r.variants.(names{i});fprintf(fid,'| %s | %d | %d | %d | %d | %d | %d | %.3f |\n', ...
        names{i},z.oracleCalls,z.acceptedSets,z.rejectedTransitions,z.verifiedSecants, ...
        z.geometryChanges,z.acceptedZeros,z.acceptedZerosPerOracle);
end
fprintf(fid,'\nTransitions: %d; component rows: %d. Rejection observed=%d; adaptive/gradient equivalent=%d; low-yield stop=%d; rank saturation=%d.\n\n', ...
    height(T),height(C),r.decision.anyClinicalRejection,r.decision.adaptiveGradientEquivalent, ...
    r.decision.lowYield,r.decision.rankSaturation);
fprintf(fid,'All endpoints used the unchanged 2%% dose-quality gate, exact matRad objective/constraints, support exclusion, and discrete BDT evaluation. No IPOPT repair was run.\n');
end

function r=empty_row()
r=struct('variant','','deathSetSize',0,'attempt',0,'prescribedDeathGlobalIx','', ...
    'patchSize',0,'rankBefore',0,'rankAfter',0,'kernelDimensionBefore',0, ...
    'kernelDimensionAfter',0,'qpExitflag',NaN,'actualDeaths',0,'actualBirths',0, ...
    'netExactZeros',0,'killedGlobalIx','','bornGlobalIx','','activeLayers',0, ...
    'objectiveBefore',NaN,'objectiveAfter',NaN,'objectiveChange',NaN, ...
    'maxConstraintViolation',NaN,'worstAdversePercent',NaN,'failureComponent','', ...
    'failureComponentChange',NaN,'supportPreserved',false,'accepted',false,'status','', ...
    'secantIdentityError',NaN,'secantIndependent',false, ...
    'oldDirectionOldKernelResidual',NaN,'oldDirectionNewKernelResidual',NaN, ...
    'oldDirectionExcluded',false,'exchangeDirectionRelativeChange',NaN, ...
    'geometryChanged',false,'provisionalBDT',NaN,'oracleCall',0,'svdSeconds',NaN, ...
    'qpSeconds',NaN,'clinicalEvaluationSeconds',NaN);
end
function r=empty_comp()
r=struct('variant','','deathSetSize',0,'attempt',0,'component','','before',NaN,'after',NaN,'change',NaN);
end
function n=count_layers(w,g,t),n=0;for i=1:numel(g),n=n+(norm(w(g(i).ix))>t);end,end
function r=numrank(C,t),if isempty(C),r=0;else,s=svd(C,'econ');if isempty(s),r=0;else,r=sum(s>t*max(1,s(1)));end,end,end
function s=index_list(ix),s=strjoin(arrayfun(@(x)sprintf('%d',x),reshape(ix,1,[]),'UniformOutput',false),';');end
function y=pct(x,r),y=100*(double(x)-double(r))/max(abs(double(r)),eps);end
function x=choose(c,a,b),if c,x=a;else,x=b;end,end
function add_paths(root),f={'config','planning','mapping','objective','delivery','reduction','robustness','evaluation','integrated'};for i=1:numel(f),addpath(fullfile(root,f{i}));end,end
function m=merge_structs(a,b),m=a;f=fieldnames(b);for i=1:numel(f),k=f{i};if isfield(m,k)&&isstruct(m.(k))&&isstruct(b.(k)),m.(k)=merge_structs(m.(k),b.(k));else,m.(k)=b.(k);end,end,end
