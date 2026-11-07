function evaluate_coreset_checkpoint(resultDir)
%EVALUATE_CORESET_CHECKPOINT Read-only validation of a saved coreset state.
here=fileparts(mfilename('fullpath'));root=fileparts(fileparts(here));add_paths(root);
artifact=load(fullfile(root,'results','task25_27_pipeline_20deg.mat'), ...
    'cfg','machine','dij','pln','map');
stage1=load(fullfile(root,'results','full_integrated_layer_20deg_job22180', ...
    'alpha_0p3','alpha_0p3_final.mat'),'final','finalDose');
files=dir(fullfile(resultDir,'coreset_iter_*.mat'));assert(~isempty(files));
iteration=zeros(numel(files),1);
for i=1:numel(files)
    token=regexp(files(i).name,'coreset_iter_(\d+)\.mat$','tokens','once');
    iteration(i)=str2double(token{1});
end
[~,lastIx]=max(iteration);
sourceCheckpoint=fullfile(files(lastIx).folder,files(lastIx).name);
saved=load(sourceCheckpoint,'stateValue');result=saved.stateValue;
patient=load(artifact.cfg.patientFile,'ct','cst');cfg=merge_structs(default_config(),artifact.cfg);
wFullRef=reshape(double(stage1.final.w),[],1);allowed=wFullRef>0;
[dijAllowed,allowedToFull]=reduce_dij_to_support(artifact.dij,allowed);
D=nominal_matrix(dijAllowed);dose=compute_dvh_metrics(full(D*result.w),patient.ct,patient.cst,cfg);
model=build_matrad_dose_model(dijAllowed,patient.cst,artifact.pln,0);
[objective,~,diagnostics]=evaluate_matrad_dose_model(model,result.w,false);
[groups,~]=build_group_index(artifact.map);z=zeros(size(wFullRef));z(allowedToFull)=result.w;
delivery=compute_delivery_metrics(z,artifact.map,groups,cfg,artifact.machine);
n=numel(dose.structures);rows=repmat(struct('structure','','type','', ...
    'referenceD98',NaN,'coresetD98',NaN,'referenceD95',NaN,'coresetD95',NaN, ...
    'referenceDmean',NaN,'coresetDmean',NaN,'referenceDmax',NaN,'coresetDmax',NaN),n,1);
for i=1:n
    a=stage1.finalDose.structures(i);b=dose.structures(i);rows(i).structure=a.name;rows(i).type=a.type;
    for q={'D98','D95','Dmean','Dmax'},x=q{1};rows(i).(['reference' x])=a.(x);rows(i).(['coreset' x])=b.(x);end
end
writetable(struct2table(rows),fullfile(resultDir,'coreset_final_structure_metrics.csv'));
last=result.history(end);summary=struct('activeSpots',nnz(result.w>1e-10), ...
    'activeLayers',last.activeLayers,'weightedRelativeError',last.weightedRelativeError, ...
    'clinicalObjective',objective,'maxConstraintViolation',diagnostics.maxConstraintViolation, ...
    'targetD95',dose.summary.targetD95Min,'targetD98',dose.summary.targetD98Min, ...
    'oarDmax',dose.summary.oarDmaxMax,'provisionalBDT',delivery.T_total, ...
    'coresetWallSeconds',result.wallSeconds,'sourceCheckpoint',sourceCheckpoint);
writetable(struct2table(summary,'AsArray',true), ...
    fullfile(resultDir,'coreset_final_validation.csv'));
save(fullfile(resultDir,'coreset_final_validation.mat'),'dose','objective', ...
    'diagnostics','delivery','summary','-v7.3');
end
function D=nominal_matrix(dij),if iscell(dij.physicalDose),D=dij.physicalDose{1};else,D=dij.physicalDose;end,end
function add_paths(root),f={'config','planning','mapping','objective','delivery','reduction','evaluation','integrated','robustness'};for i=1:numel(f),addpath(fullfile(root,f{i}));end,end
function m=merge_structs(a,b),m=a;f=fieldnames(b);for i=1:numel(f),k=f{i};if isfield(m,k)&&isstruct(m.(k))&&isstruct(b.(k)),m.(k)=merge_structs(m.(k),b.(k));else,m.(k)=b.(k);end,end,end
