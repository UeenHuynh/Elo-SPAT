function extract_bidirectional_results(resultDir)
%EXTRACT_BIDIRECTIONAL_RESULTS Read-only clinical report for frozen job 22312.

% This function never calls an optimizer.  It recomputes nominal dose metrics
% from the final frozen weights and writes a per-structure comparison against
% the accepted Stage-1 alphaEL=0.3 plan.

here = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(here));
addpath(fullfile(root,'evaluation'));
addpath(fullfile(root,'mapping'));
addpath(fullfile(root,'delivery'));
addpath(fullfile(root,'config'));

completeFile = fullfile(resultDir,'bidirectional_support_complete.mat');
stage1File = fullfile(root,'results','full_integrated_layer_20deg_job22180', ...
    'alpha_0p3','alpha_0p3_final.mat');
artifactFile = fullfile(root,'results','task25_27_pipeline_20deg.mat');

s = load(completeFile,'output');
stage1 = load(stage1File,'finalDose','supportMask');
source = load(artifactFile,'cfg','dij','map');
patient = load(source.cfg.patientFile,'ct','cst');
cfg = merge_structs(default_config(),source.cfg);

assert(~isempty(s.output.final),'Frozen job does not contain final polish.');
w = reshape(double(s.output.final.w),[],1);
allowed = reshape(logical(stage1.supportMask),[],1);
assert(numel(w)==numel(allowed));
excludedMaxWeight = max(abs(w(~allowed)));
assert(excludedMaxWeight <= 1e-10*max(1,max(abs(w))), ...
    'A Stage-1-excluded spot was reactivated.');

D = source.dij.physicalDose;
if iscell(D),D=D{1};end
finalDose = compute_dvh_metrics(full(D*w),patient.ct,patient.cst,cfg);
rows = dose_rows(stage1.finalDose,finalDose);
writetable(struct2table(rows),fullfile(resultDir, ...
    'final_ipopt_dose_vs_stage1.csv'));

[groups,~] = build_group_index(source.map);
activeLayers = 0;
for g = 1:numel(groups)
    activeLayers = activeLayers + any(w(groups(g).ix)>1e-10);
end
summary = table(s.output.reference.objective,s.output.finalSummary.objective, ...
    s.output.reference.delivery.T_total,s.output.finalSummary.bdt, ...
    nnz(allowed),nnz(w>1e-10),748,activeLayers,excludedMaxWeight, ...
    s.output.finalSummary.supportPreserved,s.output.finalSummary.solverAccepted, ...
    s.output.finalSummary.planAccepted,s.output.finalSummary.polishWallSeconds, ...
    'VariableNames',{'stage1Objective','finalObjective','stage1BDT','finalBDT', ...
    'stage1ActiveSpots','finalActiveSpots','stage1ActiveLayers', ...
    'finalActiveLayers','maxExcludedWeight','supportPreserved', ...
    'solverAccepted','planAccepted','polishWallSeconds'});
writetable(summary,fullfile(resultDir,'final_ipopt_summary.csv'));
save(fullfile(resultDir,'final_clinical_metrics.mat'),'finalDose','rows','summary');
fprintf('BIDIRECTIONAL_RESULT_EXTRACTION_COMPLETE\n');
disp(summary);
end

function rows = dose_rows(reference,current)
n = numel(reference.structures);
assert(numel(current.structures)==n);
rows = repmat(struct('structure','','type','','stage1D98',NaN,'finalD98',NaN, ...
    'stage1D95',NaN,'finalD95',NaN,'stage1Dmean',NaN,'finalDmean',NaN, ...
    'stage1Dmax',NaN,'finalDmax',NaN,'adverseChangePercent',NaN, ...
    'structurePass',false),n,1);
for i = 1:n
    a = reference.structures(i); b = current.structures(i);
    assert(strcmp(a.name,b.name));
    rows(i).structure = a.name; rows(i).type = a.type;
    rows(i).stage1D98 = a.D98; rows(i).finalD98 = b.D98;
    rows(i).stage1D95 = a.D95; rows(i).finalD95 = b.D95;
    rows(i).stage1Dmean = a.Dmean; rows(i).finalDmean = b.Dmean;
    rows(i).stage1Dmax = a.Dmax; rows(i).finalDmax = b.Dmax;
    if a.isTarget
        adverse = max([0,-percent_change(b.D95,a.D95), ...
            -percent_change(b.D98,a.D98)]);
        pass = adverse <= 2;
    elseif strcmpi(a.type,'OAR') && ~strcmpi(a.name,'BODY') && ...
            ~strcmpi(a.name,'EXTERNAL')
        dMean = b.Dmean-a.Dmean; dMax = b.Dmax-a.Dmax;
        pass = dMean<=max(0.2,0.02*max(abs(a.Dmean),1)) && ...
            dMax<=max(0.2,0.02*max(abs(a.Dmax),1));
        adverse = max([0,percent_change(b.Dmean,a.Dmean), ...
            percent_change(b.Dmax,a.Dmax)]);
    else
        adverse = 0; pass = true;
    end
    rows(i).adverseChangePercent = adverse;
    rows(i).structurePass = pass;
end
end

function value = percent_change(a,b)
value = 100*(double(a)-double(b))/max(abs(double(b)),eps);
end

function merged = merge_structs(a,b)
merged = a;
fields = fieldnames(b);
for i = 1:numel(fields)
    key = fields{i};
    if isfield(merged,key) && isstruct(merged.(key)) && isstruct(b.(key))
        merged.(key) = merge_structs(merged.(key),b.(key));
    else
        merged.(key) = b.(key);
    end
end
end
