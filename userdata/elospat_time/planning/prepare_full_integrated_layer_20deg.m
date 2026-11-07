function outputs = prepare_full_integrated_layer_20deg(outputDir, alphaGrid)
% Compute full-problem dimensionless lambda scaling without optimization.

narginchk(0,2);
root=fileparts(fileparts(mfilename('fullpath'))); add_module_paths(root);
if nargin<1 || isempty(outputDir), outputDir=fullfile(root,'results', ...
        ['full_integrated_layer_setup_20deg_' datestr(now,'yyyymmdd_HHMMSS')]); end
if nargin<2 || isempty(alphaGrid), alphaGrid=[0.01;0.03;0.1;0.3;1]; end
alphaGrid=reshape(double(alphaGrid),[],1); create_empty_output_directory(outputDir);
artifactFile=fullfile(root,'results','task25_27_pipeline_20deg.mat');
source=load(artifactFile,'cfg','machine','denseResult','dij','pln','map','groups');
patient=load(source.cfg.patientFile,'cst');
cfg=merge_structs_local(default_config(),source.cfg); cfg.mode='integrated';
cfg.integrated.componentStage='layer'; cfg.lambdaEL=0; cfg.lambdaB=0;
cfg.lambdaS=0; cfg.lambdaES=0; cfg.lambdaT=0;
cfg.integrated.stationarityStoppingEnabled=false; cfg=configure_integrated_stage(cfg);
w0=reshape(double(source.denseResult.w),[],1);
assert(numel(w0)==double(source.dij.totalNumOfBixels));
[groups,cpGroups]=build_group_index(source.map);
doseModel=build_matrad_dose_model(source.dij,patient.cst,source.pln, ...
    cfg.integrated.constraintPenalty);
context=build_integrated_context(groups,cpGroups,source.map,source.machine,cfg);
[objective0,parts0]=evaluate_integrated_objective(w0,doseModel,context,cfg);
fullFdose0=parts0.dose; fullREL0=parts0.layerL21half;
lambdaScaleFull=fullFdose0/max(fullREL0,eps);
lambdaEL=alphaGrid*lambdaScaleFull;
scalingTable=table(alphaGrid,lambdaEL); writetable(scalingTable, ...
    fullfile(outputDir,'full_lambda_scaling.csv'));
setup=struct('artifactFile',artifactFile,'nWeights',numel(w0), ...
    'nLayers',numel(groups),'objective0',objective0,'fullFdose0',fullFdose0, ...
    'fullREL0',fullREL0,'lambdaScaleFull',lambdaScaleFull, ...
    'alphaGrid',alphaGrid,'lambdaEL',lambdaEL, ...
    'lambdaS',0,'lambdaES',0,'lambdaT',0,'lambdaB',0, ...
    'optimizationIterationsRun',0);
save(fullfile(outputDir,'full_layer_setup_complete.mat'),'setup','scalingTable','-v7.3');
outputs=struct('setup',setup,'scalingTable',scalingTable,'outputDir',outputDir);
fprintf('FULL SETUP ONLY: Fdose=%.12g R_EL=%.12g scale=%.12g; no optimization run.\n', ...
    fullFdose0,fullREL0,lambdaScaleFull);
end

function add_module_paths(root), f={'config','mapping','objective','delivery','integrated'}; for i=1:numel(f), addpath(fullfile(root,f{i})); end, end
function create_empty_output_directory(d), if exist(d,'dir')==7, e=dir(d); e=e(~ismember({e.name},{'.','..'})); if ~isempty(e), error('Output directory must be empty: %s',d); end, else, mkdir(d); end, end
function merged=merge_structs_local(base,override), merged=base; n=fieldnames(override); for i=1:numel(n), k=n{i}; if isfield(merged,k)&&isstruct(merged.(k))&&isstruct(override.(k)), merged.(k)=merge_structs_local(merged.(k),override.(k)); else, merged.(k)=override.(k); end, end, end
