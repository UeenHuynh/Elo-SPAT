function outputs = prepare_full_integrated_spot_20deg(outputDir, alphaSGrid)
% Calibrate dimensionless smooth-L0 spot scaling without optimization.

narginchk(0,2);
root=fileparts(fileparts(mfilename('fullpath'))); add_module_paths(root);
if nargin<1 || isempty(outputDir), outputDir=fullfile(root,'results', ...
        ['full_integrated_spot_setup_20deg_' datestr(now,'yyyymmdd_HHMMSS')]); end
if nargin<2 || isempty(alphaSGrid), alphaSGrid=[0.01;0.03;0.1]; end
alphaSGrid=reshape(double(alphaSGrid),[],1); create_empty_output_directory(outputDir);

artifactFile=fullfile(root,'results','task25_27_pipeline_20deg.mat');
frozenFile=fullfile(root,'results','full_integrated_layer_20deg_job22180', ...
    'alpha_0p3','alpha_0p3_final.mat');
expectedHash='9c77e38f26f43515fc1515a29523b94ec2c52265df74e9686f0196efc7097f18';
verify_sha256(frozenFile,expectedHash);
source=load(artifactFile,'cfg','machine','denseResult','dij','pln','map','groups');
layer=load(frozenFile,'record','final','finalDose','finalDelivery','cfgTrial');
patient=load(source.cfg.patientFile,'cst');
assert(abs(layer.record.alpha-0.3)<=1e-12 && layer.record.planAccepted, ...
    'Frozen Stage-1 reference is not the accepted alphaEL=0.3 plan.');

cfg=merge_structs_local(default_config(),source.cfg); cfg.mode='integrated';
cfg.integrated.componentStage='spot'; cfg.lambdaEL=layer.record.lambdaEL;
cfg.lambdaB=0; cfg.lambdaES=0; cfg.lambdaT=0; cfg.lambdaS=0;
cfg.integrated.stationarityStoppingEnabled=false; cfg=configure_integrated_stage(cfg);
w0=reshape(double(layer.final.w),[],1);
assert(numel(w0)==double(source.dij.totalNumOfBixels));
[groups,cpGroups]=build_group_index(source.map);
doseModel=build_matrad_dose_model(source.dij,patient.cst,source.pln, ...
    cfg.integrated.constraintPenalty);
context=build_integrated_context(groups,cpGroups,source.map,source.machine,cfg);
[objective0,parts0]=evaluate_integrated_objective(w0,doseModel,context,cfg);

spotScale=parts0.dose/max(parts0.spotSmoothL0,eps);
lambdaS=alphaSGrid*spotScale;
spotTermAtReference=lambdaS*parts0.spotSmoothL0;
relativeToDose=spotTermAtReference/max(parts0.dose,eps);
scalingTable=table(alphaSGrid,lambdaS,spotTermAtReference,relativeToDose);
wPositive=w0(w0>0);
setup=struct('artifactFile',artifactFile,'frozenFile',frozenFile, ...
    'frozenSha256',expectedHash,'nWeights',numel(w0),'nLayers',numel(groups), ...
    'alphaEL',layer.record.alpha,'lambdaEL',layer.record.lambdaEL, ...
    'objectiveLayerReference',objective0,'FdoseReference',parts0.dose, ...
    'RELReference',parts0.layerL21half,'RspotReference',parts0.spotSmoothL0, ...
    'epsS',cfg.epsS,'spotScale',spotScale,'alphaSGrid',alphaSGrid, ...
    'lambdaS',lambdaS,'exactZeroSpots',nnz(w0==0), ...
    'positiveSpots',nnz(w0>0),'activeSpots',nnz(w0>cfg.spotThreshold), ...
    'subMmuSpots',nnz(w0>0 & w0<cfg.mmuThreshold), ...
    'minimumPositiveWeight',min(wPositive),'mmuThreshold',cfg.mmuThreshold, ...
    'lambdaB',0,'lambdaES',0,'lambdaT',0,'optimizationIterationsRun',0);
writetable(scalingTable,fullfile(outputDir,'spot_lambda_scaling.csv'));
save(fullfile(outputDir,'spot_setup_complete.mat'),'setup','scalingTable','-v7.3');
outputs=struct('setup',setup,'scalingTable',scalingTable,'outputDir',outputDir);
fprintf(['SPOT SETUP ONLY: Fdose=%.12g Rspot=%.12g scale=%.12g; ' ...
    'zero=%d active=%d subMMU=%d minPositive=%.6g; no optimization run.\n'], ...
    parts0.dose,parts0.spotSmoothL0,spotScale,setup.exactZeroSpots, ...
    setup.activeSpots,setup.subMmuSpots,setup.minimumPositiveWeight);
end

function verify_sha256(file,expected)
[status,out]=system(sprintf('sha256sum "%s"',file));
assert(status==0,'Unable to checksum frozen Stage-1 artifact.');
actual=strtok(strtrim(out));
assert(strcmpi(actual,expected),'Frozen Stage-1 checksum mismatch.');
end
function add_module_paths(root), f={'config','mapping','objective','delivery','integrated'}; for i=1:numel(f), addpath(fullfile(root,f{i})); end, end
function create_empty_output_directory(d), if exist(d,'dir')==7, e=dir(d); e=e(~ismember({e.name},{'.','..'})); if ~isempty(e), error('Output directory must be empty: %s',d); end, else, mkdir(d); end, end
function merged=merge_structs_local(base,override), merged=base; n=fieldnames(override); for i=1:numel(n), k=n{i}; if isfield(merged,k)&&isstruct(merged.(k))&&isstruct(override.(k)), merged.(k)=merge_structs_local(merged.(k),override.(k)); else, merged.(k)=override.(k); end, end, end
