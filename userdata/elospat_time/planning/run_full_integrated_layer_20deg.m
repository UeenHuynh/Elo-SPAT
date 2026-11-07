function outputs = run_full_integrated_layer_20deg(outputDir, alphaGrid)
% Full 20-degree, layer-only integrated RMS-ELO-SPAT experiment.
%
% This runner reuses the preserved Dense Arc artifact. It does not rebuild
% dij, does not reduce its columns, and computes lambdaEL independently for
% the full problem as alpha*Fdose(w0)/R_EL(w0).

narginchk(0, 2);
root = fileparts(fileparts(mfilename('fullpath')));
add_module_paths(root);
if nargin < 1 || isempty(outputDir)
    outputDir = fullfile(root, 'results', ...
        ['full_integrated_layer_20deg_' datestr(now,'yyyymmdd_HHMMSS')]);
end
if nargin < 2 || isempty(alphaGrid)
    % Dimensionless, predeclared logarithmic coverage; not fitted to alpha=1.
    alphaGrid = [0.01; 0.03; 0.1; 0.3; 1];
end
alphaGrid = reshape(double(alphaGrid), [], 1);
assert(all(isfinite(alphaGrid) & alphaGrid > 0) && ...
    numel(unique(alphaGrid)) == numel(alphaGrid));
create_or_reuse_output_directory(outputDir);

artifactFile = fullfile(root, 'results', 'task25_27_pipeline_20deg.mat');
scalingFile = fullfile(root, 'results', ...
    'full_integrated_layer_setup_20deg_job22168', ...
    'full_layer_setup_complete.mat');
if exist(scalingFile,'file') ~= 2
    error('Missing full-problem scaling artifact from job 22168: %s',scalingFile);
end
scalingSource = load(scalingFile,'setup');
source = load(artifactFile, 'cfg','machine','denseResult','dij','pln','map','groups');
patient = load(source.cfg.patientFile, 'ct', 'cst');
cfg = merge_structs_local(default_config(), source.cfg);
cfg.mode = 'integrated'; cfg.arcSpacing = 20;
cfg.integrated.componentStage = 'layer';
cfg.lambdaB = 0; cfg.lambdaS = 0; cfg.lambdaES = 0; cfg.lambdaT = 0;
cfg.integrated.stationarityStoppingEnabled = false;
cfg.integrated.maxIterations = 50; % checkpoint block, not total budget
cfg.integrated.minimumIterations = 10;
cfg.integrated.verbose = true;
cfg.captureConsoleOutput = false;
cfg = configure_integrated_stage(cfg);
assert_layer_only(cfg);

nFull = double(source.dij.totalNumOfBixels);
w0 = reshape(double(source.denseResult.w), [], 1);
assert(numel(w0) == nFull && numel(source.map) == nFull);
[groups, cpGroups] = build_group_index(source.map);
assert(numel(groups) == numel(source.groups));
doseModel = build_matrad_dose_model(source.dij, patient.cst, source.pln, ...
    cfg.integrated.constraintPenalty);
context = build_integrated_context(groups, cpGroups, source.map, source.machine, cfg);
cfgZero = cfg; cfgZero.lambdaEL = 0;
[~, initialParts] = evaluate_integrated_objective(w0, doseModel, context, cfgZero);
fullFdose0 = initialParts.dose;
fullREL0 = initialParts.layerL21half;
lambdaScaleFull = fullFdose0 / max(fullREL0, eps);
scalingRelativeError = abs(lambdaScaleFull-scalingSource.setup.lambdaScaleFull) / ...
    max(1,abs(scalingSource.setup.lambdaScaleFull));
assert(scalingRelativeError <= 1e-11, ...
    'Full lambda scaling disagrees with job 22168 by %.3e.',scalingRelativeError);
lambdaGrid = zeros(size(alphaGrid));
for i = 1:numel(alphaGrid)
    q = find(abs(scalingSource.setup.alphaGrid-alphaGrid(i))<=1e-12,1);
    assert(~isempty(q),'Alpha %.12g is absent from job 22168 scaling.',alphaGrid(i));
    lambdaGrid(i) = scalingSource.setup.lambdaEL(q);
end
rule = default_stabilization_rule();
maxTotalIterations = 600;

D = nominal_matrix(source.dij);
denseDose = dose_metrics(D, w0, patient, cfg);
denseDelivery = compute_delivery_metrics(w0, source.map, groups, cfg, source.machine);
setup = struct('artifactFile',artifactFile,'nWeights',nFull, ...
    'nLayers',numel(groups),'alphaGrid',alphaGrid,'lambdaGrid',lambdaGrid, ...
    'fullFdose0',fullFdose0,'fullREL0',fullREL0, ...
    'lambdaScaleFull',lambdaScaleFull,'scalingSource',scalingFile, ...
    'scalingRelativeError',scalingRelativeError,'checkpointInterval',50, ...
    'maxTotalIterations',maxTotalIterations,'rule',rule, ...
    'penalties',struct('lambdaS',0,'lambdaES',0,'lambdaT',0,'lambdaB',0));
save(fullfile(outputDir, 'full_layer_setup.mat'), 'setup','-v7.3');
writetable(table(alphaGrid,lambdaGrid), ...
    fullfile(outputDir, 'full_lambda_scaling.csv'));

records = repmat(empty_record(), numel(alphaGrid), 1);
for a = 1:numel(alphaGrid)
    alpha = alphaGrid(a); token = alpha_token(alpha);
    cfgTrial = cfg; cfgTrial.lambdaEL = lambdaGrid(a);
    [state,bdtHistory,stabilization,resumeFile] = ...
        load_latest_checkpoint(outputDir,token,alpha,cfgTrial.lambdaEL);
    try
        while isempty(state) || (state.iterations < maxTotalIterations && ...
                ~stabilization.numerically_stabilized)
            if isempty(state)
                state = run_integrated_optimizer(doseModel,context,w0,cfgTrial);
            else
                state = run_integrated_optimizer(doseModel,context,state.w,cfgTrial,state);
            end
            delivery = compute_delivery_metrics(state.w,source.map,groups,cfgTrial,source.machine);
            bdtHistory.iteration(end+1,1) = state.iterations;
            bdtHistory.value(end+1,1) = delivery.T_total;
            stabilization = evaluate_numerical_stabilization( ...
                state.history,bdtHistory,rule);
            optimizer = state; %#ok<NASGU>
            save(fullfile(outputDir,sprintf('alpha_%s_iter_%04d.mat', ...
                token,state.iterations)),'optimizer','bdtHistory', ...
                'stabilization','cfgTrial','alpha','scalingFile','-v7.3');
            if strcmp(state.status,'line_search_failed') || ...
                    strcmp(state.status,'numerical_stagnation')
                break;
            end
        end

        record = empty_record(); record.alpha = alpha;
        record.lambdaEL = cfgTrial.lambdaEL; record.iterations = state.iterations;
        record.optimizerStatus = state.status; record.stoppingReason = state.stoppingReason;
        record.objective = state.objective;
        record.normalizedStationarity = state.finalNormalizedStationarity;
        record.numericallyStabilized = stabilization.numerically_stabilized;
        record.stabilizationReason = stabilization.reason;
        record.exactZeroLayers = state.history.exactZeroLayers(end);
        record.activeLayers = state.history.exactActiveLayers(end);
        record.activeSpots = state.history.activeSpots(end);
        record.integratedBDT = bdtHistory.value(end);
        record.bdtChangeFromDensePercent = pct(record.integratedBDT,denseDelivery.T_total);
        record.maxConstraintViolation = state.history.maxConstraintViolation(end);

        if stabilization.numerically_stabilized
            extracted = extract_integrated_support(state.w,w0,context,cfgTrial);
            supportMask = reshape(logical(extracted.mask),[],1);
            final = run_final_reoptimization(source.dij,patient.cst,source.pln, ...
                extracted.w,cfgTrial,supportMask,extracted.layerMask);
            supportPreserved = ~any(abs(final.w(~supportMask)) > ...
                rule.supportLeakageTolerance);
            finalDose = dose_metrics(D,final.w,patient,cfgTrial);
            finalDelivery = compute_delivery_metrics(final.w,source.map,groups,cfgTrial,source.machine);
            finalSolver = compute_solver_metrics(final.resultGUI,struct(),final.optimizer);
            [dosePass,doseTable,worstAdverse] = dose_quality_detail( ...
                finalDose,denseDose,cfg.evaluation.ctvD95TolerancePercent);
            ipoptAcceptable = finalSolver.converged || ...
                finalSolver.acceptableButNotConverged || ...
                contains(lower(finalSolver.statusMessage),'acceptable');
            acceptance = evaluate_plan_acceptance(dosePass, ...
                record.maxConstraintViolation,supportPreserved,ipoptAcceptable,rule);
            record.planAccepted = acceptance.plan_accepted;
            record.doseQualityPass = dosePass;
            record.supportPreserved = supportPreserved;
            record.finalIpoptStatus = finalSolver.statusMessage;
            record.finalBDT = finalDelivery.T_total;
            record.bdtChangeFromDensePercent = pct(record.finalBDT,denseDelivery.T_total);
            record.finalActiveLayers = finalDelivery.nActiveLayers;
            record.finalActiveSpots = finalDelivery.nActiveSpots;
            targetIx = find([finalDose.structures.isTarget],1);
            if ~isempty(targetIx)
                record.targetD95 = finalDose.structures(targetIx).D95;
                record.targetD98 = finalDose.structures(targetIx).D98;
            end
            record.worstAdverseDoseChangePercent = worstAdverse;
            writetable(doseTable,fullfile(outputDir,sprintf( ...
                'alpha_%s_final_ipopt_dose.csv',token)));
            save(fullfile(outputDir,sprintf('alpha_%s_final.mat',token)), ...
                'record','optimizer','bdtHistory','stabilization','extracted', ...
                'supportMask','final','finalDose','finalDelivery','finalSolver', ...
                'acceptance','cfgTrial','scalingFile','resumeFile','-v7.3');
        end
        records(a) = record;
    catch err
        records(a) = empty_record(); records(a).alpha = alpha;
        records(a).lambdaEL = cfgTrial.lambdaEL;
        records(a).runStatus = ['failed: ' err.identifier];
        errorInfo = struct('identifier',err.identifier,'message',err.message, ...
            'stack',err.stack); %#ok<NASGU>
        save(fullfile(outputDir,sprintf('alpha_%s_failure.mat',token)), ...
            'errorInfo','cfgTrial','-v7.3');
        warning('Full layer trial alpha=%g failed: %s',alpha,err.message);
    end
    writetable(struct2table(records(1:a),'AsArray',true), ...
        fullfile(outputDir,'full_layer_summary.csv'));
end

save(fullfile(outputDir,'full_layer_complete.mat'), ...
    'records','setup','denseDose','denseDelivery','-v7.3');
outputs = struct('records',records,'setup',setup,'outputDir',outputDir);
end

function assert_layer_only(cfg)
assert(strcmpi(cfg.integrated.componentStage,'layer'));
assert(cfg.lambdaB==0 && cfg.lambdaS==0 && cfg.lambdaES==0 && cfg.lambdaT==0);
assert(~cfg.integrated.stationarityStoppingEnabled);
end

function r = empty_record()
r = struct('alpha',NaN,'lambdaEL',NaN,'iterations',NaN, ...
    'optimizerStatus','','stoppingReason','','objective',NaN, ...
    'normalizedStationarity',NaN,'numericallyStabilized',false, ...
    'stabilizationReason','','exactZeroLayers',NaN,'activeLayers',NaN, ...
    'activeSpots',NaN,'integratedBDT',NaN,'maxConstraintViolation',NaN, ...
    'bdtChangeFromDensePercent',NaN,'finalActiveLayers',NaN, ...
    'finalActiveSpots',NaN,'targetD95',NaN,'targetD98',NaN, ...
    'planAccepted',false,'doseQualityPass',false,'supportPreserved',false, ...
    'finalIpoptStatus','','finalBDT',NaN, ...
    'worstAdverseDoseChangePercent',NaN,'runStatus','completed');
end

function [pass,tbl,worstAdverse] = dose_quality_detail(candidate,dense,tolerance)
n=numel(dense.structures); rows=repmat(struct('structure','','type','', ...
    'denseD98',NaN,'finalD98',NaN,'denseD95',NaN,'finalD95',NaN, ...
    'denseDmean',NaN,'finalDmean',NaN,'denseDmax',NaN,'finalDmax',NaN, ...
    'adverseChangePercent',NaN,'structurePass',false),n,1);
worstAdverse=0;
for i=1:n
    d=dense.structures(i); f=candidate.structures(i);
    rows(i).structure=d.name; rows(i).type=d.type;
    rows(i).denseD98=d.D98; rows(i).finalD98=f.D98;
    rows(i).denseD95=d.D95; rows(i).finalD95=f.D95;
    rows(i).denseDmean=d.Dmean; rows(i).finalDmean=f.Dmean;
    rows(i).denseDmax=d.Dmax; rows(i).finalDmax=f.Dmax;
    if d.isTarget
        adverse=max([0,-pct(f.D95,d.D95),-pct(f.D98,d.D98)]);
    elseif strcmpi(d.type,'OAR')
        adverse=max([0,pct(f.Dmean,d.Dmean),pct(f.Dmax,d.Dmax)]);
    else
        adverse=0;
    end
    rows(i).adverseChangePercent=adverse;
    rows(i).structurePass=adverse<=tolerance; worstAdverse=max(worstAdverse,adverse);
end
tbl=struct2table(rows); pass=all([rows.structurePass]);
end

function value=pct(x,ref), value=100*(double(x)-double(ref))/max(abs(double(ref)),eps); end
function D=nominal_matrix(dij), if iscell(dij.physicalDose), D=dij.physicalDose{1}; else, D=dij.physicalDose; end, end
function m=dose_metrics(D,w,patient,cfg), m=compute_dvh_metrics(full(D*reshape(double(w),[],1)),patient.ct,patient.cst,cfg); end
function token=alpha_token(alpha), token=strrep(sprintf('%.4g',alpha),'.','p'); end
function add_module_paths(root), f={'config','planning','mapping','objective','delivery','reduction','robustness','evaluation','integrated'}; for i=1:numel(f), addpath(fullfile(root,f{i})); end, end
function create_or_reuse_output_directory(d), if exist(d,'dir')~=7, mkdir(d); end, end
function merged=merge_structs_local(base,override), merged=base; n=fieldnames(override); for i=1:numel(n), k=n{i}; if isfield(merged,k)&&isstruct(merged.(k))&&isstruct(override.(k)), merged.(k)=merge_structs_local(merged.(k),override.(k)); else, merged.(k)=override.(k); end, end, end

function [state,bdtHistory,stabilization,resumeFile] = ...
        load_latest_checkpoint(outputDir,token,alpha,lambdaEL)
state=[]; bdtHistory=struct('iteration',zeros(0,1),'value',zeros(0,1));
stabilization=struct('numerically_stabilized',false,'reason','not_started');
resumeFile='';
files=dir(fullfile(outputDir,sprintf('alpha_%s_iter_*.mat',token)));
if isempty(files), return; end
iterations=zeros(numel(files),1);
for i=1:numel(files)
    parsed=regexp(files(i).name,'_iter_(\d+)\.mat$','tokens','once');
    if ~isempty(parsed), iterations(i)=str2double(parsed{1}); end
end
[~,q]=max(iterations); resumeFile=fullfile(files(q).folder,files(q).name);
saved=load(resumeFile,'optimizer','bdtHistory','stabilization','cfgTrial','alpha');
assert(abs(saved.alpha-alpha)<=1e-12,'Checkpoint alpha mismatch.');
assert(abs(saved.cfgTrial.lambdaEL-lambdaEL)<=1e-12*max(1,abs(lambdaEL)), ...
    'Checkpoint lambdaEL mismatch.');
state=saved.optimizer; bdtHistory=saved.bdtHistory; stabilization=saved.stabilization;
fprintf('RESUME alpha=%g from iteration %d: %s\n',alpha,state.iterations,resumeFile);
end
