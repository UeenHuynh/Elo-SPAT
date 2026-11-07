function [pln, stf, dij, resultGUI, meta] = build_impt_plan(ct, cst, cfg)
% build_impt_plan  IMPT baseline planner (no RMS-ELO-SPAT novelty).
%
% This function replicates matRad's standard IMPT workflow to create a
% fair, reproducible baseline for comparison with RMS-ELO-SPAT methods.
%
% BASELINE-ONLY: This function does NOT implement:
%   - ELO-SPAT L2,1/2 group sparsity
%   - Energy sequencing penalty
%   - BDT-aware optimization
%   - Spot/layer filtering
%   - MMU filtering
%   - Robust sparse optimization
%
% INPUTS
%   ct  - matRad CT structure (ct.resolution, ct.x, ct.y, ct.z)
%   cst - matRad clinical structure table (objectives unchanged)
%   cfg - Configuration struct with fields:
%         .imptAngles      - [2-5] gantry angles (degrees)
%         .radiationMode   - 'protons'
%         .machine         - matRad machine name
%         .bioModel        - biological model ('none', 'constRBE', etc.)
%         .numOfFractions  - number of fractions (default: 1)
%         .bixelWidth      - bixel width in mm (default: 5)
%         .doseGrid        - dose grid resolution (default: ct.resolution)
%         .isoCenter       - isocenter [x,y,z] (default: auto)
%
% OUTPUTS
%   pln       - matRad plan structure
%   stf       - matRad steering information
%   dij       - matRad dose influence matrix
%   resultGUI - matRad optimization result (weights, dose, solver info)
%   meta      - metadata struct for reproducibility and comparison:
%               .methodName, .patientFile, .gantryAngles, .couchAngles,
%               .bixelWidth, .doseGrid, .machine, .bioModel, .numOfFractions,
%               .nBixels, .nBeams, .runtime, .optimizerName, .optimizerOptions,
%               .optimizerResultInfo, .converged, .statusMessage, .timestamp
%
% CONVERGENCE
%   - converged = true:  IPOPT status 0 or 1 (solved)
%   - converged = false: max iterations, infeasible, or error
%
% USAGE
%   cfg = default_config();
%   load('PROSTATE.mat', 'ct', 'cst');
%   [pln, stf, dij, resultGUI, meta] = build_impt_plan(ct, cst, cfg);
%
% =========================================================================
% RESEARCH-ONLY DISCLAIMER:
%   This function is part of RMS-ELO-SPAT, a research-only module.
%   Not for clinical use. matRad is used as a research sandbox only.
% =========================================================================

% -------------------------------------------------------------------------
% Initialize metadata
% -------------------------------------------------------------------------
meta = struct();
meta.methodName = 'IMPT baseline';
meta.timestamp = datetime('now');
startTime = tic;
pln = struct();
stf = struct();
dij = struct();
resultGUI = struct();
optimizer = [];
totalRays = 0;
totalBixels = 0;
currentStage = 'initialization';
lastSuccessfulStage = '';
runCtx = struct();
persistEnabled = false;
consoleCaptureState = struct('active', false, 'logFile', '', 'startedAt', '');

% -------------------------------------------------------------------------
% Input validation
% -------------------------------------------------------------------------
narginchk(3, 3);

if ~isstruct(ct) || ~isfield(ct, 'resolution')
    error('build_impt_plan:invalidInput', ...
        'ct must be a valid matRad CT structure with ct.resolution field.');
end

if ~iscell(cst) || size(cst, 2) < 6
    error('build_impt_plan:invalidInput', ...
        'cst must be a valid matRad cell array with at least 6 columns.');
end

if ~isstruct(cfg) || ~isfield(cfg, 'imptAngles')
    error('build_impt_plan:invalidInput', ...
        'cfg must be a valid configuration struct with cfg.imptAngles field.');
end

% Validate IMPT angles (2-5 angles for baseline)
nAngles = numel(cfg.imptAngles);
if nAngles < 2 || nAngles > 5
    error('build_impt_plan:invalidAngles', ...
        'cfg.imptAngles must contain 2 to 5 gantry angles (got %d).', nAngles);
end

% Set defaults for required fields (before fprintf)
if ~isfield(cfg, 'radiationMode') || isempty(cfg.radiationMode)
    cfg.radiationMode = 'protons';
end

if ~isfield(cfg, 'machine') || isempty(cfg.machine)
    cfg.machine = 'Generic';
end

if ~isfield(cfg, 'bioModel') || isempty(cfg.bioModel)
    cfg.bioModel = 'constRBE';
end

if ~isfield(cfg, 'captureConsoleOutput') || isempty(cfg.captureConsoleOutput)
    cfg.captureConsoleOutput = true;
end

% Set defaults for optional fields
if ~isfield(cfg, 'numOfFractions') || isempty(cfg.numOfFractions)
    cfg.numOfFractions = 1;
end

if ~isfield(cfg, 'bixelWidth') || isempty(cfg.bixelWidth)
    cfg.bixelWidth = 5;  % mm, matRad default
end

if ~isfield(cfg, 'doseGrid') || isempty(cfg.doseGrid)
    cfg.doseGrid = ct.resolution;  % use CT resolution
else
    % Handle empty fields in doseGrid struct
    if ~isfield(cfg.doseGrid, 'x') || isempty(cfg.doseGrid.x)
        cfg.doseGrid.x = ct.resolution.x;
    end
    if ~isfield(cfg.doseGrid, 'y') || isempty(cfg.doseGrid.y)
        cfg.doseGrid.y = ct.resolution.y;
    end
    if ~isfield(cfg.doseGrid, 'z') || isempty(cfg.doseGrid.z)
        cfg.doseGrid.z = ct.resolution.z;
    end
end

% Store patient file if available
if isfield(cfg, 'patientFile')
    meta.patientFile = cfg.patientFile;
else
    meta.patientFile = 'unknown';
end

[runCtx, persistEnabled] = safe_create_run_context(cfg, 'impt_baseline');
if persistEnabled
    consoleCaptureState = safe_start_console_capture(runCtx, cfg.captureConsoleOutput);
    consoleCaptureCleanup = onCleanup(@() safe_stop_console_capture(consoleCaptureState)); %#ok<NASGU>
    safe_append_run_log(runCtx, sprintf('build_impt_plan started for patient %s', meta.patientFile), 'INFO');
    meta.runId = runCtx.runId;
else
    meta.runId = '';
end
meta.outputPaths = build_output_paths(runCtx);

% -------------------------------------------------------------------------
% Configure matRad plan structure (standard IMPT workflow)
% -------------------------------------------------------------------------
fprintf('\n=== build_impt_plan: IMPT Baseline (no RMS-ELO-SPAT novelty) ===\n');
fprintf('  Gantry angles: [%s] degrees\n', ...
    strjoin(arrayfun(@num2str, cfg.imptAngles, 'UniformOutput', false), ', '));
fprintf('  Radiation mode: %s\n', cfg.radiationMode);
fprintf('  Machine: %s\n', cfg.machine);
fprintf('  Biological model: %s\n', cfg.bioModel);
fprintf('  Number of fractions: %d\n', cfg.numOfFractions);
fprintf('  Bixel width: %.1f mm\n', cfg.bixelWidth);

pln = struct();

% Radiation modality
pln.radiationMode = cfg.radiationMode;  % 'protons'
pln.machine       = cfg.machine;        % 'Generic' or machine name

% Biological model (matRad convention: pln.bioModel, not pln.propOpt.bioOptimization)
if isfield(cfg, 'bioModel') && ~isempty(cfg.bioModel)
    pln.bioModel = cfg.bioModel;
else
    pln.bioModel = 'none';
end

% Optimization quantity (for constRBE, optimize RBE-weighted dose)
if strcmp(pln.bioModel, 'constRBE')
    pln.propOpt.quantityOpt = 'RBExDose';
else
    pln.propOpt.quantityOpt = 'physicalDose';
end

% Scenario: nominal only for baseline (matRad convention: string, not struct)
pln.multScen = 'nomScen';

% Number of fractions
pln.numOfFractions = cfg.numOfFractions;

% Beam geometry — IMPT uses fixed gantry angles with couch = 0
pln.propStf.gantryAngles = cfg.imptAngles(:)';  % row vector
pln.propStf.couchAngles  = zeros(1, nAngles);   % coplanar beams
pln.propStf.bixelWidth   = cfg.bixelWidth;      % mm

% Isocenter (use matRad_getIsoCenter for robustness)
if isfield(cfg, 'isoCenter') && ~isempty(cfg.isoCenter)
    pln.propStf.isoCenter = cfg.isoCenter;
    fprintf('  Isocenter: [%.1f, %.1f, %.1f] mm (from config)\n', ...
        pln.propStf.isoCenter(1), pln.propStf.isoCenter(2), pln.propStf.isoCenter(3));
else
    % matRad standard: compute isocenter from target structures
    pln.propStf.isoCenter = matRad_getIsoCenter(cst, ct, 0);
    fprintf('  Isocenter: [%.1f, %.1f, %.1f] mm (auto-computed)\n', ...
        pln.propStf.isoCenter(1), pln.propStf.isoCenter(2), pln.propStf.isoCenter(3));
end

% Dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = cfg.doseGrid.x;
pln.propDoseCalc.doseGrid.resolution.y = cfg.doseGrid.y;
pln.propDoseCalc.doseGrid.resolution.z = cfg.doseGrid.z;

% Dose calculation engine (explicit for reproducibility)
if isfield(cfg, 'doseCalc') && isfield(cfg.doseCalc, 'engine')
    pln.propDoseCalc.engine = cfg.doseCalc.engine;
else
    pln.propDoseCalc.engine = 'HongPB';  % matRad proton pencil beam (Hong 2014)
end
pln.propDoseCalc.calcLET = 0;      % do not calculate LET for baseline

% Optimization settings (explicit for reproducibility)
pln.propOpt.optimizer = 'IPOPT';   % matRad default optimizer
pln.propOpt.runDAO = 0;            % no direct aperture optimization
pln.propSeq.runSequencing = 0;     % no sequencing for baseline

% Store metadata (before optimization, so we have config even if it fails)
meta.gantryAngles = pln.propStf.gantryAngles;
meta.couchAngles = pln.propStf.couchAngles;
meta.bixelWidth = pln.propStf.bixelWidth;
meta.isoCenter = pln.propStf.isoCenter;
meta.doseGrid = cfg.doseGrid;
meta.machine = pln.machine;
meta.bioModel = pln.bioModel;
meta.quantityOpt = pln.propOpt.quantityOpt;
meta.numOfFractions = pln.numOfFractions;
meta.optimizerName = pln.propOpt.optimizer;
meta.doseEngine = pln.propDoseCalc.engine;
pln.rmsPersist = build_output_paths(runCtx);

fprintf('  Optimizer: IPOPT (matRad defaults, no custom tuning for baseline)\n');
currentStage = 'plan_configured';
lastSuccessfulStage = 'plan_configured';
safe_save_run_summary(runCtx, build_impt_summary( ...
    meta, runCtx, currentStage, lastSuccessfulStage, pln, stf, dij, resultGUI, ...
    totalRays, totalBixels, 'running', ''));

% -------------------------------------------------------------------------
% Generate steering information (matRad standard workflow)
% -------------------------------------------------------------------------
fprintf('\n--- Generating beam steering information (matRad_generateStf) ---\n');
currentStage = 'stf_generation';
safe_append_run_log(runCtx, 'Starting STF generation.', 'INFO');
safe_save_run_summary(runCtx, build_impt_summary( ...
    meta, runCtx, currentStage, lastSuccessfulStage, pln, stf, dij, resultGUI, ...
    totalRays, totalBixels, 'running', ''));

try
    stf = matRad_generateStf(ct, cst, pln);
    fprintf('[OK] Steering information generated successfully.\n');
    fprintf('     Number of beams: %d\n', numel(stf));
    
    % Count total rays and bixels
    totalRays = 0;
    totalBixels = 0;
    for b = 1:numel(stf)
        totalRays = totalRays + numel(stf(b).ray);
        for r = 1:numel(stf(b).ray)
            totalBixels = totalBixels + numel(stf(b).ray(r).energy);
        end
    end
    fprintf('     Total rays: %d\n', totalRays);
    fprintf('     Total bixels (spots): %d\n', totalBixels);
    
    meta.nBeams = numel(stf);
    meta.nBixels = totalBixels;
    lastSuccessfulStage = 'stf_generation';
    safe_append_run_log(runCtx, sprintf( ...
        'STF generation completed: %d beams, %d rays, %d bixels.', ...
        numel(stf), totalRays, totalBixels), 'INFO');
    safe_save_checkpoint(runCtx, 'stf_summary', struct( ...
        'nBeams', numel(stf), ...
        'totalRays', totalRays, ...
        'totalBixels', totalBixels, ...
        'gantryAngles', pln.propStf.gantryAngles, ...
        'couchAngles', pln.propStf.couchAngles), 1);
    safe_save_run_summary(runCtx, build_impt_summary( ...
        meta, runCtx, currentStage, lastSuccessfulStage, pln, stf, dij, resultGUI, ...
        totalRays, totalBixels, 'running', ''));
    
catch ME
    meta.converged = false;
    meta.exitFlag = NaN;
    meta.exitStatus = NaN;
    meta.exitMsg = ME.message;
    meta.exitMessage = ME.message;
    meta.statusMessage = sprintf('STF generation failed: %s', ME.message);
    meta.errorMessage = ME.message;
    meta.runtime = toc(startTime);
    
    fprintf('[FAILED] matRad_generateStf error: %s\n', ME.message);
    safe_append_run_log(runCtx, meta.statusMessage, 'ERROR');
    safe_save_checkpoint(runCtx, 'stf_failure', struct( ...
        'errorMessage', ME.message, ...
        'runtimeSeconds', meta.runtime, ...
        'currentStage', currentStage), 1);
    safe_save_run_summary(runCtx, build_impt_summary( ...
        meta, runCtx, currentStage, lastSuccessfulStage, pln, stf, dij, resultGUI, ...
        totalRays, totalBixels, 'failed', ME.message));
    
    % Return empty outputs with metadata
    pln = struct();
    stf = struct();
    dij = struct();
    resultGUI = struct( ...
        'converged', false, ...
        'exitFlag', NaN, ...
        'exitStatus', NaN, ...
        'exitMsg', ME.message, ...
        'exitMessage', ME.message, ...
        'statusMessage', meta.statusMessage);
    return;
end

% -------------------------------------------------------------------------
% Calculate dose influence matrix (matRad standard workflow)
% -------------------------------------------------------------------------
fprintf('\n--- Computing dose influence matrix (matRad_calcDoseInfluence) ---\n');
currentStage = 'dij_calculation';
safe_append_run_log(runCtx, 'Starting dij calculation.', 'INFO');
safe_save_checkpoint(runCtx, 'before_dij', struct( ...
    'nBeams', meta.nBeams, ...
    'totalRays', totalRays, ...
    'totalBixels', totalBixels, ...
    'doseGrid', cfg.doseGrid, ...
    'machine', pln.machine), 1);
safe_save_run_summary(runCtx, build_impt_summary( ...
    meta, runCtx, currentStage, lastSuccessfulStage, pln, stf, dij, resultGUI, ...
    totalRays, totalBixels, 'running', ''));

try
    dij = matRad_calcDoseInfluence(ct, cst, stf, pln);
    fprintf('[OK] Dose influence matrix computed successfully.\n');
    fprintf('     Total number of bixels: %d\n', dij.totalNumOfBixels);
    fprintf('     Number of voxels: %d\n', dij.doseGrid.numOfVoxels);
    
    % Check for multi-scenario dij
    if iscell(dij.physicalDose)
        nScenarios = numel(dij.physicalDose);
        fprintf('     Multi-scenario dij detected: %d scenarios\n', nScenarios);
    else
        fprintf('     Single-scenario (nominal) dij\n');
    end
    
    % Update metadata
    meta.nBixels = dij.totalNumOfBixels;
    lastSuccessfulStage = 'dij_calculation';
    safe_append_run_log(runCtx, sprintf( ...
        'Dij calculation completed: %d bixels, %d voxels.', ...
        dij.totalNumOfBixels, dij.doseGrid.numOfVoxels), 'INFO');
    safe_save_checkpoint(runCtx, 'dij_summary', struct( ...
        'totalNumOfBixels', dij.totalNumOfBixels, ...
        'numOfVoxels', dij.doseGrid.numOfVoxels, ...
        'isMultiScenario', iscell(dij.physicalDose)), 1);
    safe_save_run_summary(runCtx, build_impt_summary( ...
        meta, runCtx, currentStage, lastSuccessfulStage, pln, stf, dij, resultGUI, ...
        totalRays, totalBixels, 'running', ''));
    
catch ME
    meta.converged = false;
    meta.exitFlag = NaN;
    meta.exitStatus = NaN;
    meta.exitMsg = ME.message;
    meta.exitMessage = ME.message;
    meta.statusMessage = sprintf('Dij calculation failed: %s', ME.message);
    meta.errorMessage = ME.message;
    meta.runtime = toc(startTime);
    
    fprintf('[FAILED] matRad_calcDoseInfluence error: %s\n', ME.message);
    safe_append_run_log(runCtx, meta.statusMessage, 'ERROR');
    safe_save_checkpoint(runCtx, 'dij_failure', struct( ...
        'errorMessage', ME.message, ...
        'runtimeSeconds', meta.runtime, ...
        'nBeams', meta.nBeams, ...
        'totalRays', totalRays, ...
        'totalBixels', totalBixels), 1);
    safe_save_run_summary(runCtx, build_impt_summary( ...
        meta, runCtx, currentStage, lastSuccessfulStage, pln, stf, dij, resultGUI, ...
        totalRays, totalBixels, 'failed', ME.message));
    
    % Return partial outputs with metadata
    dij = struct();
    resultGUI = struct( ...
        'converged', false, ...
        'exitFlag', NaN, ...
        'exitStatus', NaN, ...
        'exitMsg', ME.message, ...
        'exitMessage', ME.message, ...
        'statusMessage', meta.statusMessage);
    return;
end

% -------------------------------------------------------------------------
% Optimize spot weights (matRad standard fluence optimization)
% -------------------------------------------------------------------------
fprintf('\n--- Running fluence optimization (matRad_fluenceOptimization) ---\n');
fprintf('    NOTE: Using matRad objectives from cst (no RMS-ELO-SPAT penalties)\n');
currentStage = 'fluence_optimization';
safe_append_run_log(runCtx, 'Starting fluence optimization.', 'INFO');
safe_save_checkpoint(runCtx, 'before_optimization', struct( ...
    'totalNumOfBixels', meta.nBixels, ...
    'numOfVoxels', safe_get_dij_voxels(dij), ...
    'optimizer', pln.propOpt.optimizer), 1);
safe_save_run_summary(runCtx, build_impt_summary( ...
    meta, runCtx, currentStage, lastSuccessfulStage, pln, stf, dij, resultGUI, ...
    totalRays, totalBixels, 'running', ''));

try
    % Get both resultGUI and optimizer info
    [resultGUI, optimizer] = matRad_fluenceOptimization(dij, cst, pln);
    fprintf('[OK] Fluence optimization completed.\n');
    
    % Store optimizer metadata (optimizer may be object, not struct)
    if exist('optimizer', 'var') && ~isempty(optimizer)
        % Try isprop first (for objects), then isfield (for structs)
        if isprop(optimizer, 'options') || isfield(optimizer, 'options')
            meta.optimizerOptions = optimizer.options;
        else
            meta.optimizerOptions = struct();
        end
        
        if isprop(optimizer, 'resultInfo') || isfield(optimizer, 'resultInfo')
            meta.optimizerResultInfo = optimizer.resultInfo;
        else
            meta.optimizerResultInfo = struct();
        end
    else
        meta.optimizerOptions = struct();
        meta.optimizerResultInfo = struct();
    end
    lastSuccessfulStage = 'fluence_optimization';
    safe_append_run_log(runCtx, 'Fluence optimization completed.', 'INFO');
    safe_save_checkpoint(runCtx, 'optimization_summary', struct( ...
        'weightsCount', safe_numel_field(resultGUI, 'w'), ...
        'hasDoseCube', isfield(resultGUI, 'physicalDose') && ~isempty(resultGUI.physicalDose), ...
        'objectiveValue', safe_numeric_field(resultGUI, 'fval')), 1);
     
catch ME
    meta.converged = false;
    meta.exitFlag = NaN;
    meta.exitStatus = NaN;
    meta.exitMsg = ME.message;
    meta.exitMessage = ME.message;
    meta.statusMessage = sprintf('Optimization failed: %s', ME.message);
    meta.errorMessage = ME.message;
    meta.runtime = toc(startTime);
    
    fprintf('[FAILED] matRad_fluenceOptimization error: %s\n', ME.message);
    safe_append_run_log(runCtx, meta.statusMessage, 'ERROR');
    safe_save_checkpoint(runCtx, 'optimization_failure', struct( ...
        'errorMessage', ME.message, ...
        'runtimeSeconds', meta.runtime, ...
        'weightsCount', safe_numel_field(resultGUI, 'w')), 1);
    safe_save_run_summary(runCtx, build_impt_summary( ...
        meta, runCtx, currentStage, lastSuccessfulStage, pln, stf, dij, resultGUI, ...
        totalRays, totalBixels, 'failed', ME.message));
    
    % Return partial outputs with metadata
    resultGUI = struct( ...
        'converged', false, ...
        'exitFlag', NaN, ...
        'exitStatus', NaN, ...
        'exitMsg', ME.message, ...
        'exitMessage', ME.message, ...
        'statusMessage', meta.statusMessage);
    return;
end

% -------------------------------------------------------------------------
% Check convergence status (strict: only solved = converged)
% -------------------------------------------------------------------------
fprintf('\n--- Checking convergence status ---\n');

% IPOPT status codes (from optimizer.resultInfo.status):
%   0, 1: solved → CONVERGED
%   -1: max iterations reached → NOT CONVERGED
%   < -1: failure (infeasibility, numerical error) → NOT CONVERGED

converged = false;
status = NaN;
statusMsg = 'unknown';

% Try to extract status from optimizer.resultInfo (optimizer may be object)
if exist('optimizer', 'var') && ~isempty(optimizer)
    % Try isprop first (for objects), then isfield (for structs)
    if (isprop(optimizer, 'resultInfo') || isfield(optimizer, 'resultInfo'))
        resultInfo = optimizer.resultInfo;
        
        if (isprop(resultInfo, 'status') || isfield(resultInfo, 'status'))
            status = resultInfo.status;
        end
        
        if (isprop(resultInfo, 'message') || isfield(resultInfo, 'message'))
            statusMsg = resultInfo.message;
        elseif (isprop(resultInfo, 'statusMsg') || isfield(resultInfo, 'statusMsg'))
            statusMsg = resultInfo.statusMsg;
        end
    end
end

% Fallback: try resultGUI.info (older matRad versions)
if isnan(status) && isfield(resultGUI, 'info')
    if isfield(resultGUI.info, 'exitflag')
        status = resultGUI.info.exitflag;
    elseif isfield(resultGUI.info, 'status')
        status = resultGUI.info.status;
    end
    
    if isfield(resultGUI.info, 'exitMsg')
        statusMsg = resultGUI.info.exitMsg;
    elseif isfield(resultGUI.info, 'message')
        statusMsg = resultGUI.info.message;
    end
end

% Determine convergence (strict criteria)
if ~isnan(status)
    if status == 0 || status == 1
        converged = true;
        fprintf('[CONVERGED] IPOPT status: %d (%s)\n', status, statusMsg);
        meta.statusMessage = sprintf('Converged (status %d)', status);
    elseif status == -1
        converged = false;
        fprintf('[NOT CONVERGED] Maximum iterations reached (status: %d)\n', status);
        fprintf('                Plan stored but NOT marked as converged.\n');
        fprintf('                Message: %s\n', statusMsg);
        meta.statusMessage = sprintf('Max iterations reached (status %d)', status);
    else
        converged = false;
        fprintf('[FAILED] Optimization failed (status: %d)\n', status);
        fprintf('         Message: %s\n', statusMsg);
        meta.statusMessage = sprintf('Failed (status %d): %s', status, statusMsg);
        
        % Extract diagnostics
        if isfield(resultGUI, 'fval')
            fprintf('         Final objective value: %.6e\n', resultGUI.fval);
        end
        if isfield(meta.optimizerResultInfo, 'iterations')
            fprintf('         Iterations completed: %d\n', meta.optimizerResultInfo.iterations);
        end
    end
else
    % No status available — check for weights as last resort
    if isfield(resultGUI, 'w') && ~isempty(resultGUI.w)
        fprintf('[WARNING] Convergence status unknown (no status code found).\n');
        fprintf('          Optimized weights are available (numel(w) = %d).\n', numel(resultGUI.w));
        fprintf('          Marking as NOT converged for safety.\n');
        converged = false;
        meta.statusMessage = 'Unknown status (no status code, weights available)';
    else
        converged = false;
        meta.statusMessage = 'Failed (no weights available)';
        fprintf('[FAILED] No weights available and no status code.\n');
    end
end

% Store convergence metadata
meta.converged = converged;
meta.exitFlag = status;
meta.exitStatus = status;
meta.exitMsg = statusMsg;
meta.exitMessage = statusMsg;

% Mirror the most useful convergence fields onto resultGUI for
% compatibility with earlier tests/documentation.
resultGUI.converged = converged;
resultGUI.exitFlag = status;
resultGUI.exitStatus = status;
resultGUI.exitMsg = statusMsg;
resultGUI.exitMessage = statusMsg;
resultGUI.statusMessage = meta.statusMessage;

% -------------------------------------------------------------------------
% Finalize metadata and output summary
% -------------------------------------------------------------------------
meta.runtime = toc(startTime);
currentStage = 'completed';
lastSuccessfulStage = 'completed';

fprintf('\n--- Output summary ---\n');
fprintf('  Method: %s\n', meta.methodName);
fprintf('  Converged: %s\n', mat2str(meta.converged));
fprintf('  Status: %s\n', meta.statusMessage);

if isfield(resultGUI, 'w')
    fprintf('  Optimized weights: %d values (resultGUI.w)\n', numel(resultGUI.w));
end

if isfield(resultGUI, 'physicalDose') && ~isempty(resultGUI.physicalDose)
    fprintf('  Dose cube: available (resultGUI.physicalDose)\n');
else
    fprintf('  Dose cube: not computed (call matRad_calcCubes if needed)\n');
end

if isfield(resultGUI, 'fval')
    fprintf('  Final objective value: %.6e\n', resultGUI.fval);
end

if isfield(meta.optimizerResultInfo, 'iterations')
    fprintf('  Solver iterations: %d\n', meta.optimizerResultInfo.iterations);
elseif isfield(resultGUI, 'info') && isfield(resultGUI.info, 'iterations')
    fprintf('  Solver iterations: %d\n', resultGUI.info.iterations);
end

fprintf('  Runtime: %.2f seconds\n', meta.runtime);

fprintf('\n[COMPLETE] build_impt_plan finished.\n');
fprintf('           Returned: pln, stf, dij, resultGUI, meta\n');
fprintf('           This is a BASELINE plan (no RMS-ELO-SPAT novelty).\n\n');
safe_append_run_log(runCtx, sprintf('build_impt_plan completed. %s', meta.statusMessage), 'INFO');
safe_save_checkpoint(runCtx, 'final_summary', struct( ...
    'converged', meta.converged, ...
    'statusMessage', meta.statusMessage, ...
    'exitFlag', meta.exitFlag, ...
    'runtimeSeconds', meta.runtime, ...
    'weightsCount', safe_numel_field(resultGUI, 'w')), 1);
safe_save_run_summary(runCtx, build_impt_summary( ...
    meta, runCtx, currentStage, lastSuccessfulStage, pln, stf, dij, resultGUI, ...
    totalRays, totalBixels, 'completed', ''));

end

function [runCtx, persistEnabled] = safe_create_run_context(cfg, runLabel)
runCtx = struct();
persistEnabled = false;

if ~isfield(cfg, 'outputDir') || isempty(cfg.outputDir)
    return;
end

try
    runCtx = create_run_context(cfg, runLabel);
    persistEnabled = true;
catch ME
    fprintf('[WARNING] Run context creation failed: %s\n', ME.message);
end
end

function captureState = safe_start_console_capture(runCtx, enabled)
captureState = struct('active', false, 'logFile', '', 'startedAt', '');

if ~isstruct(runCtx) || isempty(fieldnames(runCtx))
    return;
end

try
    captureState = start_run_console_capture(runCtx, enabled);
catch ME
    safe_append_run_log(runCtx, sprintf('Console capture unavailable: %s', ME.message), 'WARN');
end
end

function safe_stop_console_capture(captureState)
if ~isstruct(captureState) || ~isfield(captureState, 'active')
    return;
end

try
    stop_run_console_capture(captureState);
catch
    % Best effort only.
end
end

function safe_append_run_log(runCtx, message, level)
if ~isstruct(runCtx) || isempty(fieldnames(runCtx))
    return;
end

try
    append_run_log(runCtx, message, level);
catch
    % Best effort only.
end
end

function safe_save_checkpoint(runCtx, stageName, checkpointData, iteration)
if ~isstruct(runCtx) || isempty(fieldnames(runCtx))
    return;
end

try
    save_checkpoint(runCtx, stageName, checkpointData, iteration);
catch ME
    safe_append_run_log(runCtx, sprintf('Checkpoint save failed for %s: %s', stageName, ME.message), 'WARN');
end
end

function safe_save_run_summary(runCtx, summaryData)
if ~isstruct(runCtx) || isempty(fieldnames(runCtx))
    return;
end

try
    save_run_summary(runCtx, summaryData);
catch ME
    safe_append_run_log(runCtx, sprintf('Summary save failed: %s', ME.message), 'WARN');
end
end

function outputPaths = build_output_paths(runCtx)
outputPaths = struct( ...
    'runId', '', ...
    'runDir', '', ...
    'logFile', '', ...
    'consoleLogFile', '', ...
    'summaryFile', '', ...
    'checkpointDir', '', ...
    'saveIntermediate', false);

if ~isstruct(runCtx) || isempty(fieldnames(runCtx))
    return;
end

outputPaths.runId = safe_struct_value(runCtx, 'runId', '');
outputPaths.runDir = safe_struct_value(runCtx, 'outputDir', '');
outputPaths.logFile = safe_struct_value(runCtx, 'logFile', '');
outputPaths.consoleLogFile = safe_struct_value(runCtx, 'consoleLogFile', '');
outputPaths.summaryFile = safe_struct_value(runCtx, 'summaryFile', '');
outputPaths.checkpointDir = safe_struct_value(runCtx, 'checkpointDir', '');
outputPaths.saveIntermediate = safe_struct_value(runCtx, 'saveIntermediate', false);
end

function summaryData = build_impt_summary(meta, runCtx, currentStage, lastSuccessfulStage, pln, stf, dij, resultGUI, totalRays, totalBixels, runStatus, errorMessage)
summaryData = struct();
summaryData.methodName = safe_struct_value(meta, 'methodName', 'IMPT baseline');
summaryData.runStatus = runStatus;
summaryData.currentStage = currentStage;
summaryData.lastSuccessfulStage = lastSuccessfulStage;
summaryData.patientFile = safe_struct_value(meta, 'patientFile', 'unknown');
summaryData.runId = safe_struct_value(meta, 'runId', '');
summaryData.outputPaths = build_output_paths(runCtx);
summaryData.statusMessage = safe_struct_value(meta, 'statusMessage', '');
summaryData.errorMessage = errorMessage;
summaryData.converged = safe_struct_value(meta, 'converged', false);
summaryData.exitFlag = safe_struct_value(meta, 'exitFlag', NaN);
summaryData.exitStatus = safe_struct_value(meta, 'exitStatus', NaN);
summaryData.exitMessage = safe_struct_value(meta, 'exitMessage', '');
summaryData.runtimeSeconds = safe_struct_value(meta, 'runtime', NaN);
summaryData.machine = safe_struct_value(meta, 'machine', '');
summaryData.bioModel = safe_struct_value(meta, 'bioModel', '');
summaryData.quantityOpt = safe_struct_value(meta, 'quantityOpt', '');
summaryData.numOfFractions = safe_struct_value(meta, 'numOfFractions', NaN);
summaryData.bixelWidth = safe_struct_value(meta, 'bixelWidth', NaN);
summaryData.doseGrid = safe_struct_value(meta, 'doseGrid', struct());
summaryData.gantryAngles = safe_struct_value(meta, 'gantryAngles', []);
summaryData.couchAngles = safe_struct_value(meta, 'couchAngles', []);
summaryData.isoCenter = safe_struct_value(meta, 'isoCenter', []);
summaryData.nConfiguredBeams = safe_configured_beam_count(pln);
summaryData.nGeneratedBeams = numel(stf);
summaryData.totalRays = totalRays;
summaryData.totalBixelsFromStf = totalBixels;
summaryData.totalBixels = safe_struct_value(meta, 'nBixels', totalBixels);
summaryData.numOfVoxels = safe_get_dij_voxels(dij);
summaryData.weightsCount = safe_numel_field(resultGUI, 'w');
summaryData.hasDoseCube = isfield(resultGUI, 'physicalDose') && ~isempty(resultGUI.physicalDose);
summaryData.objectiveValue = safe_numeric_field(resultGUI, 'fval');
summaryData.solverIterations = extract_solver_iterations(meta, resultGUI);
end

function value = safe_struct_value(s, fieldName, defaultValue)
value = defaultValue;

if ~isstruct(s) || ~isfield(s, fieldName)
    return;
end

fieldValue = s.(fieldName);
if isempty(fieldValue)
    return;
end

value = fieldValue;
end

function voxels = safe_get_dij_voxels(dij)
voxels = 0;

if isstruct(dij) && isfield(dij, 'doseGrid') && isfield(dij.doseGrid, 'numOfVoxels') ...
        && ~isempty(dij.doseGrid.numOfVoxels)
    voxels = dij.doseGrid.numOfVoxels;
end
end

function count = safe_configured_beam_count(pln)
count = 0;

if isstruct(pln) && isfield(pln, 'propStf') && isstruct(pln.propStf) ...
        && isfield(pln.propStf, 'gantryAngles') && ~isempty(pln.propStf.gantryAngles)
    count = numel(pln.propStf.gantryAngles);
end
end

function count = safe_numel_field(s, fieldName)
count = 0;

if isstruct(s) && isfield(s, fieldName) && ~isempty(s.(fieldName))
    count = numel(s.(fieldName));
end
end

function value = safe_numeric_field(s, fieldName)
value = NaN;

if isstruct(s) && isfield(s, fieldName) && isnumeric(s.(fieldName)) ...
        && isscalar(s.(fieldName)) && ~isempty(s.(fieldName))
    value = s.(fieldName);
end
end

function iterations = extract_solver_iterations(meta, resultGUI)
iterations = NaN;

if isstruct(meta) && isfield(meta, 'optimizerResultInfo') ...
        && isstruct(meta.optimizerResultInfo) ...
        && isfield(meta.optimizerResultInfo, 'iterations') ...
        && ~isempty(meta.optimizerResultInfo.iterations)
    iterations = meta.optimizerResultInfo.iterations;
    return;
end

if isstruct(resultGUI) && isfield(resultGUI, 'info') && isstruct(resultGUI.info) ...
        && isfield(resultGUI.info, 'iterations') && ~isempty(resultGUI.info.iterations)
    iterations = resultGUI.info.iterations;
end
end
