function [resultGUI, meta] = run_dense_arc_baseline(dij, cst, pln)
% run_dense_arc_baseline  Dense Arc baseline optimization (Method 2).
%
% This function runs matRad fluence optimization on the full discrete arc
% dose influence matrix without any sparsification or filtering. It serves
% as the "Dense Arc" baseline method for comparison with RMS-ELO-SPAT.
%
% BASELINE-ONLY: This function does NOT implement:
%   - ELO-SPAT L₂,₁/₂ group sparsity
%   - Energy sequencing penalty
%   - BDT-aware optimization
%   - Spot/layer filtering
%   - MMU filtering
%   - Robust sparse optimization
%
% INPUTS
%   dij - matRad dose influence matrix (from build_discrete_arc_plan)
%         Must contain dij.totalNumOfBixels > 0
%   cst - matRad clinical structure table (objectives unchanged)
%   pln - matRad plan structure (from build_discrete_arc_plan)
%         Must contain pln.propStf.gantryAngles (many control points)
%
% OUTPUTS
%   resultGUI - matRad optimization result structure containing:
%               .w - optimized spot weights (N x 1 vector)
%               .info - solver status (exit flag, iterations, runtime)
%               .fval - final objective value
%               .physicalDose - dose cube (if computed)
%               .converged, .exitFlag, .statusMessage - compatibility fields
%   meta      - metadata struct for reproducibility and comparison:
%               .methodName - 'Dense Arc baseline'
%               .nControlPoints - number of control points in arc
%               .gantryAngles - array of gantry angles (degrees)
%               .arcSpacing - angular spacing (degrees, if uniform)
%               .nBixels - total number of bixels
%               .runtime - optimization runtime in seconds
%               .optimizerName - optimizer name ('IPOPT')
%               .optimizerOptions - optimizer options struct
%               .optimizerResultInfo - solver result info struct
%               .converged - boolean convergence flag (strict criteria)
%               .exitStatus - IPOPT exit flag
%               .exitFlag - compatibility alias for exitStatus
%               .exitMessage - IPOPT exit message string
%               .statusMessage - human-readable convergence status
%               .timestamp - MATLAB datetime object
%
% CONVERGENCE
%   - converged = true:  IPOPT status 0 or 1 (solved)
%   - converged = false: max iterations, infeasible, or error
%
% USAGE
%   % Step 1: Build discrete arc geometry
%   cfg = default_config();
%   cfg.arcSpacing = 20;  % 20-degree spacing → 18 control points
%   load('PROSTATE.mat', 'ct', 'cst');
%   [pln, stf, dij] = build_discrete_arc_plan(ct, cst, cfg);
%
%   % Step 2: Run Dense Arc baseline optimization
%   [resultGUI, meta] = run_dense_arc_baseline(dij, cst, pln);
%
%   if meta.converged
%       fprintf('Dense Arc baseline converged successfully.\n');
%       fprintf('Control points: %d\n', meta.nControlPoints);
%       fprintf('Total bixels: %d\n', meta.nBixels);
%       fprintf('Runtime: %.2f seconds\n', meta.runtime);
%   else
%       fprintf('Dense Arc baseline did not converge: %s\n', meta.statusMessage);
%   end
%
% ALGORITHM
%   1. Input validation (check dij, cst, pln structures)
%   2. Extract arc geometry metadata (control points, gantry angles)
%   3. Call matRad_fluenceOptimization(dij, cst, pln) with matRad objectives
%   4. Check convergence status (IPOPT exit flag)
%   5. Build metadata struct with all parameters and solver status
%   6. Return resultGUI and meta structures
%
% REQUIREMENTS
%   Implements Requirement 2 (Discrete Proton Arc Surrogate) and
%   Requirement 13 (Comparison Baseline Methods):
%   - Method 2: Dense Arc baseline
%   - Optimizes all bixels across all control points without filtering
%   - Uses matRad objectives from cst (no RMS-ELO-SPAT penalties)
%   - Serves as fair comparison baseline for RMS-ELO-SPAT methods
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
meta.methodName = 'Dense Arc baseline';
meta.timestamp = datetime('now');
startTime = tic;
resultGUI = struct();
optimizer = [];
currentStage = 'initialization';
lastSuccessfulStage = '';
runCtx = struct();
consoleCaptureState = struct('active', false, 'logFile', '', 'startedAt', '');

% -------------------------------------------------------------------------
% Input validation
% -------------------------------------------------------------------------
narginchk(3, 3);

if ~isstruct(dij) || ~isfield(dij, 'totalNumOfBixels')
    error('run_dense_arc_baseline:invalidInput', ...
        'dij must be a valid matRad dose influence matrix with dij.totalNumOfBixels field.');
end

if dij.totalNumOfBixels <= 0
    error('run_dense_arc_baseline:invalidDij', ...
        'dij.totalNumOfBixels must be > 0 (got %d).', dij.totalNumOfBixels);
end

if ~iscell(cst) || size(cst, 2) < 6
    error('run_dense_arc_baseline:invalidInput', ...
        'cst must be a valid matRad cell array with at least 6 columns.');
end

if ~isstruct(pln) || ~isfield(pln, 'propStf') || ~isfield(pln.propStf, 'gantryAngles')
    error('run_dense_arc_baseline:invalidInput', ...
        'pln must be a valid matRad plan structure with pln.propStf.gantryAngles field.');
end

% -------------------------------------------------------------------------
% Extract arc geometry metadata
% -------------------------------------------------------------------------
[runCtx, hasPersistContext] = restore_run_context_from_pln(pln);
if hasPersistContext
    consoleCaptureState = safe_start_console_capture(runCtx, true);
    consoleCaptureCleanup = onCleanup(@() safe_stop_console_capture(consoleCaptureState)); %#ok<NASGU>
    safe_append_run_log(runCtx, 'run_dense_arc_baseline attached to existing arc run context.', 'INFO');
    meta.runId = safe_struct_value(runCtx, 'runId', '');
else
    meta.runId = '';
end
meta.outputPaths = build_output_paths(runCtx);

gantryAngles = pln.propStf.gantryAngles;
nControlPoints = numel(gantryAngles);

if nControlPoints < 2
    error('run_dense_arc_baseline:invalidArc', ...
        'Arc must have at least 2 control points (got %d).', nControlPoints);
end

% Detect arc spacing (if uniform)
if nControlPoints >= 2
    spacings = diff(gantryAngles);
    if all(abs(spacings - spacings(1)) < 1e-6)
        arcSpacing = spacings(1);
    else
        arcSpacing = NaN;  % Non-uniform spacing
    end
else
    arcSpacing = NaN;
end

% Store metadata
meta.nControlPoints = nControlPoints;
meta.gantryAngles = gantryAngles;
meta.arcSpacing = arcSpacing;
meta.nBixels = dij.totalNumOfBixels;
meta.runStatus = 'running';

% -------------------------------------------------------------------------
% Console output: method header
% -------------------------------------------------------------------------
fprintf('\n=== run_dense_arc_baseline: Dense Arc Baseline (Method 2) ===\n');
fprintf('  Control points: %d\n', nControlPoints);
fprintf('  Gantry angles: [%s] degrees\n', ...
    strjoin(arrayfun(@(x) sprintf('%.1f', x), gantryAngles(1:min(5, end)), 'UniformOutput', false), ', '));
if nControlPoints > 5
    fprintf('                 ... (%d total angles)\n', nControlPoints);
end

if ~isnan(arcSpacing)
    fprintf('  Arc spacing: %.1f degrees (uniform)\n', arcSpacing);
else
    fprintf('  Arc spacing: non-uniform\n');
end

fprintf('  Total bixels: %d\n', dij.totalNumOfBixels);
fprintf('  Optimizer: IPOPT (matRad defaults, no RMS-ELO-SPAT penalties)\n');
currentStage = 'plan_received';
lastSuccessfulStage = 'plan_received';
safe_save_run_summary(runCtx, build_dense_arc_summary( ...
    meta, runCtx, currentStage, lastSuccessfulStage, pln, dij, resultGUI, ...
    'running', '', toc(startTime)));

% Check for multi-scenario dij
if iscell(dij.physicalDose)
    nScenarios = numel(dij.physicalDose);
    fprintf('  Multi-scenario dij detected: %d scenarios\n', nScenarios);
else
    fprintf('  Single-scenario (nominal) dij\n');
end

% -------------------------------------------------------------------------
% Run fluence optimization (matRad standard workflow)
% -------------------------------------------------------------------------
fprintf('\n--- Running fluence optimization (matRad_fluenceOptimization) ---\n');
fprintf('    NOTE: Using matRad objectives from cst (no RMS-ELO-SPAT penalties)\n');
fprintf('    NOTE: Optimizing ALL bixels across ALL control points (no filtering)\n');
currentStage = 'fluence_optimization';
safe_append_run_log(runCtx, sprintf( ...
    'Dense Arc optimization started for %d control points and %d bixels.', ...
    nControlPoints, dij.totalNumOfBixels), 'INFO');
safe_save_checkpoint(runCtx, 'dense_arc_before_optimization', struct( ...
    'nControlPoints', nControlPoints, ...
    'nBixels', dij.totalNumOfBixels, ...
    'arcSpacing', arcSpacing), 1);
safe_save_run_summary(runCtx, build_dense_arc_summary( ...
    meta, runCtx, currentStage, lastSuccessfulStage, pln, dij, resultGUI, ...
    'running', '', toc(startTime)));

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
    safe_append_run_log(runCtx, 'Dense Arc fluence optimization completed.', 'INFO');
    safe_save_checkpoint(runCtx, 'dense_arc_optimization_summary', struct( ...
        'weightsCount', safe_numel_field(resultGUI, 'w'), ...
        'objectiveValue', safe_numeric_field(resultGUI, 'fval'), ...
        'hasDoseCube', isfield(resultGUI, 'physicalDose') && ~isempty(resultGUI.physicalDose)), 1);
     
catch ME
    meta.converged = false;
    meta.exitFlag = NaN;
    meta.exitStatus = NaN;
    meta.exitMsg = ME.message;
    meta.exitMessage = ME.message;
    meta.statusMessage = sprintf('Optimization failed: %s', ME.message);
    meta.errorMessage = ME.message;
    meta.runtime = toc(startTime);
    meta.runStatus = 'failed';
    
    fprintf('[FAILED] matRad_fluenceOptimization error: %s\n', ME.message);
    safe_append_run_log(runCtx, meta.statusMessage, 'ERROR');
    safe_save_checkpoint(runCtx, 'dense_arc_optimization_failure', struct( ...
        'errorMessage', ME.message, ...
        'runtimeSeconds', meta.runtime, ...
        'nControlPoints', nControlPoints, ...
        'nBixels', dij.totalNumOfBixels), 1);
    
    % Return partial outputs with metadata
    resultGUI = struct( ...
        'converged', false, ...
        'exitFlag', NaN, ...
        'exitStatus', NaN, ...
        'exitMsg', ME.message, ...
        'exitMessage', ME.message, ...
        'statusMessage', meta.statusMessage);
    safe_save_run_summary(runCtx, build_dense_arc_summary( ...
        meta, runCtx, currentStage, lastSuccessfulStage, pln, dij, resultGUI, ...
        'failed', ME.message, toc(startTime)));
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
meta.optimizerName = 'IPOPT';
meta.runStatus = ternary_status(converged, status);

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
    
    % Verify weights match dij dimensions
    if numel(resultGUI.w) == dij.totalNumOfBixels
        fprintf('  Weights dimension check: PASS (matches dij.totalNumOfBixels)\n');
    else
        fprintf('  [WARNING] Weights dimension mismatch: %d weights, %d bixels\n', ...
            numel(resultGUI.w), dij.totalNumOfBixels);
    end
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

fprintf('\n[COMPLETE] run_dense_arc_baseline finished.\n');
fprintf('           Returned: resultGUI, meta\n');
fprintf('           This is a BASELINE method (no RMS-ELO-SPAT novelty).\n');
fprintf('           All %d bixels across %d control points were optimized.\n', ...
    meta.nBixels, meta.nControlPoints);
fprintf('           Use this result as Method 2 comparison baseline.\n\n');
safe_append_run_log(runCtx, sprintf('Dense Arc baseline finished. %s', meta.statusMessage), 'INFO');
safe_save_checkpoint(runCtx, 'dense_arc_final_summary', struct( ...
    'converged', meta.converged, ...
    'exitFlag', meta.exitFlag, ...
    'statusMessage', meta.statusMessage, ...
    'runtimeSeconds', meta.runtime, ...
    'weightsCount', safe_numel_field(resultGUI, 'w')), 1);
safe_save_run_summary(runCtx, build_dense_arc_summary( ...
    meta, runCtx, currentStage, lastSuccessfulStage, pln, dij, resultGUI, ...
    meta.runStatus, '', toc(startTime)));

end

function [runCtx, hasPersistContext] = restore_run_context_from_pln(pln)
runCtx = struct();
hasPersistContext = false;

if ~isstruct(pln) || ~isfield(pln, 'rmsPersist') || ~isstruct(pln.rmsPersist)
    return;
end

persist = pln.rmsPersist;
runDir = safe_struct_value(persist, 'runDir', '');
logFile = safe_struct_value(persist, 'logFile', '');
summaryFile = safe_struct_value(persist, 'summaryFile', '');
checkpointDir = safe_struct_value(persist, 'checkpointDir', '');

if isempty(runDir) || isempty(logFile) || isempty(summaryFile) || isempty(checkpointDir)
    return;
end

runCtx = struct();
runCtx.runId = safe_struct_value(persist, 'runId', '');
runCtx.label = 'dense_arc_baseline';
runCtx.outputDir = runDir;
runCtx.checkpointDir = checkpointDir;
runCtx.logFile = logFile;
runCtx.consoleLogFile = safe_struct_value(persist, 'consoleLogFile', '');
runCtx.summaryFile = summaryFile;
runCtx.configSnapshotFile = fullfile(runDir, 'config_snapshot.mat');
runCtx.saveIntermediate = logical(safe_struct_value(persist, 'saveIntermediate', true));
runCtx.createdAt = datestr(now, 'yyyy-mm-dd HH:MM:SS.FFF');
hasPersistContext = true;
end

function captureState = safe_start_console_capture(runCtx, enabled)
captureState = struct('active', false, 'logFile', '', 'startedAt', '');

if ~enabled || ~isstruct(runCtx) || isempty(fieldnames(runCtx))
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

function summaryData = build_dense_arc_summary(meta, runCtx, currentStage, lastSuccessfulStage, pln, dij, resultGUI, runStatus, errorMessage, runtimeSeconds)
summaryData = struct();
summaryData.methodName = safe_struct_value(meta, 'methodName', 'Dense Arc baseline');
summaryData.runStatus = runStatus;
summaryData.currentStage = currentStage;
summaryData.lastSuccessfulStage = lastSuccessfulStage;
summaryData.runId = safe_struct_value(meta, 'runId', '');
summaryData.outputPaths = build_output_paths(runCtx);
summaryData.errorMessage = errorMessage;
summaryData.statusMessage = safe_struct_value(meta, 'statusMessage', '');
summaryData.runtimeSeconds = runtimeSeconds;
summaryData.nControlPoints = safe_struct_value(meta, 'nControlPoints', 0);
summaryData.gantryAngles = safe_struct_value(meta, 'gantryAngles', []);
summaryData.arcSpacing = safe_struct_value(meta, 'arcSpacing', NaN);
summaryData.nBixels = safe_struct_value(meta, 'nBixels', safe_get_total_bixels(dij));
summaryData.numOfVoxels = safe_get_dij_voxels(dij);
summaryData.converged = safe_struct_value(meta, 'converged', false);
summaryData.exitFlag = safe_struct_value(meta, 'exitFlag', NaN);
summaryData.exitStatus = safe_struct_value(meta, 'exitStatus', NaN);
summaryData.exitMessage = safe_struct_value(meta, 'exitMessage', '');
summaryData.weightsCount = safe_numel_field(resultGUI, 'w');
summaryData.objectiveValue = safe_numeric_field(resultGUI, 'fval');
summaryData.hasDoseCube = isfield(resultGUI, 'physicalDose') && ~isempty(resultGUI.physicalDose);
summaryData.optimizerName = safe_struct_value(meta, 'optimizerName', 'IPOPT');
summaryData.numConfiguredBeams = safe_configured_beam_count(pln);
summaryData.acceptableButNotConverged = (~summaryData.converged) && isequal(summaryData.exitStatus, -1);
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

function totalBixels = safe_get_total_bixels(dij)
totalBixels = 0;

if isstruct(dij) && isfield(dij, 'totalNumOfBixels') && ~isempty(dij.totalNumOfBixels)
    totalBixels = dij.totalNumOfBixels;
end
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

function runStatus = ternary_status(converged, status)
if converged
    runStatus = 'completed';
elseif isequal(status, -1)
    runStatus = 'completed_nonconverged';
else
    runStatus = 'failed';
end
end
