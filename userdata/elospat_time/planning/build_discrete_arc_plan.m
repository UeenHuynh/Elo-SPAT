function [pln, stf, dij] = build_discrete_arc_plan(ct, cst, cfg)
% build_discrete_arc_plan  Builds a discrete proton arc geometry and computes dij.
%
% This function creates a discrete proton arc therapy plan using many
% coplanar control points distributed around the patient. It serves as the
% geometry and dij foundation for Dense Arc and RMS-ELO-SPAT workflows.
%
% INPUTS
%   ct  - matRad CT structure
%   cst - matRad structure/objective table
%   cfg - RMS-ELO-SPAT configuration struct
%
% OUTPUTS
%   pln - matRad plan structure with configured arc geometry
%   stf - matRad steering information
%   dij - matRad dose influence matrix
%
% RESEARCH-ONLY DISCLAIMER:
%   This function is part of RMS-ELO-SPAT, a research-only module.
%   Not for clinical use. matRad is used as a research sandbox only.

% -------------------------------------------------------------------------
% Input validation
% -------------------------------------------------------------------------
narginchk(3, 3);
startTime = tic;

pln = struct();
stf = struct();
dij = struct();
totalRays = 0;
totalBixels = 0;
estimatedMemoryMB = NaN;
currentStage = 'initialization';
lastSuccessfulStage = '';
runCtx = struct();
consoleCaptureState = struct('active', false, 'logFile', '', 'startedAt', '');

if ~isstruct(ct) || ~isfield(ct, 'resolution')
    error('build_discrete_arc_plan:invalidInput', ...
        'ct must be a valid matRad CT structure with ct.resolution field.');
end

if ~iscell(cst) || size(cst, 2) < 6
    error('build_discrete_arc_plan:invalidInput', ...
        'cst must be a valid matRad cell array with at least 6 columns.');
end

if ~isstruct(cfg) || ~isfield(cfg, 'arcSpacing')
    error('build_discrete_arc_plan:invalidInput', ...
        'cfg must be a valid configuration struct with cfg.arcSpacing field.');
end

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

if ~isfield(cfg, 'numOfFractions') || isempty(cfg.numOfFractions)
    cfg.numOfFractions = 1;
end

if ~isfield(cfg, 'bixelWidth') || isempty(cfg.bixelWidth)
    cfg.bixelWidth = 5;
end

if ~isfield(cfg, 'doseGrid') || isempty(cfg.doseGrid)
    cfg.doseGrid = ct.resolution;
else
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

% -------------------------------------------------------------------------
% Validate arc spacing
% -------------------------------------------------------------------------
arcSpacing = cfg.arcSpacing;
if mod(360, arcSpacing) ~= 0
    error('build_discrete_arc_plan:invalidArcSpacing', ...
        ['Arc spacing %.2f degrees does not divide 360 degrees evenly. ', ...
         'Supported default spacings: 20, 10, 5, 4, 3, 2 degrees. ', ...
         'For custom nonuniform gantry angles, use build_impt_plan instead.'], ...
        arcSpacing);
end

[runCtx, persistEnabled] = safe_create_run_context(cfg, sprintf('discrete_arc_%gdeg', arcSpacing));
if persistEnabled
    consoleCaptureState = safe_start_console_capture(runCtx, cfg.captureConsoleOutput);
    consoleCaptureCleanup = onCleanup(@() safe_stop_console_capture(consoleCaptureState)); %#ok<NASGU>
    safe_append_run_log(runCtx, sprintf('build_discrete_arc_plan started at %.2f-degree spacing.', arcSpacing), 'INFO');
end

% -------------------------------------------------------------------------
% Step 1: Generate discrete arc geometry
% -------------------------------------------------------------------------
fprintf('\n=== build_discrete_arc_plan: Configuring discrete proton arc ===\n');
fprintf('  Arc spacing: %.2f degrees\n', arcSpacing);

gantryAngles = 0:arcSpacing:(360-arcSpacing);
nControlPoints = numel(gantryAngles);
couchAngles = zeros(1, nControlPoints);

fprintf('  Number of control points: %d\n', nControlPoints);
fprintf('  Gantry angles: [%s] degrees\n', ...
    strjoin(arrayfun(@(x) sprintf('%.1f', x), gantryAngles, 'UniformOutput', false), ', '));
fprintf('  Radiation mode: %s\n', cfg.radiationMode);
fprintf('  Machine: %s\n', cfg.machine);
fprintf('  Biological model: %s\n', cfg.bioModel);

% -------------------------------------------------------------------------
% Step 2: Configure matRad plan structure
% -------------------------------------------------------------------------
pln = struct();
pln.radiationMode = cfg.radiationMode;
pln.machine = cfg.machine;

if isfield(cfg, 'bioModel') && ~isempty(cfg.bioModel)
    pln.bioModel = cfg.bioModel;
else
    pln.bioModel = 'none';
end

if strcmp(pln.bioModel, 'constRBE')
    pln.propOpt.quantityOpt = 'RBExDose';
else
    pln.propOpt.quantityOpt = 'physicalDose';
end

pln.multScen = 'nomScen';
pln.propStf.gantryAngles = gantryAngles(:)';
pln.propStf.couchAngles = couchAngles;
pln.propStf.bixelWidth = cfg.bixelWidth;

if isfield(cfg, 'isoCenter') && ~isempty(cfg.isoCenter)
    pln.propStf.isoCenter = cfg.isoCenter;
else
    pln.propStf.isoCenter = matRad_getIsoCenter(cst, ct, 0);
end

pln.propDoseCalc.doseGrid.resolution.x = cfg.doseGrid.x;
pln.propDoseCalc.doseGrid.resolution.y = cfg.doseGrid.y;
pln.propDoseCalc.doseGrid.resolution.z = cfg.doseGrid.z;

if isfield(cfg, 'doseCalc') && isfield(cfg.doseCalc, 'engine')
    pln.propDoseCalc.engine = cfg.doseCalc.engine;
else
    pln.propDoseCalc.engine = 'HongPB';
end
pln.propDoseCalc.calcLET = 0;
pln.numOfFractions = cfg.numOfFractions;
pln.rmsPersist = build_output_paths(runCtx);

currentStage = 'plan_configured';
lastSuccessfulStage = 'plan_configured';
safe_save_run_summary(runCtx, build_arc_summary( ...
    runCtx, currentStage, lastSuccessfulStage, pln, stf, dij, cfg, arcSpacing, ...
    totalRays, totalBixels, estimatedMemoryMB, 'running', '', toc(startTime)));

% -------------------------------------------------------------------------
% Step 3: Generate steering information
% -------------------------------------------------------------------------
fprintf('\n--- Generating beam steering information (matRad_generateStf) ---\n');
currentStage = 'stf_generation';
safe_append_run_log(runCtx, 'Starting STF generation for discrete arc.', 'INFO');
safe_save_run_summary(runCtx, build_arc_summary( ...
    runCtx, currentStage, lastSuccessfulStage, pln, stf, dij, cfg, arcSpacing, ...
    totalRays, totalBixels, estimatedMemoryMB, 'running', '', toc(startTime)));

try
    stf = matRad_generateStf(ct, cst, pln);
    fprintf('[OK] Steering information generated successfully.\n');
    fprintf('     Number of beams: %d\n', numel(stf));

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

    lastSuccessfulStage = 'stf_generation';
    safe_append_run_log(runCtx, sprintf( ...
        'Discrete arc STF completed: %d beams, %d rays, %d bixels.', ...
        numel(stf), totalRays, totalBixels), 'INFO');
    safe_save_checkpoint(runCtx, 'stf_summary', struct( ...
        'nBeams', numel(stf), ...
        'totalRays', totalRays, ...
        'totalBixels', totalBixels, ...
        'gantryAngles', pln.propStf.gantryAngles, ...
        'couchAngles', pln.propStf.couchAngles), 1);
    safe_save_run_summary(runCtx, build_arc_summary( ...
        runCtx, currentStage, lastSuccessfulStage, pln, stf, dij, cfg, arcSpacing, ...
        totalRays, totalBixels, estimatedMemoryMB, 'running', '', toc(startTime)));
catch ME
    safe_append_run_log(runCtx, sprintf('STF generation failed: %s', ME.message), 'ERROR');
    safe_save_checkpoint(runCtx, 'stf_failure', struct( ...
        'errorMessage', ME.message, ...
        'runtimeSeconds', toc(startTime), ...
        'currentStage', currentStage), 1);
    safe_save_run_summary(runCtx, build_arc_summary( ...
        runCtx, currentStage, lastSuccessfulStage, pln, stf, dij, cfg, arcSpacing, ...
        totalRays, totalBixels, estimatedMemoryMB, 'failed', ME.message, toc(startTime)));
    error('build_discrete_arc_plan:stfGenerationFailed', ...
        'matRad_generateStf failed: %s', ME.message);
end

% -------------------------------------------------------------------------
% Step 4: Calculate dose influence matrix
% -------------------------------------------------------------------------
fprintf('\n--- Computing dose influence matrix (matRad_calcDoseInfluence) ---\n');
fprintf('    NOTE: This may require significant memory for dense arcs.\n');
fprintf('          If memory is insufficient, try coarser arc spacing.\n');

currentStage = 'dij_calculation';
safe_append_run_log(runCtx, 'Starting discrete arc dij calculation.', 'INFO');
safe_save_checkpoint(runCtx, 'before_dij', struct( ...
    'nBeams', numel(stf), ...
    'totalRays', totalRays, ...
    'totalBixels', totalBixels, ...
    'arcSpacing', arcSpacing, ...
    'doseGrid', cfg.doseGrid), 1);
safe_save_run_summary(runCtx, build_arc_summary( ...
    runCtx, currentStage, lastSuccessfulStage, pln, stf, dij, cfg, arcSpacing, ...
    totalRays, totalBixels, estimatedMemoryMB, 'running', '', toc(startTime)));

try
    dij = matRad_calcDoseInfluence(ct, cst, stf, pln);
    fprintf('[OK] Dose influence matrix computed successfully.\n');
    fprintf('     Total number of bixels: %d\n', dij.totalNumOfBixels);
    fprintf('     Number of voxels: %d\n', dij.doseGrid.numOfVoxels);

    if iscell(dij.physicalDose)
        nScenarios = numel(dij.physicalDose);
        fprintf('     Multi-scenario dij detected: %d scenarios\n', nScenarios);
    else
        fprintf('     Single-scenario (nominal) dij\n');
    end

    estimatedMemoryMB = estimate_dij_memory_mb(dij);
    lastSuccessfulStage = 'dij_calculation';
    safe_append_run_log(runCtx, sprintf( ...
        'Discrete arc dij completed: %d bixels, %d voxels, estimated %.1f MB.', ...
        dij.totalNumOfBixels, dij.doseGrid.numOfVoxels, estimatedMemoryMB), 'INFO');
    safe_save_checkpoint(runCtx, 'dij_summary', struct( ...
        'totalNumOfBixels', dij.totalNumOfBixels, ...
        'numOfVoxels', dij.doseGrid.numOfVoxels, ...
        'estimatedMemoryMB', estimatedMemoryMB, ...
        'isMultiScenario', iscell(dij.physicalDose)), 1);
    safe_save_run_summary(runCtx, build_arc_summary( ...
        runCtx, currentStage, lastSuccessfulStage, pln, stf, dij, cfg, arcSpacing, ...
        totalRays, totalBixels, estimatedMemoryMB, 'running', '', toc(startTime)));
catch ME
    safe_append_run_log(runCtx, sprintf('Dij calculation failed: %s', ME.message), 'ERROR');
    safe_save_checkpoint(runCtx, 'dij_failure', struct( ...
        'errorMessage', ME.message, ...
        'runtimeSeconds', toc(startTime), ...
        'nBeams', numel(stf), ...
        'totalRays', totalRays, ...
        'totalBixels', totalBixels, ...
        'arcSpacing', arcSpacing), 1);
    safe_save_run_summary(runCtx, build_arc_summary( ...
        runCtx, currentStage, lastSuccessfulStage, pln, stf, dij, cfg, arcSpacing, ...
        totalRays, totalBixels, estimatedMemoryMB, 'failed', ME.message, toc(startTime)));

    if contains(ME.message, 'memory', 'IgnoreCase', true) || ...
            contains(ME.message, 'Out of memory', 'IgnoreCase', true)
        error('build_discrete_arc_plan:memoryExceeded', ...
            ['Insufficient memory to compute dose influence matrix.\n', ...
             'Current arc spacing: %.2f degrees (%d control points)\n', ...
             'Recommendations:\n', ...
             '  1. Use coarser arc spacing (e.g., 20 instead of 10 degrees)\n', ...
             '  2. Reduce dose grid resolution in ct.resolution\n', ...
             '  3. Use a machine with more RAM\n', ...
             'Original error: %s'], ...
            arcSpacing, nControlPoints, ME.message);
    else
        error('build_discrete_arc_plan:dijCalculationFailed', ...
            'matRad_calcDoseInfluence failed: %s', ME.message);
    end
end

% -------------------------------------------------------------------------
% Step 5: Validate dij
% -------------------------------------------------------------------------
fprintf('\n--- Validating dose influence matrix ---\n');
if ~isfield(dij, 'totalNumOfBixels') || dij.totalNumOfBixels <= 0
    error('build_discrete_arc_plan:invalidDij', ...
        'Dose influence matrix is invalid: totalNumOfBixels = %d (expected > 0)', ...
        dij.totalNumOfBixels);
end

fprintf('[OK] Dose influence matrix validation passed.\n');
fprintf('     dij.totalNumOfBixels = %d (> 0)\n', dij.totalNumOfBixels);
lastSuccessfulStage = 'dij_validated';
currentStage = 'dij_validated';
safe_append_run_log(runCtx, 'Discrete arc dij validation passed.', 'INFO');
safe_save_run_summary(runCtx, build_arc_summary( ...
    runCtx, currentStage, lastSuccessfulStage, pln, stf, dij, cfg, arcSpacing, ...
    totalRays, totalBixels, estimatedMemoryMB, 'running', '', toc(startTime)));

% -------------------------------------------------------------------------
% Step 6: Report output summary
% -------------------------------------------------------------------------
fprintf('\n--- Output summary ---\n');
fprintf('  Control points: %d\n', nControlPoints);
fprintf('  Arc spacing: %.2f degrees\n', arcSpacing);
fprintf('  Total bixels: %d\n', dij.totalNumOfBixels);
fprintf('  Total rays: %d\n', totalRays);

estimatedMemoryMB = estimate_dij_memory_mb(dij);
fprintf('  Estimated dij memory: ~%.1f MB\n', estimatedMemoryMB);

fprintf('\n[COMPLETE] build_discrete_arc_plan finished successfully.\n');
fprintf('           pln, stf, dij structures returned.\n');
fprintf('           Ready for optimization (run_dense_arc_baseline or RMS-ELO-SPAT).\n\n');

currentStage = 'completed';
lastSuccessfulStage = 'completed';
safe_append_run_log(runCtx, 'build_discrete_arc_plan completed successfully.', 'INFO');
safe_save_checkpoint(runCtx, 'final_summary', struct( ...
    'nBeams', numel(stf), ...
    'totalRays', totalRays, ...
    'totalNumOfBixels', dij.totalNumOfBixels, ...
    'estimatedMemoryMB', estimatedMemoryMB, ...
    'runtimeSeconds', toc(startTime)), 1);
safe_save_run_summary(runCtx, build_arc_summary( ...
    runCtx, currentStage, lastSuccessfulStage, pln, stf, dij, cfg, arcSpacing, ...
    totalRays, totalBixels, estimatedMemoryMB, 'completed', '', toc(startTime)));

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

function summaryData = build_arc_summary(runCtx, currentStage, lastSuccessfulStage, pln, stf, dij, cfg, arcSpacing, totalRays, totalBixels, estimatedMemoryMB, runStatus, errorMessage, runtimeSeconds)
summaryData = struct();
summaryData.methodName = 'Discrete proton arc surrogate';
summaryData.runStatus = runStatus;
summaryData.currentStage = currentStage;
summaryData.lastSuccessfulStage = lastSuccessfulStage;
summaryData.patientFile = safe_struct_value(cfg, 'patientFile', 'unknown');
summaryData.arcSpacing = arcSpacing;
summaryData.runtimeSeconds = runtimeSeconds;
summaryData.errorMessage = errorMessage;
summaryData.machine = safe_struct_value(cfg, 'machine', '');
summaryData.bioModel = safe_struct_value(cfg, 'bioModel', '');
summaryData.numOfFractions = safe_struct_value(cfg, 'numOfFractions', NaN);
summaryData.bixelWidth = safe_struct_value(cfg, 'bixelWidth', NaN);
summaryData.doseGrid = safe_struct_value(cfg, 'doseGrid', struct());
summaryData.outputPaths = build_output_paths(runCtx);
summaryData.nConfiguredBeams = safe_configured_beam_count(pln);
summaryData.nGeneratedBeams = numel(stf);
summaryData.gantryAngles = safe_gantry_angles(pln);
summaryData.couchAngles = safe_couch_angles(pln);
summaryData.totalRays = totalRays;
summaryData.totalBixelsFromStf = totalBixels;
summaryData.totalBixels = safe_get_total_bixels(dij, totalBixels);
summaryData.numOfVoxels = safe_get_dij_voxels(dij);
summaryData.estimatedDijMemoryMB = estimatedMemoryMB;
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

function count = safe_configured_beam_count(pln)
count = 0;

if isstruct(pln) && isfield(pln, 'propStf') && isstruct(pln.propStf) ...
        && isfield(pln.propStf, 'gantryAngles') && ~isempty(pln.propStf.gantryAngles)
    count = numel(pln.propStf.gantryAngles);
end
end

function gantryAngles = safe_gantry_angles(pln)
gantryAngles = [];

if isstruct(pln) && isfield(pln, 'propStf') && isstruct(pln.propStf) ...
        && isfield(pln.propStf, 'gantryAngles') && ~isempty(pln.propStf.gantryAngles)
    gantryAngles = pln.propStf.gantryAngles;
end
end

function couchAngles = safe_couch_angles(pln)
couchAngles = [];

if isstruct(pln) && isfield(pln, 'propStf') && isstruct(pln.propStf) ...
        && isfield(pln.propStf, 'couchAngles') && ~isempty(pln.propStf.couchAngles)
    couchAngles = pln.propStf.couchAngles;
end
end

function totalBixels = safe_get_total_bixels(dij, fallbackValue)
totalBixels = fallbackValue;

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

function estimatedMemoryMB = estimate_dij_memory_mb(dij)
estimatedMemoryMB = NaN;

if ~isstruct(dij) || ~isfield(dij, 'doseGrid') || ~isfield(dij.doseGrid, 'numOfVoxels') ...
        || ~isfield(dij, 'totalNumOfBixels') || isempty(dij.doseGrid.numOfVoxels) ...
        || isempty(dij.totalNumOfBixels)
    return;
end

estimatedNonzeros = dij.doseGrid.numOfVoxels * dij.totalNumOfBixels * 0.01;
if isfield(dij, 'physicalDose') && iscell(dij.physicalDose)
    nScenarios = numel(dij.physicalDose);
else
    nScenarios = 1;
end

estimatedMemoryMB = (estimatedNonzeros * 8 * nScenarios) / (1024^2);
end
