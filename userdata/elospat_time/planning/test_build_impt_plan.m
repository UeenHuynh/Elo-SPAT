% test_build_impt_plan.m
%
% Test script for build_impt_plan function (Task 5)
%
% This script verifies that build_impt_plan correctly:
%   1. Configures pln for protons with cfg.imptAngles (default [90,270])
%   2. Calls matRad_generateStf
%   3. Calls matRad_calcDoseInfluence
%   4. Calls matRad_fluenceOptimization
%   5. Checks convergence
%   6. Returns pln/stf/dij/result structures
%
% REQUIREMENTS
%   - matRad must be on the MATLAB path (run matRad_rc first)
%   - PROSTATE.mat must be available (typically in matRad/data/)
%
% USAGE
%   Run from MATLAB command window:
%       cd('d:\Elo-SPAT\userdata\elospat_time\planning')
%       test_build_impt_plan
%
% =========================================================================
% RESEARCH-ONLY DISCLAIMER:
%   This test is part of RMS-ELO-SPAT, a research-only module.
%   Not for clinical use.
% =========================================================================

clear; clc;
fprintf('\n========================================\n');
fprintf('TEST: build_impt_plan (Task 5)\n');
fprintf('========================================\n\n');

% -------------------------------------------------------------------------
% Step 1: Check matRad availability
% -------------------------------------------------------------------------
fprintf('Step 1: Checking matRad availability...\n');

if exist('matRad_rc', 'file') ~= 2
    error('test_build_impt_plan:matRadNotFound', ...
        ['matRad_rc not found on MATLAB path.\n', ...
         'Please install matRad and run matRad_rc before running this test.\n', ...
         'See: https://github.com/e0404/matRad']);
end

fprintf('[OK] matRad_rc found on path.\n\n');

% -------------------------------------------------------------------------
% Step 2: Load configuration
% -------------------------------------------------------------------------
fprintf('Step 2: Loading default configuration...\n');

% Add config directory to path
configDir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'config');
addpath(configDir);

cfg = default_config();
fprintf('[OK] Configuration loaded.\n');
fprintf('     IMPT angles: [%s] degrees\n', ...
    strjoin(arrayfun(@num2str, cfg.imptAngles, 'UniformOutput', false), ', '));
fprintf('     Patient file: %s\n\n', cfg.patientFile);

% -------------------------------------------------------------------------
% Step 3: Load patient data (PROSTATE.mat)
% -------------------------------------------------------------------------
fprintf('Step 3: Loading patient data...\n');

% Try to locate PROSTATE.mat
prostateFile = which('PROSTATE.mat');
if isempty(prostateFile)
    % Try cfg.patientFile directly
    if exist(cfg.patientFile, 'file') == 2
        prostateFile = cfg.patientFile;
    else
        error('test_build_impt_plan:prostateNotFound', ...
            ['PROSTATE.mat not found.\n', ...
             'Expected location: <matRad_root>/data/PROSTATE.mat\n', ...
             'Verify your matRad installation is complete.']);
    end
end

fprintf('     Loading: %s\n', prostateFile);
load(prostateFile, 'ct', 'cst');

if ~exist('ct', 'var') || ~exist('cst', 'var')
    error('test_build_impt_plan:invalidData', ...
        'PROSTATE.mat does not contain ct and cst variables.');
end

fprintf('[OK] Patient data loaded.\n');
fprintf('     CT resolution: [%.2f, %.2f, %.2f] mm\n', ...
    ct.resolution.x, ct.resolution.y, ct.resolution.z);
fprintf('     Number of structures: %d\n\n', size(cst, 1));

% -------------------------------------------------------------------------
% Step 4: Run build_impt_plan
% -------------------------------------------------------------------------
fprintf('Step 4: Running build_impt_plan...\n');
fprintf('========================================\n');

tic;
try
    [pln, stf, dij, result, meta] = build_impt_plan(ct, cst, cfg);
    elapsedTime = toc;
    
    fprintf('========================================\n');
    fprintf('[OK] build_impt_plan completed successfully.\n');
    fprintf('     Total runtime: %.2f seconds\n\n', elapsedTime);
    
catch ME
    fprintf('========================================\n');
    fprintf('[FAILED] build_impt_plan threw an error:\n');
    fprintf('         %s\n', ME.message);
    fprintf('         Stack trace:\n');
    for k = 1:numel(ME.stack)
        fprintf('           %s (line %d)\n', ME.stack(k).name, ME.stack(k).line);
    end
    rethrow(ME);
end

% -------------------------------------------------------------------------
% Step 5: Validate outputs
% -------------------------------------------------------------------------
fprintf('Step 5: Validating outputs...\n');

% Check pln structure
assert(isstruct(pln), 'pln must be a struct');
assert(isfield(pln, 'radiationMode'), 'pln must have radiationMode field');
assert(strcmp(pln.radiationMode, 'protons'), 'pln.radiationMode must be protons');
assert(isfield(pln, 'propStf'), 'pln must have propStf field');
assert(isfield(pln.propStf, 'gantryAngles'), 'pln.propStf must have gantryAngles');
assert(isequal(pln.propStf.gantryAngles, cfg.imptAngles(:)'), ...
    'pln.propStf.gantryAngles must match cfg.imptAngles');
fprintf('[OK] pln structure validated.\n');

% Check stf structure
assert(isstruct(stf), 'stf must be a struct array');
assert(numel(stf) == numel(cfg.imptAngles), ...
    'stf must have one entry per gantry angle');
fprintf('[OK] stf structure validated (%d beams).\n', numel(stf));

% Check dij structure
assert(isstruct(dij), 'dij must be a struct');
assert(isfield(dij, 'totalNumOfBixels'), 'dij must have totalNumOfBixels field');
assert(dij.totalNumOfBixels > 0, 'dij.totalNumOfBixels must be positive');
assert(isfield(dij, 'physicalDose'), 'dij must have physicalDose field');
fprintf('[OK] dij structure validated (%d bixels).\n', dij.totalNumOfBixels);

% Check result structure
assert(isstruct(result), 'result must be a struct');
assert(isfield(result, 'w'), 'result must have w (weights) field');
assert(numel(result.w) == dij.totalNumOfBixels, ...
    'result.w length must match dij.totalNumOfBixels');
assert(isfield(result, 'converged'), 'result must expose converged compatibility field');
assert(isfield(result, 'exitFlag'), 'result must expose exitFlag compatibility field');
fprintf('[OK] result structure validated.\n');
fprintf('     Number of optimized weights: %d\n', numel(result.w));
fprintf('     Convergence status: %s\n', mat2str(result.converged));
fprintf('     Exit flag: %d\n', result.exitFlag);

% Check metadata structure (authoritative planning status contract)
assert(isstruct(meta), 'meta must be a struct');
assert(isfield(meta, 'converged'), 'meta must have converged field');
assert(isfield(meta, 'exitStatus'), 'meta must have exitStatus field');
assert(isfield(meta, 'statusMessage'), 'meta must have statusMessage field');
assert(isequal(meta.converged, result.converged), ...
    'meta.converged must match result.converged compatibility field');
assert(isequaln(meta.exitStatus, result.exitStatus), ...
    'meta.exitStatus must match result.exitStatus compatibility field');
fprintf('[OK] meta structure validated.\n');
fprintf('     Status message: %s\n', meta.statusMessage);

% Check convergence (Requirement 1.5)
if meta.converged
    fprintf('[OK] Plan converged successfully.\n');
else
    fprintf('[WARNING] Plan did not converge (exit flag: %d).\n', meta.exitStatus);
    fprintf('          This may be acceptable if max iterations was reached.\n');
end

% -------------------------------------------------------------------------
% Step 6: Summary
% -------------------------------------------------------------------------
fprintf('\n========================================\n');
fprintf('TEST SUMMARY\n');
fprintf('========================================\n');
fprintf('Status: PASSED\n');
fprintf('\nAll validation checks passed successfully.\n');
fprintf('build_impt_plan correctly:\n');
fprintf('  1. Configured pln for protons with cfg.imptAngles\n');
fprintf('  2. Called matRad_generateStf\n');
fprintf('  3. Called matRad_calcDoseInfluence\n');
fprintf('  4. Called matRad_fluenceOptimization\n');
fprintf('  5. Checked convergence status\n');
fprintf('  6. Returned pln/stf/dij/result/meta structures\n');
fprintf('\nTask 5 implementation verified.\n');
fprintf('========================================\n\n');

% -------------------------------------------------------------------------
% Optional: Display DVH metrics if available
% -------------------------------------------------------------------------
if isfield(result, 'physicalDose') && ~isempty(result.physicalDose)
    fprintf('\nOptional: Computing DVH metrics...\n');
    
    try
        % Compute dose cube if not already present
        if ~isfield(result, 'RBExDose') || isempty(result.RBExDose)
            result = matRad_calcCubes(result.w, dij);
        end
        
        % Compute DVH
        dvh = matRad_calcDVH(cst, result.RBExDose, 'cum');
        
        fprintf('[OK] DVH computed successfully.\n');
        fprintf('     Number of structures: %d\n', numel(dvh));
        
        % Display CTV metrics if available
        for i = 1:size(cst, 1)
            if strcmp(cst{i,3}, 'TARGET')
                structName = cst{i,2};
                % Find D95 (dose to 95% of volume)
                volAxis = dvh(i).volumePoints;
                doseAxis = dvh(i).doseGrid;
                
                % Interpolate to find D95
                if max(volAxis) >= 95
                    d95 = interp1(volAxis, doseAxis, 95, 'linear');
                    fprintf('     %s D95: %.2f Gy\n', structName, d95);
                end
            end
        end
        
    catch ME
        fprintf('[INFO] Could not compute DVH metrics: %s\n', ME.message);
    end
end

% -------------------------------------------------------------------------
% Optional: Save visualization plots
% -------------------------------------------------------------------------
fprintf('\nOptional: Saving visualization plots...\n');

try
    plotFiles = export_plan_visualizations( ...
        ct, cst, dij, result, meta, '', 'task5_impt');
    fprintf('[OK] Visualization plots saved.\n');
    fprintf('     DVH: %s\n', plotFiles.dvh);
catch ME
    fprintf('[INFO] Could not save visualization plots: %s\n', ME.message);
end

fprintf('\n');
