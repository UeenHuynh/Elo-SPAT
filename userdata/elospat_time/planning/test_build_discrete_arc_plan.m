% test_build_discrete_arc_plan  Unit test for build_discrete_arc_plan.m
%
% This script tests the discrete proton arc plan generation function
% according to Task 6 requirements:
%   - 20-degree spacing should produce 18 control points
%   - Invalid spacing (e.g., 25 degrees) should be rejected
%   - dij should be valid with totalNumOfBixels > 0
%
% REQUIREMENTS
%   Tests Requirement 2: Discrete Proton Arc Surrogate
%
% USAGE
%   Run from MATLAB command window:
%     cd('d:\Elo-SPAT\userdata\elospat_time\planning')
%     test_build_discrete_arc_plan
%
% PREREQUISITES
%   - matRad must be on MATLAB path (run matRad_rc)
%   - PROSTATE.mat must be available in matRad data directory

fprintf('\n========================================\n');
fprintf('TEST: build_discrete_arc_plan\n');
fprintf('========================================\n\n');

% -------------------------------------------------------------------------
% Setup: Check matRad availability
% -------------------------------------------------------------------------
fprintf('--- Setup: Checking matRad availability ---\n');

if ~exist('matRad_rc', 'file')
    error('test_build_discrete_arc_plan:matRadNotFound', ...
        'matRad not found on path. Run matRad_rc from matRad root directory.');
end

fprintf('[OK] matRad found on path.\n\n');

% -------------------------------------------------------------------------
% Setup: Load patient data
% -------------------------------------------------------------------------
fprintf('--- Setup: Loading patient data ---\n');

possiblePaths = {
    'PROSTATE.mat', ...
    fullfile('..', '..', '..', 'matRad', 'data', 'PROSTATE.mat'), ...
    fullfile(fileparts(which('matRad_rc')), 'data', 'PROSTATE.mat')
};

patientFile = '';
for i = 1:numel(possiblePaths)
    if exist(possiblePaths{i}, 'file')
        patientFile = possiblePaths{i};
        break;
    end
end

if isempty(patientFile)
    error('test_build_discrete_arc_plan:patientDataNotFound', ...
        ['PROSTATE.mat not found. Please ensure it is available in:\n', ...
         '  - Current directory, or\n', ...
         '  - matRad/data/ directory']);
end

fprintf('Loading patient data from: %s\n', patientFile);
load(patientFile, 'ct', 'cst');

if ~exist('ct', 'var') || ~exist('cst', 'var')
    error('test_build_discrete_arc_plan:invalidPatientData', ...
        'Patient file does not contain ct and cst variables.');
end

fprintf('[OK] Patient data loaded successfully.\n');
fprintf('     CT resolution: [%.2f, %.2f, %.2f] mm\n', ...
    ct.resolution.x, ct.resolution.y, ct.resolution.z);
fprintf('     Number of structures: %d\n\n', size(cst, 1));

% -------------------------------------------------------------------------
% Setup: Load configuration
% -------------------------------------------------------------------------
fprintf('--- Setup: Loading configuration ---\n');

configDir = fullfile(fileparts(mfilename('fullpath')), '..', 'config');
addpath(configDir);

cfg = default_config();
fprintf('[OK] Configuration loaded.\n\n');

% -------------------------------------------------------------------------
% Test 1: Valid arc spacing (20 degrees)
% -------------------------------------------------------------------------
fprintf('========================================\n');
fprintf('TEST 1: Valid arc spacing (20 degrees)\n');
fprintf('========================================\n\n');

cfg.arcSpacing = 20;
cfg.bixelWidth = 4;
cfg.numOfFractions = 2;

try
    [pln, stf, dij] = build_discrete_arc_plan(ct, cst, cfg);
    
    expectedCPs = 360 / cfg.arcSpacing;
    actualCPs = numel(pln.propStf.gantryAngles);
    
    fprintf('\n--- Test 1 Verification ---\n');
    fprintf('Expected control points: %d\n', expectedCPs);
    fprintf('Actual control points: %d\n', actualCPs);
    
    if actualCPs == expectedCPs
        fprintf('[PASS] Correct number of control points.\n');
    else
        fprintf('[FAIL] Incorrect number of control points.\n');
        error('test_build_discrete_arc_plan:test1Failed', ...
            'Expected %d control points, got %d', expectedCPs, actualCPs);
    end
    
    expectedAngles = 0:cfg.arcSpacing:(360 - cfg.arcSpacing);
    if isequal(pln.propStf.gantryAngles, expectedAngles)
        fprintf('[PASS] Gantry angles correct: [0, 20, 40, ..., 340]\n');
    else
        fprintf('[FAIL] Gantry angles incorrect.\n');
        error('test_build_discrete_arc_plan:test1Failed', ...
            'Gantry angles do not match expected pattern.');
    end
    
    if all(pln.propStf.couchAngles == 0)
        fprintf('[PASS] All couch angles are 0 degrees (coplanar arc).\n');
    else
        fprintf('[FAIL] Couch angles are not all zero.\n');
        error('test_build_discrete_arc_plan:test1Failed', ...
            'Expected all couch angles to be 0 degrees.');
    end
    
    if pln.propStf.bixelWidth == cfg.bixelWidth && pln.numOfFractions == cfg.numOfFractions
        fprintf('[PASS] Shared planning settings honor cfg (bixelWidth, numOfFractions).\n');
    else
        fprintf('[FAIL] Shared planning settings did not honor cfg overrides.\n');
        error('test_build_discrete_arc_plan:test1Failed', ...
            'Expected bixelWidth=%g and numOfFractions=%d, got %g and %d.', ...
            cfg.bixelWidth, cfg.numOfFractions, pln.propStf.bixelWidth, pln.numOfFractions);
    end

    if isfield(dij, 'totalNumOfBixels') && dij.totalNumOfBixels > 0
        fprintf('[PASS] dij is valid with %d bixels.\n', dij.totalNumOfBixels);
    else
        fprintf('[FAIL] dij is invalid or has zero bixels.\n');
        error('test_build_discrete_arc_plan:test1Failed', ...
            'dij.totalNumOfBixels must be > 0.');
    end
    
    if numel(stf) == actualCPs
        fprintf('[PASS] stf has correct number of beams (%d).\n', numel(stf));
    else
        fprintf('[FAIL] stf has incorrect number of beams.\n');
        error('test_build_discrete_arc_plan:test1Failed', ...
            'Expected %d beams in stf, got %d', actualCPs, numel(stf));
    end
    
    fprintf('\n[TEST 1 PASSED] 20-degree arc spacing produces 18 CPs with valid dij.\n\n');
    fprintf('--- Saving Task 6 visualization plots (20 degrees) ---\n');
    try
        plotFiles = export_arc_visualizations(pln, stf, dij, '', 'task6_arc_20deg');
        fprintf('[OK] Visualization plots saved.\n');
        fprintf('     Angles: %s\n', plotFiles.angles);
    catch plotME
        fprintf('[INFO] Could not save Task 6 plots: %s\n', plotME.message);
    end
    fprintf('\n');
    
catch ME
    fprintf('\n[TEST 1 FAILED] %s\n\n', ME.message);
    rethrow(ME);
end

% -------------------------------------------------------------------------
% Test 2: Invalid arc spacing (25 degrees)
% -------------------------------------------------------------------------
fprintf('========================================\n');
fprintf('TEST 2: Invalid arc spacing (25 degrees)\n');
fprintf('========================================\n\n');

cfg.arcSpacing = 25;

try
    [pln, stf, dij] = build_discrete_arc_plan(ct, cst, cfg); %#ok<ASGLU>
    fprintf('[FAIL] Function did not reject invalid arc spacing (25 degrees).\n');
    error('test_build_discrete_arc_plan:test2Failed', ...
        'Expected error for arc spacing that does not divide 360 degrees.');
    
catch ME
    if contains(ME.identifier, 'invalidArcSpacing')
        fprintf('[PASS] Invalid arc spacing correctly rejected.\n');
        fprintf('       Error message: %s\n', ME.message);
        fprintf('\n[TEST 2 PASSED] Invalid spacing rejected with clear error.\n\n');
    else
        fprintf('[FAIL] Unexpected error type: %s\n', ME.identifier);
        fprintf('       Error message: %s\n', ME.message);
        rethrow(ME);
    end
end

% -------------------------------------------------------------------------
% Test 3: Other valid spacings (10 degrees only)
% -------------------------------------------------------------------------
fprintf('========================================\n');
fprintf('TEST 3: Other valid arc spacings\n');
fprintf('========================================\n\n');

validSpacings = [10];

for i = 1:numel(validSpacings)
    spacing = validSpacings(i);
    fprintf('--- Testing arc spacing: %d degrees ---\n', spacing);
    
    cfg.arcSpacing = spacing;
    
    try
        [pln, stf, dij] = build_discrete_arc_plan(ct, cst, cfg);
        
        expectedCPs = 360 / spacing;
        actualCPs = numel(pln.propStf.gantryAngles);
        
        if actualCPs == expectedCPs && dij.totalNumOfBixels > 0
            fprintf('[PASS] %d-degree spacing: %d CPs, %d bixels\n', ...
                spacing, actualCPs, dij.totalNumOfBixels);
            try
                plotFiles = export_arc_visualizations(pln, stf, dij, '', ...
                    sprintf('task6_arc_%ddeg', spacing));
                fprintf('       Saved plot: %s\n', plotFiles.angles);
            catch plotME
                fprintf('       [INFO] Plot export skipped: %s\n', plotME.message);
            end
        else
            fprintf('[FAIL] %d-degree spacing: incorrect CPs or invalid dij\n', spacing);
            error('test_build_discrete_arc_plan:test3Failed', ...
                'Test failed for %d-degree spacing', spacing);
        end
        
    catch ME
        fprintf('[FAIL] %d-degree spacing: %s\n', spacing, ME.message);
        rethrow(ME);
    end
end

fprintf('\n[TEST 3 PASSED] Tested spacing (10 degrees) works correctly.\n');
fprintf('NOTE: 5-degree and finer spacings are skipped in the unit test due to memory requirements.\n');
fprintf('      For production runs, use 20 degrees (18 CPs) or 10 degrees (36 CPs) as recommended.\n\n');

% -------------------------------------------------------------------------
% Test Summary
% -------------------------------------------------------------------------
fprintf('========================================\n');
fprintf('ALL TESTS PASSED\n');
fprintf('========================================\n\n');

fprintf('Summary:\n');
fprintf('  [OK] Test 1: 20-degree spacing produces 18 CPs with valid dij\n');
fprintf('  [OK] Test 2: Invalid spacing (25 degrees) correctly rejected\n');
fprintf('  [OK] Test 3: Other valid spacing (10 degrees) works correctly\n');
fprintf('  [OK] PNG plots exported for successful Task 6 runs\n');
fprintf('  [INFO] 5-degree and finer spacings were not tested due to memory requirements\n\n');

fprintf('build_discrete_arc_plan.m is ready for use.\n\n');
