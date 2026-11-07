function validate_config(cfg, machine)
% validate_config  Validates all RMS-ELO-SPAT configuration fields.
%
% Checks that all required fields exist, numeric values fall within valid
% ranges, string fields match allowed values, the patient file is
% resolvable, and the arc spacing divides 360 exactly.
%
% Throws a descriptive error (with the offending field name) on the first
% invalid input encountered.
%
% INPUTS
%   cfg     - configuration struct, as returned by default_config()
%   machine - machine timing model struct, as returned by
%             default_machine_model().  Pass [] to skip machine validation.
%
% USAGE
%   cfg     = default_config();
%   machine = default_machine_model();
%   validate_config(cfg, machine);   % throws on any invalid field
%
% =========================================================================
% RESEARCH-ONLY DISCLAIMER:
%   This validator is part of RMS-ELO-SPAT, a research-only module.
%   Not for clinical use.
% =========================================================================

narginchk(1, 2);
if nargin < 2
    machine = [];
end

% -------------------------------------------------------------------------
% 1. Top-level required fields and their types
% -------------------------------------------------------------------------
validateField_string(cfg, 'mode',          {'wrapper','integrated'});
validateField_string(cfg, 'patientFile',   {});          % any non-empty string
validateField_string(cfg, 'radiationMode', {'protons'});
validateField_string(cfg, 'machine',       {});          % any non-empty string
validateField_string(cfg, 'bioModel',      {'none','constRBE','MCN','HEL'});
validateField_positive_integer(cfg, 'numOfFractions');
validateField_positive_scalar(cfg, 'bixelWidth');

requireField(cfg, 'doseGrid');
if ~isstruct(cfg.doseGrid)
    error('validate_config:invalidField', ...
        'cfg.doseGrid must be a struct with fields x, y, z.');
end
validateDoseGridField(cfg.doseGrid, 'x');
validateDoseGridField(cfg.doseGrid, 'y');
validateDoseGridField(cfg.doseGrid, 'z');

requireField(cfg, 'isoCenter');
if ~isempty(cfg.isoCenter)
    if ~isnumeric(cfg.isoCenter) || numel(cfg.isoCenter) ~= 3 || any(~isfinite(cfg.isoCenter))
        error('validate_config:invalidField', ...
            'cfg.isoCenter must be empty or a finite 3-element numeric vector [x y z].');
    end
end

% -------------------------------------------------------------------------
% 2. Patient file must exist on disk or be resolvable via MATLAB path
% -------------------------------------------------------------------------
if ~ischar(cfg.patientFile) || isempty(strtrim(cfg.patientFile))
    error('validate_config:invalidField', ...
        'cfg.patientFile must be a non-empty string.');
end

% Accept absolute paths, relative paths, or filenames on MATLAB path
fileFound = exist(cfg.patientFile, 'file') == 2;
if ~fileFound
    error('validate_config:fileNotFound', ...
        ['cfg.patientFile not found: ''%s''. ', ...
         'Ensure the patient .mat file exists on disk or is on the MATLAB path.'], ...
        cfg.patientFile);
end

% -------------------------------------------------------------------------
% 3. Geometry — IMPT angles
% -------------------------------------------------------------------------
requireField(cfg, 'imptAngles');
if ~isnumeric(cfg.imptAngles) || ~isvector(cfg.imptAngles)
    error('validate_config:invalidField', ...
        'cfg.imptAngles must be a numeric vector.');
end
nAngles = numel(cfg.imptAngles);
if nAngles < 2 || nAngles > 5
    error('validate_config:outOfRange', ...
        'cfg.imptAngles must contain 2 to 5 angles (got %d).', nAngles);
end
if any(cfg.imptAngles < 0) || any(cfg.imptAngles >= 360)
    error('validate_config:outOfRange', ...
        'cfg.imptAngles values must be in [0, 360).');
end

% -------------------------------------------------------------------------
% 4. Geometry — arc spacing
% -------------------------------------------------------------------------
requireField(cfg, 'arcSpacing');
validateField_scalar_numeric(cfg, 'arcSpacing', 1, 180);

if mod(360, cfg.arcSpacing) ~= 0
    error('validate_config:invalidArcSpacing', ...
        ['cfg.arcSpacing = %g does not divide 360 exactly. ', ...
         'Use a value such that 360/arcSpacing is an integer ', ...
         '(e.g. 1,2,3,4,5,6,8,9,10,12,15,18,20,24,30,36,40,45,60,72,90,120,180).'], ...
        cfg.arcSpacing);
end

requireField(cfg, 'couchAngles');
if ~isnumeric(cfg.couchAngles)
    error('validate_config:invalidField', ...
        'cfg.couchAngles must be numeric (scalar or vector).');
end

% -------------------------------------------------------------------------
% 5. Penalty weights — all must be non-negative scalars
% -------------------------------------------------------------------------
for fname = {'lambdaEL','lambdaB','lambdaES','lambdaS','lambdaT'}
    requireField(cfg, fname{1});
    validateField_scalar_numeric(cfg, fname{1}, 0, Inf);
end

% -------------------------------------------------------------------------
% 6. Smoothing parameters — valid ranges from spec
% -------------------------------------------------------------------------
validateField_scalar_numeric(cfg, 'epsEL', 1e-10, 1e-6);
validateField_scalar_numeric(cfg, 'epsB',  1e-10, 1e-6);
validateField_scalar_numeric(cfg, 'epsZ',  1e-8,  1e-4);
validateField_scalar_numeric(cfg, 'epsS',  1e-5,  1e-1);

% -------------------------------------------------------------------------
% 7. Thresholds
% -------------------------------------------------------------------------
validateField_scalar_numeric(cfg, 'layerThreshold', 1e-6,  1e-4);
validateField_scalar_numeric(cfg, 'spotThreshold',  1e-8,  1e-4);
validateField_scalar_numeric(cfg, 'mmuThreshold',   0,     Inf);

% -------------------------------------------------------------------------
% 8. Robustness sub-struct
% -------------------------------------------------------------------------
requireSubField(cfg, 'robust', 'enabled');
if ~(islogical(cfg.robust.enabled) || ...
        (isnumeric(cfg.robust.enabled) && isscalar(cfg.robust.enabled) && any(cfg.robust.enabled == [0 1])))
    error('validate_config:invalidField', ...
        'cfg.robust.enabled must be logical true/false (or numeric 0/1).');
end

requireSubField(cfg, 'robust', 'mode');
if ~ischar(cfg.robust.mode) || ~ismember(cfg.robust.mode, {'COWC','OWC'})
    error('validate_config:invalidField', ...
        'cfg.robust.mode must be ''COWC'' or ''OWC'' (got ''%s'').', ...
        cfg.robust.mode);
end

requireSubField(cfg, 'robust', 'setupMM');
if ~isnumeric(cfg.robust.setupMM) || numel(cfg.robust.setupMM) ~= 3 ...
        || any(cfg.robust.setupMM < 0)
    error('validate_config:invalidField', ...
        'cfg.robust.setupMM must be a non-negative 3-element numeric vector [x y z] in mm.');
end

requireSubField(cfg, 'robust', 'rangePercent');
if ~isnumeric(cfg.robust.rangePercent) || ~isscalar(cfg.robust.rangePercent) ...
        || cfg.robust.rangePercent < 0
    error('validate_config:invalidField', ...
        'cfg.robust.rangePercent must be a non-negative scalar (got %g).', ...
        cfg.robust.rangePercent);
end

% -------------------------------------------------------------------------
% 9. Wrapper loop settings
% -------------------------------------------------------------------------
requireSubField(cfg, 'wrapper', 'maxLayerIterations');
requireSubField(cfg, 'wrapper', 'maxSpotIterations');
requireSubField(cfg, 'wrapper', 'relTol');

validateSubField_positive_integer(cfg, 'wrapper', 'maxLayerIterations');
validateSubField_positive_integer(cfg, 'wrapper', 'maxSpotIterations');
validateSubField_scalar_numeric_sub(cfg, 'wrapper', 'relTol', 1e-6, 1e-2);

% -------------------------------------------------------------------------
% 9b. Integrated MATLAB optimizer
% -------------------------------------------------------------------------
requireSubField(cfg, 'integrated', 'componentStage');
if ~ischar(cfg.integrated.componentStage) || ...
        ~ismember(cfg.integrated.componentStage, {'layer','spot','full'})
    error('validate_config:invalidField', ...
        'cfg.integrated.componentStage must be layer, spot, or full.');
end
for fname = {'maxIterations','maxBacktracking','printEvery', ...
        'minimumIterations','stagnationWindow','spotMMUpdatePeriod'}
    validateSubField_positive_integer(cfg, 'integrated', fname{1});
end
for fname = {'initialStep','minimumStep','maximumStep', ...
        'relativeStepTolerance','relativeObjectiveTolerance', ...
        'stationarityTolerance'}
    requireSubField(cfg, 'integrated', fname{1});
    value = cfg.integrated.(fname{1});
    if ~isnumeric(value) || ~isscalar(value) || ~isfinite(value) || value <= 0
        error('validate_config:invalidField', ...
            'cfg.integrated.%s must be positive and finite.', fname{1});
    end
end
for fname = {'stagnationRelativeStepTolerance', ...
        'stagnationRelativeObjectiveTolerance'}
    requireSubField(cfg, 'integrated', fname{1});
    value = cfg.integrated.(fname{1});
    if ~isnumeric(value) || ~isscalar(value) || ~isfinite(value) || value < 0
        error('validate_config:invalidField', ...
            'cfg.integrated.%s must be non-negative and finite.', fname{1});
    end
end
if cfg.integrated.minimumIterations > cfg.integrated.maxIterations
    error('validate_config:invalidField', ...
        'minimumIterations must not exceed maxIterations.');
end
for fname = {'backtrackFactor','backtrackGrowth'}
    requireSubField(cfg, 'integrated', fname{1});
    value = cfg.integrated.(fname{1});
    if ~isnumeric(value) || ~isscalar(value) || value <= 0 || value >= 1
        error('validate_config:invalidField', ...
            'cfg.integrated.%s must lie strictly between 0 and 1.', fname{1});
    end
end
for fname = {'sufficientDecreaseTolerance','constraintPenalty'}
    requireSubField(cfg, 'integrated', fname{1});
    value = cfg.integrated.(fname{1});
    if ~isnumeric(value) || ~isscalar(value) || ~isfinite(value) || value < 0
        error('validate_config:invalidField', ...
            'cfg.integrated.%s must be non-negative and finite.', fname{1});
    end
end
requireSubField(cfg, 'integrated', 'verbose');
if ~(islogical(cfg.integrated.verbose) || ...
        (isnumeric(cfg.integrated.verbose) && isscalar(cfg.integrated.verbose) && ...
        any(cfg.integrated.verbose == [0 1])))
    error('validate_config:invalidField', ...
        'cfg.integrated.verbose must be logical.');
end
requireSubField(cfg, 'integrated', 'enforceFinalMinimumWeight');
if ~(islogical(cfg.integrated.enforceFinalMinimumWeight) || ...
        (isnumeric(cfg.integrated.enforceFinalMinimumWeight) && ...
        isscalar(cfg.integrated.enforceFinalMinimumWeight) && ...
        any(cfg.integrated.enforceFinalMinimumWeight == [0 1])))
    error('validate_config:invalidField', ...
        'cfg.integrated.enforceFinalMinimumWeight must be logical.');
end
requireSubField(cfg, 'integrated', 'stationarityStoppingEnabled');
if ~(islogical(cfg.integrated.stationarityStoppingEnabled) || ...
        (isnumeric(cfg.integrated.stationarityStoppingEnabled) && ...
        isscalar(cfg.integrated.stationarityStoppingEnabled) && ...
        any(cfg.integrated.stationarityStoppingEnabled == [0 1])))
    error('validate_config:invalidIntegratedStationarityStopping', ...
        'cfg.integrated.stationarityStoppingEnabled must be logical.');
end
if cfg.integrated.minimumStep > cfg.integrated.initialStep || ...
        cfg.integrated.initialStep > cfg.integrated.maximumStep
    error('validate_config:invalidField', ...
        'Integrated step sizes must satisfy minimum <= initial <= maximum.');
end

% -------------------------------------------------------------------------
% 10. Evaluation tolerances
% -------------------------------------------------------------------------
requireSubField(cfg, 'evaluation', 'ctvD95TolerancePercent');
requireSubField(cfg, 'evaluation', 'robustD95TolerancePercent');
validateSubField_scalar_numeric_sub(cfg, 'evaluation', 'ctvD95TolerancePercent',    0, 100);
validateSubField_scalar_numeric_sub(cfg, 'evaluation', 'robustD95TolerancePercent', 0, 100);

% -------------------------------------------------------------------------
% 11. Output settings
% -------------------------------------------------------------------------
requireField(cfg, 'outputDir');
if ~ischar(cfg.outputDir) || isempty(strtrim(cfg.outputDir))
    error('validate_config:invalidField', ...
        'cfg.outputDir must be a non-empty string path.');
end

requireField(cfg, 'saveIntermediate');
if ~(islogical(cfg.saveIntermediate) || ...
        (isnumeric(cfg.saveIntermediate) && isscalar(cfg.saveIntermediate) && any(cfg.saveIntermediate == [0 1])))
    error('validate_config:invalidField', ...
        'cfg.saveIntermediate must be logical true/false (or numeric 0/1).');
end

requireField(cfg, 'captureConsoleOutput');
if ~(islogical(cfg.captureConsoleOutput) || ...
        (isnumeric(cfg.captureConsoleOutput) && isscalar(cfg.captureConsoleOutput) && any(cfg.captureConsoleOutput == [0 1])))
    error('validate_config:invalidField', ...
        'cfg.captureConsoleOutput must be logical true/false (or numeric 0/1).');
end

requireField(cfg, 'randomSeed');
if ~isnumeric(cfg.randomSeed) || ~isscalar(cfg.randomSeed) ...
        || cfg.randomSeed < 0 || cfg.randomSeed ~= floor(cfg.randomSeed)
    error('validate_config:invalidField', ...
        'cfg.randomSeed must be a non-negative integer (got %g).', ...
        cfg.randomSeed);
end

% -------------------------------------------------------------------------
% 12. Machine model (optional — pass [] to skip)
% -------------------------------------------------------------------------
if ~isempty(machine)
    validate_machine_model(machine);
end

end % validate_config

% =========================================================================
% LOCAL HELPER FUNCTIONS
% =========================================================================

function requireField(cfg, fname)
% Throws if cfg.(fname) does not exist.
if ~isfield(cfg, fname)
    error('validate_config:missingField', ...
        'Required field cfg.%s is missing.', fname);
end
end

function requireSubField(cfg, parent, fname)
% Throws if cfg.(parent).(fname) does not exist.
if ~isfield(cfg, parent)
    error('validate_config:missingField', ...
        'Required sub-struct cfg.%s is missing.', parent);
end
if ~isfield(cfg.(parent), fname)
    error('validate_config:missingField', ...
        'Required field cfg.%s.%s is missing.', parent, fname);
end
end

function validateField_string(cfg, fname, allowedValues)
% Validates that cfg.(fname) is a non-empty string.
% If allowedValues is non-empty, the value must be a member.
requireField(cfg, fname);
val = cfg.(fname);
if ~ischar(val) || isempty(strtrim(val))
    error('validate_config:invalidField', ...
        'cfg.%s must be a non-empty string.', fname);
end
if ~isempty(allowedValues) && ~ismember(val, allowedValues)
    error('validate_config:invalidValue', ...
        'cfg.%s = ''%s'' is not in the allowed set {%s}.', ...
        fname, val, strjoin(allowedValues, ', '));
end
end

function validateField_scalar_numeric(cfg, fname, lo, hi)
% Validates cfg.(fname) is a finite scalar in [lo, hi].
requireField(cfg, fname);
val = cfg.(fname);
if ~isnumeric(val) || ~isscalar(val)
    error('validate_config:invalidField', ...
        'cfg.%s must be a numeric scalar.', fname);
end
if ~isfinite(val) && ~isinf(hi)
    error('validate_config:invalidField', ...
        'cfg.%s must be a finite numeric scalar.', fname);
end
if val < lo || val > hi
    error('validate_config:outOfRange', ...
        'cfg.%s = %g is out of valid range [%g, %g].', fname, val, lo, hi);
end
end

function validateField_positive_integer(cfg, fname)
% Validates cfg.(fname) is a positive integer.
requireField(cfg, fname);
val = cfg.(fname);
if ~isnumeric(val) || ~isscalar(val) || val < 1 || val ~= floor(val)
    error('validate_config:invalidField', ...
        'cfg.%s must be a positive integer (got %g).', fname, val);
end
end

function validateField_positive_scalar(cfg, fname)
% Validates cfg.(fname) is a strictly positive finite scalar.
requireField(cfg, fname);
val = cfg.(fname);
if ~isnumeric(val) || ~isscalar(val) || ~isfinite(val) || val <= 0
    error('validate_config:invalidField', ...
        'cfg.%s must be a positive finite scalar (got %g).', fname, val);
end
end

function validateDoseGridField(doseGrid, fname)
% Validates doseGrid.<fname> is empty or a positive finite scalar.
if ~isfield(doseGrid, fname)
    error('validate_config:missingField', ...
        'Required field cfg.doseGrid.%s is missing.', fname);
end
val = doseGrid.(fname);
if isempty(val)
    return;
end
if ~isnumeric(val) || ~isscalar(val) || ~isfinite(val) || val <= 0
    error('validate_config:invalidField', ...
        'cfg.doseGrid.%s must be empty or a positive finite scalar.', fname);
end
end

function validateSubField_scalar_numeric_sub(cfg, parent, fname, lo, hi)
% Validates cfg.(parent).(fname) is a finite scalar in [lo, hi].
requireSubField(cfg, parent, fname);
val = cfg.(parent).(fname);
if ~isnumeric(val) || ~isscalar(val)
    error('validate_config:invalidField', ...
        'cfg.%s.%s must be a numeric scalar.', parent, fname);
end
if val < lo || val > hi
    error('validate_config:outOfRange', ...
        'cfg.%s.%s = %g is out of valid range [%g, %g].', ...
        parent, fname, val, lo, hi);
end
end

function validateSubField_positive_integer(cfg, parent, fname)
% Validates cfg.(parent).(fname) is a positive integer.
requireSubField(cfg, parent, fname);
val = cfg.(parent).(fname);
if ~isnumeric(val) || ~isscalar(val) || val < 1 || val ~= floor(val)
    error('validate_config:invalidField', ...
        'cfg.%s.%s must be a positive integer (got %g).', parent, fname, val);
end
end

% =========================================================================
% MACHINE MODEL VALIDATION (sub-function)
% =========================================================================

function validate_machine_model(machine)
% validate_machine_model  Validates the machine timing model struct.

% --- Energy switching ---
requireMachineField(machine, 'energy', 'upBase');
requireMachineField(machine, 'energy', 'downBase');
requireMachineField(machine, 'energy', 'perMeVUp');
requireMachineField(machine, 'energy', 'perMeVDown');
validateMachineScalar(machine.energy.upBase,    'machine.energy.upBase',    0, Inf);
validateMachineScalar(machine.energy.downBase,  'machine.energy.downBase',  0, Inf);
validateMachineScalar(machine.energy.perMeVUp,  'machine.energy.perMeVUp',  0, Inf);
validateMachineScalar(machine.energy.perMeVDown,'machine.energy.perMeVDown',0, Inf);

% --- Spot switching ---
requireMachineField(machine, 'spot', 'switchTime');
validateMachineScalar(machine.spot.switchTime, 'machine.spot.switchTime', 0, Inf);

% --- Irradiation ---
requireMachineField(machine, 'irradiation', 'muRate');
validateMachineScalar(machine.irradiation.muRate, 'machine.irradiation.muRate', 0, Inf);
if machine.irradiation.muRate <= 0
    error('validate_config:invalidField', ...
        'machine.irradiation.muRate must be strictly positive.');
end

requireMachineField(machine, 'irradiation', 'weightUnit');
if ~ischar(machine.irradiation.weightUnit)
    error('validate_config:invalidField', ...
        'machine.irradiation.weightUnit must be a string.');
end

% --- Gantry ---
requireMachineField(machine, 'gantry', 'maxSpeed');
requireMachineField(machine, 'gantry', 'maxAccel');
requireMachineField(machine, 'gantry', 'toleranceWindow');
validateMachineScalar(machine.gantry.maxSpeed,        'machine.gantry.maxSpeed',        0, Inf);
validateMachineScalar(machine.gantry.maxAccel,        'machine.gantry.maxAccel',        0, Inf);
validateMachineScalar(machine.gantry.toleranceWindow, 'machine.gantry.toleranceWindow', 0, 360);
if machine.gantry.maxSpeed <= 0
    error('validate_config:invalidField', ...
        'machine.gantry.maxSpeed must be strictly positive.');
end
if machine.gantry.maxAccel <= 0
    error('validate_config:invalidField', ...
        'machine.gantry.maxAccel must be strictly positive.');
end

% --- Buffer ---
requireMachineField(machine, 'buffer', 'perLayer');
requireMachineField(machine, 'buffer', 'perControlPoint');
validateMachineScalar(machine.buffer.perLayer,        'machine.buffer.perLayer',        0, Inf);
validateMachineScalar(machine.buffer.perControlPoint, 'machine.buffer.perControlPoint', 0, Inf);

% --- Validation metadata ---
if ~isfield(machine, 'validation') || ~isfield(machine.validation, 'status')
    error('validate_config:missingField', ...
        'machine.validation.status is missing.');
end
if ~ischar(machine.validation.status)
    error('validate_config:invalidField', ...
        'machine.validation.status must be a string.');
end
allowedStatus = {'proxy_unvalidated','literature_sourced','vendor_spec','validated'};
if ~ismember(machine.validation.status, allowedStatus)
    error('validate_config:invalidValue', ...
        ['machine.validation.status = ''%s'' is not recognised. ', ...
         'Allowed values: %s.'], ...
        machine.validation.status, strjoin(allowedStatus, ', '));
end

end % validate_machine_model

function requireMachineField(machine, parent, fname)
if ~isfield(machine, parent)
    error('validate_config:missingField', ...
        'Required machine sub-struct machine.%s is missing.', parent);
end
if ~isfield(machine.(parent), fname)
    error('validate_config:missingField', ...
        'Required field machine.%s.%s is missing.', parent, fname);
end
end

function validateMachineScalar(val, label, lo, hi)
if ~isnumeric(val) || ~isscalar(val)
    error('validate_config:invalidField', ...
        '%s must be a numeric scalar.', label);
end
if val < lo || val > hi
    error('validate_config:outOfRange', ...
        '%s = %g is out of valid range [%g, %g].', label, val, lo, hi);
end
end
