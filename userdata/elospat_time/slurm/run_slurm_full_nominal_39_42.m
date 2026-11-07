function run_slurm_full_nominal_39_42()
% run_slurm_full_nominal_39_42  SLURM entrypoint for nominal Tasks 39-42.

rootDir = getenv_required('ELOSPAT_ROOT');
matRadDir = getenv_required('MATRAD_ROOT');
nominalRoot = getenv_required('NOMINAL_ROOT');
task39ArtifactPointer = getenv_required('TASK39_NOMINAL_ARTIFACT_POINTER');
arcSpacing = getenv_numeric('ARC_SPACING', 20);

cd(rootDir);
addpath(matRadDir);
matRad_rc;
addpath(rootDir);
addpath(fullfile(rootDir, 'planning'));
addpath(fullfile(rootDir, 'evaluation'));

ensure_dir(nominalRoot);

fprintf('\nSLURM full nominal Tasks 39-42 entrypoint\n');
fprintf('Root: %s\n', rootDir);
fprintf('Nominal root: %s\n', nominalRoot);
fprintf('Task 39 arc spacing: %.2f\n\n', arcSpacing);

report = struct();
report.timestamp = datetime('now');
report.arcSpacing = arcSpacing;
report.nominalRoot = nominalRoot;

% Task 39 nominal integration: full nominal pipeline, robust disabled.
task39Out = fullfile(nominalRoot, 'task39');
ensure_dir(task39Out);
cfg = struct();
cfg.arcSpacing = arcSpacing;
cfg.outputDir = task39Out;
cfg.captureConsoleOutput = false;
cfg.runtime = struct('headless', true, 'forceSoftwareOpenGL', false);
cfg.robust = struct( ...
    'enabled', false, ...
    'mode', 'COWC', ...
    'setupMM', [3 3 3], ...
    'rangePercent', 3.5);

fprintf('\n===== Nominal Task 39 =====\n');
task39Outputs = run_rms_elospat_main(cfg);
task39ArtifactFile = task39Outputs.files.nominalArtifactFile;
write_pointer(task39ArtifactPointer, task39ArtifactFile);
report.task39 = struct( ...
    'status', 'completed', ...
    'outputDir', task39Out, ...
    'nominalArtifactFile', task39ArtifactFile, ...
    'matFile', task39Outputs.files.matFile);

% Task 40 nominal parameter sweep.
fprintf('\n===== Nominal Task 40 =====\n');
task40Out = fullfile(nominalRoot, 'task40');
ensure_dir(task40Out);
task40Outputs = run_task40_parameter_sweep([], [], task40Out);
report.task40 = struct( ...
    'status', 'completed', ...
    'outputDir', task40Out, ...
    'files', task40Outputs.files);

% Task 41 nominal arc-spacing validation.
fprintf('\n===== Nominal Task 41 =====\n');
task41Out = fullfile(nominalRoot, 'task41');
ensure_dir(task41Out);
task41Outputs = run_task41_arc_spacing_validation([], [], task41Out);
report.task41 = struct( ...
    'status', 'completed', ...
    'outputDir', task41Out, ...
    'files', task41Outputs.files);

% Task 42 documentation/reporting check for the nominal branch.
fprintf('\n===== Nominal Task 42 =====\n');
task42Out = fullfile(nominalRoot, 'task42');
ensure_dir(task42Out);
task42Outputs = run_task42_documentation_check(task42Out);
report.task42 = struct( ...
    'status', 'completed', ...
    'outputDir', task42Out, ...
    'files', task42Outputs.files, ...
    'passed', task42Outputs.report.passed);

reportFile = fullfile(nominalRoot, 'full_nominal_39_42_report.mat');
save(reportFile, 'report', '-v7.3');

fprintf('\nFull nominal Tasks 39-42 complete.\n');
fprintf('Report: %s\n', reportFile);
fprintf('Task 39 nominal artifact pointer: %s\n', task39ArtifactPointer);
fprintf('Task 39 nominal artifact: %s\n\n', task39ArtifactFile);
end

function value = getenv_required(name)
value = getenv(name);
if isempty(value)
    error('run_slurm_full_nominal_39_42:missingEnv', ...
        'Required environment variable %s is not set.', name);
end
end

function value = getenv_numeric(name, defaultValue)
raw = getenv(name);
if isempty(raw)
    value = defaultValue;
    return;
end
value = str2double(raw);
if ~isfinite(value)
    error('run_slurm_full_nominal_39_42:invalidNumericEnv', ...
        'Environment variable %s must be numeric, got: %s', name, raw);
end
end

function ensure_dir(pathValue)
if exist(pathValue, 'dir') ~= 7
    mkdir(pathValue);
end
end

function write_pointer(pointerFile, artifactFile)
fid = fopen(pointerFile, 'w');
if fid < 0
    error('run_slurm_full_nominal_39_42:cannotWritePointer', ...
        'Could not write artifact pointer: %s', pointerFile);
end
cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid, '%s\n', artifactFile);
end
