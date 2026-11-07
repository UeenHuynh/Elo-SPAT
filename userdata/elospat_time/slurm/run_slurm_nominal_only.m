function run_slurm_nominal_only()
% run_slurm_nominal_only  SLURM entrypoint for nominal-only RMS-ELO-SPAT.

rootDir = getenv_required('ELOSPAT_ROOT');
matRadDir = getenv_required('MATRAD_ROOT');
nominalOut = getenv_required('NOMINAL_OUT');
artifactPointer = getenv_required('NOMINAL_ARTIFACT_POINTER');

arcSpacing = getenv_numeric('ARC_SPACING', 20);

cd(rootDir);
addpath(matRadDir);
matRad_rc;
addpath(rootDir);
addpath(fullfile(rootDir, 'planning'));

if exist(nominalOut, 'dir') ~= 7
    mkdir(nominalOut);
end

cfg = struct();
cfg.arcSpacing = arcSpacing;
cfg.outputDir = nominalOut;
cfg.captureConsoleOutput = false;
cfg.runtime = struct( ...
    'headless', true, ...
    'forceSoftwareOpenGL', false);
cfg.robust = struct( ...
    'enabled', false, ...
    'mode', 'COWC', ...
    'setupMM', [3 3 3], ...
    'rangePercent', 3.5);

fprintf('\nSLURM nominal-only entrypoint\n');
fprintf('Root: %s\n', rootDir);
fprintf('Output: %s\n', nominalOut);
fprintf('Arc spacing: %.2f\n\n', arcSpacing);

outputs = run_rms_elospat_main(cfg);
nominalArtifactFile = outputs.files.nominalArtifactFile;

fid = fopen(artifactPointer, 'w');
if fid < 0
    error('run_slurm_nominal_only:cannotWritePointer', ...
        'Could not write nominal artifact pointer: %s', artifactPointer);
end
cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid, '%s\n', nominalArtifactFile);

fprintf('\nNominal artifact pointer written: %s\n', artifactPointer);
fprintf('Nominal artifact: %s\n\n', nominalArtifactFile);
end

function value = getenv_required(name)
value = getenv(name);
if isempty(value)
    error('run_slurm_nominal_only:missingEnv', ...
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
    error('run_slurm_nominal_only:invalidNumericEnv', ...
        'Environment variable %s must be numeric, got: %s', name, raw);
end
end
