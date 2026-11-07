function run_slurm_robust_from_nominal()
% run_slurm_robust_from_nominal  SLURM entrypoint for robust evaluation.

rootDir = getenv_required('ELOSPAT_ROOT');
matRadDir = getenv_required('MATRAD_ROOT');
robustOut = getenv_required('ROBUST_OUT');
artifactPointer = getenv_required('NOMINAL_ARTIFACT_POINTER');
saveFullDij = getenv_logical('SAVE_FULL_DIJ', false);

cd(rootDir);
addpath(matRadDir);
matRad_rc;
addpath(rootDir);
addpath(fullfile(rootDir, 'planning'));

if exist(robustOut, 'dir') ~= 7
    mkdir(robustOut);
end

nominalArtifactFile = read_artifact_pointer(artifactPointer);
if exist(nominalArtifactFile, 'file') ~= 2
    error('run_slurm_robust_from_nominal:missingNominalArtifact', ...
        'Nominal artifact does not exist: %s', nominalArtifactFile);
end

fprintf('\nSLURM robust-from-nominal entrypoint\n');
fprintf('Root: %s\n', rootDir);
fprintf('Nominal artifact: %s\n', nominalArtifactFile);
fprintf('Output: %s\n', robustOut);
fprintf('Save full robust dij: %s\n\n', mat2str(saveFullDij));

outputs = run_task28_30_from_nominal(nominalArtifactFile, robustOut, saveFullDij);

report = struct();
report.nominalArtifactFile = nominalArtifactFile;
report.robustOutputs = outputs;
report.timestamp = datetime('now');
report.saveFullDij = saveFullDij;

reportFile = fullfile(robustOut, 'robust_from_nominal_report.mat');
save(reportFile, 'report', '-v7.3');

fprintf('\nRobust-from-nominal report saved: %s\n', reportFile);
fprintf('Robust MAT: %s\n', outputs.matFile);
fprintf('Robust warm-start source: %s\n\n', outputs.warmStartSource);
end

function value = getenv_required(name)
value = getenv(name);
if isempty(value)
    error('run_slurm_robust_from_nominal:missingEnv', ...
        'Required environment variable %s is not set.', name);
end
end

function value = getenv_logical(name, defaultValue)
raw = getenv(name);
if isempty(raw)
    value = defaultValue;
    return;
end

raw = lower(strtrim(raw));
if ismember(raw, {'1', 'true', 'yes', 'on'})
    value = true;
elseif ismember(raw, {'0', 'false', 'no', 'off'})
    value = false;
else
    error('run_slurm_robust_from_nominal:invalidLogicalEnv', ...
        'Environment variable %s must be true/false, got: %s', name, raw);
end
end

function nominalArtifactFile = read_artifact_pointer(pointerFile)
if exist(pointerFile, 'file') ~= 2
    error('run_slurm_robust_from_nominal:missingPointer', ...
        'Nominal artifact pointer does not exist: %s', pointerFile);
end

fid = fopen(pointerFile, 'r');
if fid < 0
    error('run_slurm_robust_from_nominal:cannotReadPointer', ...
        'Could not read nominal artifact pointer: %s', pointerFile);
end
cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>

line = fgetl(fid);
if ~ischar(line)
    line = '';
end
nominalArtifactFile = strtrim(line);
if isempty(nominalArtifactFile)
    error('run_slurm_robust_from_nominal:emptyPointer', ...
        'Nominal artifact pointer is empty: %s', pointerFile);
end
end
