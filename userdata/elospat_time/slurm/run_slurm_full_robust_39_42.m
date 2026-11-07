function run_slurm_full_robust_39_42()
% run_slurm_full_robust_39_42  SLURM entrypoint for robust Tasks 39-42.

rootDir = getenv_required('ELOSPAT_ROOT');
matRadDir = getenv_required('MATRAD_ROOT');
nominalRoot = getenv_required('NOMINAL_ROOT');
robustRoot = getenv_required('ROBUST_ROOT');
task39ArtifactPointer = getenv_required('TASK39_NOMINAL_ARTIFACT_POINTER');
saveFullDij = getenv_logical('SAVE_FULL_DIJ', false);

cd(rootDir);
addpath(matRadDir);
matRad_rc;
addpath(rootDir);
addpath(fullfile(rootDir, 'planning'));
addpath(fullfile(rootDir, 'evaluation'));

ensure_dir(robustRoot);

fprintf('\nSLURM full robust Tasks 39-42 entrypoint\n');
fprintf('Root: %s\n', rootDir);
fprintf('Nominal root: %s\n', nominalRoot);
fprintf('Robust root: %s\n', robustRoot);
fprintf('Save full robust dij: %s\n\n', mat2str(saveFullDij));

report = struct();
report.timestamp = datetime('now');
report.nominalRoot = nominalRoot;
report.robustRoot = robustRoot;
report.saveFullDij = saveFullDij;

% Task 39 robust evaluation from the primary nominal artifact.
fprintf('\n===== Robust Task 39 =====\n');
task39Artifact = read_pointer(task39ArtifactPointer);
task39RobustOut = fullfile(robustRoot, 'task39');
ensure_dir(task39RobustOut);
report.task39 = empty_record();
report.task39.nominalArtifactFile = task39Artifact;
report.task39.outputDir = task39RobustOut;
try
    task39RobustOutputs = run_task28_30_from_nominal(task39Artifact, task39RobustOut, saveFullDij);
    report.task39 = build_robust_record(task39Artifact, task39RobustOut, task39RobustOutputs);
catch ME
    report.task39.status = 'failed';
    report.task39.errorMessage = ME.message;
    fprintf(2, 'Robust Task 39 FAILED: %s\n', ME.message);
end

% Task 40 robust follow-up: evaluate every successful nominal sweep artifact.
fprintf('\n===== Robust Task 40 =====\n');
task40NominalRoot = fullfile(nominalRoot, 'task40');
task40RobustRoot = fullfile(robustRoot, 'task40');
task40Artifacts = find_nominal_artifacts(task40NominalRoot);
task40RobustRecords = run_robust_artifact_set(task40Artifacts, task40RobustRoot, saveFullDij);
report.task40 = struct( ...
    'nominalRoot', task40NominalRoot, ...
    'robustRoot', task40RobustRoot, ...
    'nArtifacts', numel(task40Artifacts), ...
    'records', task40RobustRecords);

% Task 41 robust follow-up: evaluate every successful arc-spacing artifact.
fprintf('\n===== Robust Task 41 =====\n');
task41NominalRoot = fullfile(nominalRoot, 'task41');
task41RobustRoot = fullfile(robustRoot, 'task41');
task41Artifacts = find_nominal_artifacts(task41NominalRoot);
task41RobustRecords = run_robust_artifact_set(task41Artifacts, task41RobustRoot, saveFullDij);
report.task41 = struct( ...
    'nominalRoot', task41NominalRoot, ...
    'robustRoot', task41RobustRoot, ...
    'nArtifacts', numel(task41Artifacts), ...
    'records', task41RobustRecords);

% Task 42 documentation/reporting check for the robust branch.
fprintf('\n===== Robust Task 42 =====\n');
task42RobustOut = fullfile(robustRoot, 'task42');
ensure_dir(task42RobustOut);
task42Outputs = run_task42_documentation_check(task42RobustOut);
report.task42 = struct( ...
    'status', 'completed', ...
    'outputDir', task42RobustOut, ...
    'files', task42Outputs.files, ...
    'passed', task42Outputs.report.passed);

reportFile = fullfile(robustRoot, 'full_robust_39_42_report.mat');
save(reportFile, 'report', '-v7.3');

fprintf('\nFull robust Tasks 39-42 complete.\n');
fprintf('Report: %s\n\n', reportFile);
end

function records = run_robust_artifact_set(artifactFiles, outputRoot, saveFullDij)
ensure_dir(outputRoot);
records = repmat(empty_record(), numel(artifactFiles), 1);

if isempty(artifactFiles)
    fprintf('No nominal artifacts found under %s\n', outputRoot);
    return;
end

for i = 1:numel(artifactFiles)
    artifactFile = artifactFiles{i};
    [~, artifactBase] = fileparts(artifactFile);
    runOut = fullfile(outputRoot, sanitize_path_token(artifactBase));
    ensure_dir(runOut);

    fprintf('\nRobust artifact %d/%d\n', i, numel(artifactFiles));
    fprintf('Artifact: %s\n', artifactFile);
    fprintf('Output: %s\n', runOut);

    records(i).nominalArtifactFile = artifactFile;
    records(i).outputDir = runOut;
    try
        outputs = run_task28_30_from_nominal(artifactFile, runOut, saveFullDij);
        records(i) = build_robust_record(artifactFile, runOut, outputs);
    catch ME
        records(i).status = 'failed';
        records(i).errorMessage = ME.message;
        fprintf(2, '  FAILED: %s\n', ME.message);
    end
end
end

function record = build_robust_record(artifactFile, outputDir, outputs)
record = empty_record();
record.status = 'completed';
record.errorMessage = '';
record.nominalArtifactFile = artifactFile;
record.outputDir = outputDir;
record.matFile = outputs.matFile;
record.dvhPlot = outputs.dvhPlot;
record.summaryPlot = outputs.summaryPlot;
record.runtimeSeconds = outputs.runtimeSeconds;
record.nScenarios = outputs.nScenarios;
record.passRate = outputs.passRate;
record.warmStartSource = outputs.warmStartSource;
end

function record = empty_record()
record = struct( ...
    'status', '', ...
    'errorMessage', '', ...
    'nominalArtifactFile', '', ...
    'outputDir', '', ...
    'matFile', '', ...
    'dvhPlot', '', ...
    'summaryPlot', '', ...
    'runtimeSeconds', NaN, ...
    'nScenarios', NaN, ...
    'passRate', NaN, ...
    'warmStartSource', '');
end

function artifactFiles = find_nominal_artifacts(rootDir)
artifactFiles = {};
if exist(rootDir, 'dir') ~= 7
    return;
end
artifactFiles = find_nominal_artifacts_recursive(rootDir, artifactFiles);
artifactFiles = sort(artifactFiles(:));
end

function artifactFiles = find_nominal_artifacts_recursive(currentDir, artifactFiles)
listing = dir(currentDir);
for i = 1:numel(listing)
    name = listing(i).name;
    if strcmp(name, '.') || strcmp(name, '..')
        continue;
    end

    fullPath = fullfile(currentDir, name);
    if listing(i).isdir
        artifactFiles = find_nominal_artifacts_recursive(fullPath, artifactFiles);
    elseif ends_with(name, '_nominal_artifact.mat')
        artifactFiles{end + 1, 1} = fullPath; %#ok<AGROW>
    end
end
end

function tf = ends_with(textValue, suffix)
textValue = char(textValue);
suffix = char(suffix);
tf = numel(textValue) >= numel(suffix) && strcmp(textValue(end-numel(suffix)+1:end), suffix);
end

function token = sanitize_path_token(value)
token = regexprep(char(value), '[^A-Za-z0-9_.-]', '_');
end

function value = getenv_required(name)
value = getenv(name);
if isempty(value)
    error('run_slurm_full_robust_39_42:missingEnv', ...
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
    error('run_slurm_full_robust_39_42:invalidLogicalEnv', ...
        'Environment variable %s must be true/false, got: %s', name, raw);
end
end

function ensure_dir(pathValue)
if exist(pathValue, 'dir') ~= 7
    mkdir(pathValue);
end
end

function artifactFile = read_pointer(pointerFile)
if exist(pointerFile, 'file') ~= 2
    error('run_slurm_full_robust_39_42:missingPointer', ...
        'Artifact pointer does not exist: %s', pointerFile);
end

fid = fopen(pointerFile, 'r');
if fid < 0
    error('run_slurm_full_robust_39_42:cannotReadPointer', ...
        'Could not read artifact pointer: %s', pointerFile);
end
cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>

line = fgetl(fid);
if ~ischar(line)
    line = '';
end
artifactFile = strtrim(line);
if isempty(artifactFile)
    error('run_slurm_full_robust_39_42:emptyPointer', ...
        'Artifact pointer is empty: %s', pointerFile);
end
if exist(artifactFile, 'file') ~= 2
    error('run_slurm_full_robust_39_42:missingArtifact', ...
        'Artifact from pointer does not exist: %s', artifactFile);
end
end
