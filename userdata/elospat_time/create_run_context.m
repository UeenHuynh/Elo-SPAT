function runCtx = create_run_context(cfg, runLabel)
% create_run_context  Creates a run-specific output folder and log context.
%
% INPUTS
%   cfg      - configuration struct with cfg.outputDir and cfg.saveIntermediate
%   runLabel - optional short label for the run directory name
%
% OUTPUT
%   runCtx - struct with runId, outputDir, checkpointDir, logFile,
%            consoleLogFile, summaryFile, configSnapshotFile,
%            saveIntermediate, createdAt, and label

narginchk(1, 2);

if nargin < 2 || isempty(runLabel)
    runLabel = 'run';
end

if ~isstruct(cfg)
    error('create_run_context:invalidInput', ...
        'cfg must be a configuration struct.');
end

if ~isfield(cfg, 'outputDir') || isempty(cfg.outputDir)
    error('create_run_context:missingField', ...
        'cfg.outputDir is required to create a run context.');
end

if ~isfield(cfg, 'saveIntermediate') || isempty(cfg.saveIntermediate)
    saveIntermediate = false;
else
    saveIntermediate = logical(cfg.saveIntermediate);
end

timestamp = datestr(now, 'yyyymmdd_HHMMSSFFF');
safeLabel = sanitize_label(runLabel);
runId = sprintf('%s_%s', safeLabel, timestamp);

baseDir = cfg.outputDir;
if exist(baseDir, 'dir') ~= 7
    mkdir(baseDir);
end

outputDir = fullfile(baseDir, runId);
checkpointDir = fullfile(outputDir, 'checkpoints');

mkdir(outputDir);
mkdir(checkpointDir);

runCtx = struct();
runCtx.runId = runId;
runCtx.label = safeLabel;
runCtx.outputDir = outputDir;
runCtx.checkpointDir = checkpointDir;
runCtx.logFile = fullfile(outputDir, 'run.log');
runCtx.consoleLogFile = fullfile(outputDir, 'console.log');
runCtx.summaryFile = fullfile(outputDir, 'summary.mat');
runCtx.configSnapshotFile = fullfile(outputDir, 'config_snapshot.mat');
runCtx.saveIntermediate = saveIntermediate;
runCtx.createdAt = datestr(now, 'yyyy-mm-dd HH:MM:SS.FFF');

save(fullfile(outputDir, 'run_context.mat'), 'runCtx');
save(runCtx.configSnapshotFile, 'cfg');
append_run_log(runCtx, sprintf('Created run context in %s', outputDir), 'INFO');
append_run_log(runCtx, sprintf('Saved config snapshot: %s', runCtx.configSnapshotFile), 'INFO');

end

function safeLabel = sanitize_label(runLabel)
if isstring(runLabel)
    runLabel = char(runLabel);
end

if ~ischar(runLabel)
    error('create_run_context:invalidInput', ...
        'runLabel must be a character vector or string scalar.');
end

safeLabel = lower(strtrim(runLabel));
if isempty(safeLabel)
    safeLabel = 'run';
end

safeLabel = regexprep(safeLabel, '[^a-zA-Z0-9_-]', '_');
safeLabel = regexprep(safeLabel, '_+', '_');
safeLabel = regexprep(safeLabel, '^_+|_+$', '');

if isempty(safeLabel)
    safeLabel = 'run';
end
end
