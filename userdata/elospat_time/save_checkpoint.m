function checkpointPath = save_checkpoint(runCtx, stageName, checkpointData, iteration)
% save_checkpoint  Saves a stage checkpoint MAT file when enabled.
%
% INPUTS
%   runCtx         - run context from create_run_context()
%   stageName      - stage identifier, e.g. 'layer_filter'
%   checkpointData - struct or value to persist
%   iteration      - optional numeric or text iteration identifier
%
% OUTPUT
%   checkpointPath - full path to the checkpoint file, or '' if skipped

narginchk(3, 4);

if ~isstruct(runCtx)
    error('save_checkpoint:invalidInput', ...
        'runCtx must be a struct returned by create_run_context().');
end

if ~isfield(runCtx, 'saveIntermediate') || ~logical(runCtx.saveIntermediate)
    checkpointPath = '';
    return;
end

if isstring(stageName)
    stageName = char(stageName);
end

if ~ischar(stageName) || isempty(strtrim(stageName))
    error('save_checkpoint:invalidInput', ...
        'stageName must be a non-empty character vector or string scalar.');
end

if nargin < 4
    iterationSuffix = '';
    iterationValue = [];
else
    iterationValue = iteration;
    iterationSuffix = format_iteration_suffix(iteration);
end

safeStage = regexprep(lower(strtrim(stageName)), '[^a-zA-Z0-9_-]', '_');
safeStage = regexprep(safeStage, '_+', '_');
safeStage = regexprep(safeStage, '^_+|_+$', '');
if isempty(safeStage)
    safeStage = 'checkpoint';
end

timestamp = datestr(now, 'yyyymmdd_HHMMSSFFF');
fileName = sprintf('%s%s_%s.mat', safeStage, iterationSuffix, timestamp);
checkpointPath = fullfile(runCtx.checkpointDir, fileName);
savedAt = datestr(now, 'yyyy-mm-dd HH:MM:SS.FFF');

save(checkpointPath, 'stageName', 'checkpointData', 'iterationValue', 'savedAt', '-v7.3');
append_run_log(runCtx, sprintf('Saved checkpoint: %s', checkpointPath), 'INFO');

end

function iterationSuffix = format_iteration_suffix(iteration)
if isstring(iteration)
    iteration = char(iteration);
end

if ischar(iteration)
    safeIter = regexprep(strtrim(iteration), '[^a-zA-Z0-9_-]', '_');
    safeIter = regexprep(safeIter, '_+', '_');
    safeIter = regexprep(safeIter, '^_+|_+$', '');
    if isempty(safeIter)
        iterationSuffix = '';
    else
        iterationSuffix = ['_' safeIter];
    end
elseif isnumeric(iteration) && isscalar(iteration) && isfinite(iteration)
    iterationSuffix = sprintf('_iter%03d', round(iteration));
else
    error('save_checkpoint:invalidInput', ...
        'iteration must be a finite scalar or a text label.');
end
end
