function summaryPath = save_run_summary(runCtx, summaryData)
% save_run_summary  Persists a lightweight run summary MAT file.
%
% INPUTS
%   runCtx      - run context from create_run_context()
%   summaryData - lightweight struct describing run status/progress
%
% OUTPUT
%   summaryPath - resolved summary MAT file path

narginchk(2, 2);

if ~isstruct(runCtx) || ~isfield(runCtx, 'outputDir') || isempty(runCtx.outputDir)
    error('save_run_summary:invalidInput', ...
        'runCtx must be a valid struct returned by create_run_context().');
end

if ~isstruct(summaryData)
    error('save_run_summary:invalidInput', ...
        'summaryData must be a struct.');
end

if isfield(runCtx, 'summaryFile') && ~isempty(runCtx.summaryFile)
    summaryPath = runCtx.summaryFile;
else
    summaryPath = fullfile(runCtx.outputDir, 'summary.mat');
end

summaryDir = fileparts(summaryPath);
if ~isempty(summaryDir) && exist(summaryDir, 'dir') ~= 7
    mkdir(summaryDir);
end

savedAt = datestr(now, 'yyyy-mm-dd HH:MM:SS.FFF');
summaryData.savedAt = savedAt;

save(summaryPath, 'summaryData', 'savedAt');

end
