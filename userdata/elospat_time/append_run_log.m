function logPath = append_run_log(target, message, level)
% append_run_log  Appends a timestamped line to a run log text file.
%
% INPUTS
%   target  - run context struct with .logFile, or a log file path string
%   message - text to append
%   level   - optional level string (default: 'INFO')
%
% OUTPUT
%   logPath - resolved log file path

narginchk(2, 3);

if nargin < 3 || isempty(level)
    level = 'INFO';
end

logPath = resolve_log_path(target);

logDir = fileparts(logPath);
if ~isempty(logDir) && exist(logDir, 'dir') ~= 7
    mkdir(logDir);
end

if isstring(message)
    message = char(message);
end

if ~ischar(message)
    error('append_run_log:invalidInput', ...
        'message must be a character vector or string scalar.');
end

if isstring(level)
    level = char(level);
end

timestamp = datestr(now, 'yyyy-mm-dd HH:MM:SS.FFF');
lines = regexp(message, '\r\n|\n|\r', 'split');

fid = fopen(logPath, 'a');
if fid == -1
    error('append_run_log:fileOpenFailed', ...
        'Could not open log file for writing: %s', logPath);
end

cleanupObj = onCleanup(@() fclose(fid));

for i = 1:numel(lines)
    fprintf(fid, '[%s] [%s] %s\n', timestamp, upper(level), lines{i});
end

clear cleanupObj;

end

function logPath = resolve_log_path(target)
if isstruct(target)
    if ~isfield(target, 'logFile') || isempty(target.logFile)
        error('append_run_log:missingField', ...
            'run context must contain a non-empty logFile field.');
    end
    logPath = target.logFile;
elseif isstring(target)
    logPath = char(target);
elseif ischar(target)
    logPath = target;
else
    error('append_run_log:invalidInput', ...
        'target must be a run context struct or a log file path.');
end
end
