function captureState = start_run_console_capture(runCtx, enabled)
% start_run_console_capture  Starts MATLAB diary capture for a run.
%
% INPUTS
%   runCtx  - run context from create_run_context()
%   enabled - optional logical flag (default: true)
%
% OUTPUT
%   captureState - struct with fields:
%                  .active
%                  .logFile
%                  .startedAt

narginchk(1, 2);

if nargin < 2 || isempty(enabled)
    enabled = true;
end

captureState = struct( ...
    'active', false, ...
    'logFile', '', ...
    'startedAt', '');

if ~enabled
    return;
end

if ~isstruct(runCtx) || ~isfield(runCtx, 'consoleLogFile') || isempty(runCtx.consoleLogFile)
    error('start_run_console_capture:invalidInput', ...
        'runCtx must contain a non-empty consoleLogFile field.');
end

consoleDir = fileparts(runCtx.consoleLogFile);
if ~isempty(consoleDir) && exist(consoleDir, 'dir') ~= 7
    mkdir(consoleDir);
end

diary(runCtx.consoleLogFile);
diary on;

captureState.active = true;
captureState.logFile = runCtx.consoleLogFile;
captureState.startedAt = datestr(now, 'yyyy-mm-dd HH:MM:SS.FFF');

end
