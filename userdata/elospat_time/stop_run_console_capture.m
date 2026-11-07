function stop_run_console_capture(captureState)
% stop_run_console_capture  Stops MATLAB diary capture if owned by the caller.

narginchk(1, 1);

if ~isstruct(captureState) || ~isfield(captureState, 'active')
    error('stop_run_console_capture:invalidInput', ...
        'captureState must be the struct returned by start_run_console_capture().');
end

if ~logical(captureState.active)
    return;
end

diary off;

end
