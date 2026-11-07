function [resultPath, runCtx] = save_plan_result(pln, stf, dij, resultGUI, meta, cfg, runCtx, artifactName)
% save_plan_result  Saves plan artifacts to a MAT file inside cfg.outputDir.
%
% INPUTS
%   pln, stf, dij, resultGUI, meta - plan artifacts to persist
%   cfg         - configuration struct with cfg.outputDir
%   runCtx      - optional run context from create_run_context()
%   artifactName - optional short artifact label (default: 'plan_result')
%
% OUTPUTS
%   resultPath - full path to the saved MAT file
%   runCtx     - run context used for the save

narginchk(6, 8);

if nargin < 7 || isempty(runCtx)
    runCtx = create_run_context(cfg, 'plan_result');
end

if nargin < 8 || isempty(artifactName)
    artifactName = 'plan_result';
end

safeName = sanitize_artifact_name(artifactName);
timestamp = datestr(now, 'yyyymmdd_HHMMSSFFF');
fileName = sprintf('%s_%s.mat', safeName, timestamp);
resultPath = fullfile(runCtx.outputDir, fileName);
savedAt = datestr(now, 'yyyy-mm-dd HH:MM:SS.FFF');

save(resultPath, 'pln', 'stf', 'dij', 'resultGUI', 'meta', 'savedAt', '-v7.3');
append_run_log(runCtx, sprintf('Saved plan result MAT: %s', resultPath), 'INFO');

end

function safeName = sanitize_artifact_name(artifactName)
if isstring(artifactName)
    artifactName = char(artifactName);
end

if ~ischar(artifactName)
    error('save_plan_result:invalidInput', ...
        'artifactName must be a character vector or string scalar.');
end

safeName = lower(strtrim(artifactName));
if isempty(safeName)
    safeName = 'plan_result';
end

safeName = regexprep(safeName, '[^a-zA-Z0-9_-]', '_');
safeName = regexprep(safeName, '_+', '_');
safeName = regexprep(safeName, '^_+|_+$', '');

if isempty(safeName)
    safeName = 'plan_result';
end
end
