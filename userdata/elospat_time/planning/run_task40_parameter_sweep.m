function outputs = run_task40_parameter_sweep(sweepOverride, baseCfgOverride, outputDir)
% run_task40_parameter_sweep  Convenience wrapper for Task 40 sweep mode.
%
% Default grid:
%   - arcSpacing in [20, 30]
%   - wrapper.maxLayerIterations in [2, 4]

narginchk(0, 3);

planningDir = fileparts(mfilename('fullpath'));
workspaceDir = fileparts(planningDir);
resultsDir = fullfile(workspaceDir, 'results');

addpath(workspaceDir);

cfg = struct();
cfg.captureConsoleOutput = false;
cfg.runtime = struct('headless', true, 'forceSoftwareOpenGL', false);
cfg.robust = struct('enabled', false);
cfg.sweep = struct();
cfg.sweep.enabled = true;
cfg.sweep.parameters = {'arcSpacing', 'wrapper.maxLayerIterations'};
cfg.sweep.values = {[20; 30], [2; 4]};

if nargin >= 1 && ~isempty(sweepOverride)
    cfg.sweep = merge_override(cfg.sweep, sweepOverride);
end

if nargin >= 2 && ~isempty(baseCfgOverride)
    cfg = merge_override(cfg, baseCfgOverride);
end
if ~isfield(cfg, 'robust') || ~isstruct(cfg.robust)
    cfg.robust = struct();
end
cfg.robust.enabled = false;

if nargin >= 3 && ~isempty(outputDir)
    cfg.sweep.outputDir = outputDir;
elseif ~isfield(cfg.sweep, 'outputDir') || isempty(cfg.sweep.outputDir)
    cfg.sweep.outputDir = fullfile(resultsDir, sprintf('task40_%s', datestr(now, 'yyyy-mm-dd_HHMMSS')));
end

fprintf('\n========================================\n');
fprintf('Task 40 Parameter Sweep Runner\n');
fprintf('========================================\n');
fprintf('Sweep output directory: %s\n\n', cfg.sweep.outputDir);

outputs = run_rms_elospat_main(cfg);
end

function merged = merge_override(base, override)
merged = base;
fields = fieldnames(override);
for i = 1:numel(fields)
    key = fields{i};
    if isfield(base, key) && isstruct(base.(key)) && isstruct(override.(key))
        merged.(key) = merge_override(base.(key), override.(key));
    else
        merged.(key) = override.(key);
    end
end
end
