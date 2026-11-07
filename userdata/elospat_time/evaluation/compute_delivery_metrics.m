function metrics = compute_delivery_metrics(varargin)
% compute_delivery_metrics  Summarize delivery complexity and BDT metrics.
%
% USAGE
%   metrics = compute_delivery_metrics(w, map, groups, cfg, machine)
%   metrics = compute_delivery_metrics(seq, machine)
%   metrics = compute_delivery_metrics(seq, machine, cfg)
%
% OUTPUT
%   metrics - struct with fields:
%             .nControlPoints
%             .nActiveLayers
%             .nActiveSpots
%             .nSwitchUp
%             .totalSwitchMagnitude
%             .T_energy, .T_spot, .T_irradiation, .T_gantry, .T_buffer
%             .T_total
%             .spotWeightMin, .spotWeightMedian, .spotWeightMax
%             .nSubMMU, .pctSubMMU
%             .mmuThreshold

narginchk(2, 5);

[seq, activeSpotWeights, cfg, machine] = normalize_inputs(varargin{:});
bdtMetrics = estimate_bdt(seq, machine);

metrics = struct();
metrics.nControlPoints = count_unique_control_points(seq);
metrics.nActiveLayers = bdtMetrics.nActiveLayers;
metrics.nActiveSpots = bdtMetrics.nActiveSpots;
metrics.nSwitchUp = bdtMetrics.nSwitchUp;
metrics.totalSwitchMagnitude = bdtMetrics.totalSwitchMagnitude;
metrics.T_energy = bdtMetrics.T_energy;
metrics.T_spot = bdtMetrics.T_spot;
metrics.T_irradiation = bdtMetrics.T_irradiation;
metrics.T_gantry = bdtMetrics.T_gantry;
metrics.T_buffer = bdtMetrics.T_buffer;
metrics.T_total = bdtMetrics.T_total;

metrics.mmuThreshold = resolve_mmu_threshold(cfg);

if isempty(activeSpotWeights)
    metrics.spotWeightMin = NaN;
    metrics.spotWeightMedian = NaN;
    metrics.spotWeightMax = NaN;
    metrics.nSubMMU = 0;
    metrics.pctSubMMU = NaN;
    return;
end

activeSpotWeights = reshape(double(activeSpotWeights), [], 1);
metrics.spotWeightMin = min(activeSpotWeights);
metrics.spotWeightMedian = median(activeSpotWeights);
metrics.spotWeightMax = max(activeSpotWeights);

if isfinite(metrics.mmuThreshold)
    metrics.nSubMMU = sum(activeSpotWeights < metrics.mmuThreshold);
    metrics.pctSubMMU = 100 * metrics.nSubMMU / numel(activeSpotWeights);
else
    metrics.nSubMMU = NaN;
    metrics.pctSubMMU = NaN;
end

end

function [seq, activeSpotWeights, cfg, machine] = normalize_inputs(varargin)
cfg = struct();

if isnumeric(varargin{1}) && nargin == 5
    w = reshape(double(varargin{1}), [], 1);
    map = varargin{2};
    groups = varargin{3};
    cfg = varargin{4};
    machine = varargin{5};

    validate_weight_mode_inputs(w, map, groups, cfg, machine);
    seq = extract_delivery_sequence(w, map, groups, cfg);
    activeSpotWeights = w(w > double(cfg.spotThreshold));
    return;
end

if isstruct(varargin{1}) && nargin >= 2
    seq = varargin{1};
    machine = varargin{2};
    if nargin >= 3
        cfg = varargin{3};
    end

    validate_sequence_mode_inputs(seq, machine, cfg);
    activeSpotWeights = collect_spot_weights(seq);
    return;
end

error('compute_delivery_metrics:invalidInput', ...
    ['Expected either (w, map, groups, cfg, machine) or ' ...
     '(seq, machine [, cfg]).']);
end

function validate_weight_mode_inputs(w, map, groups, cfg, machine)
if ~isnumeric(w) || ~isvector(w) || isempty(w)
    error('compute_delivery_metrics:invalidInput', ...
        'w must be a non-empty numeric vector.');
end
if ~isstruct(map)
    error('compute_delivery_metrics:invalidInput', ...
        'map must be a struct array.');
end
if ~isstruct(groups)
    error('compute_delivery_metrics:invalidInput', ...
        'groups must be a struct array.');
end
if ~isstruct(cfg) || ~isfield(cfg, 'spotThreshold') || ~isfield(cfg, 'layerThreshold')
    error('compute_delivery_metrics:invalidInput', ...
        'cfg.spotThreshold and cfg.layerThreshold are required.');
end
if ~isstruct(machine)
    error('compute_delivery_metrics:invalidInput', ...
        'machine must be a struct.');
end
end

function validate_sequence_mode_inputs(seq, machine, cfg)
if ~isstruct(seq)
    error('compute_delivery_metrics:invalidInput', ...
        'seq must be a struct array.');
end
if ~isstruct(machine)
    error('compute_delivery_metrics:invalidInput', ...
        'machine must be a struct.');
end
if ~isstruct(cfg)
    error('compute_delivery_metrics:invalidInput', ...
        'cfg must be a struct when provided.');
end
end

function activeSpotWeights = collect_spot_weights(seq)
if isempty(seq)
    activeSpotWeights = zeros(0, 1);
    return;
end

activeSpotWeights = zeros(0, 1);
for i = 1:numel(seq)
    if ~isfield(seq(i), 'spotWeights') || isempty(seq(i).spotWeights)
        continue;
    end
    activeSpotWeights = [activeSpotWeights; reshape(double(seq(i).spotWeights), [], 1)]; %#ok<AGROW>
end
end

function nControlPoints = count_unique_control_points(seq)
if isempty(seq)
    nControlPoints = 0;
else
    nControlPoints = numel(unique([seq.cpIndex]', 'stable'));
end
end

function mmuThreshold = resolve_mmu_threshold(cfg)
mmuThreshold = NaN;
if isfield(cfg, 'mmuThreshold') && isnumeric(cfg.mmuThreshold) && ...
        isscalar(cfg.mmuThreshold) && isfinite(cfg.mmuThreshold) && cfg.mmuThreshold >= 0
    mmuThreshold = double(cfg.mmuThreshold);
end
end
