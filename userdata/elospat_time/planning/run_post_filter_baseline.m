function [result, spotMask] = run_post_filter_baseline(dij, cst, pln, wInit, varargin)
% run_post_filter_baseline  Threshold-based post-filter baseline (Method 4).
%
% USAGE
%   [result, spotMask] = run_post_filter_baseline(dij, cst, pln, wInit, cfg)
%   [result, spotMask] = run_post_filter_baseline(dij, cst, pln, wInit, map, groups, cfg, machine)

narginchk(5, 8);

[map, groups, cfg, machine] = normalize_optional_inputs(varargin{:});
validate_inputs(dij, cst, pln, wInit, cfg);

nBixels = double(dij.totalNumOfBixels);
wInit = reshape(double(wInit), [], 1);
threshold = resolve_post_filter_threshold(cfg);

spotMask = wInit >= threshold;
if ~any(spotMask)
    [~, maxIx] = max(wInit);
    spotMask(maxIx) = true;
end

[dijReduced, ~, wReduced] = cut_dij_by_mask(dij, spotMask, [], wInit);
[resultGui, optimizer] = reoptimize_reduced_problem(dijReduced, cst, pln, wReduced);
reducedWeights = extract_result_weights(resultGui, sum(spotMask));

wFull = zeros(nBixels, 1);
wFull(spotMask) = reducedWeights;

result = struct();
result.w = wFull;
result.resultGUI = resultGui;
result.optimizer = optimizer;
result.spotMask = reshape(logical(spotMask), [], 1);
result.thresholdUsed = threshold;
result.retainedFraction = sum(spotMask) / nBixels;
result.reductionSummary = struct( ...
    'nOriginalSpots', nBixels, ...
    'nRetainedSpots', sum(spotMask), ...
    'nRemovedSpots', nBixels - sum(spotMask));

if ~isempty(map) && ~isempty(groups) && ~isempty(machine)
    result.bdtMetrics = compute_delivery_metrics(wFull, map, groups, cfg, machine);
else
    result.bdtMetrics = struct();
end
end

function [map, groups, cfg, machine] = normalize_optional_inputs(varargin)
map = [];
groups = [];
machine = [];

if numel(varargin) == 1
    cfg = varargin{1};
    return;
end

if numel(varargin) == 4
    map = varargin{1};
    groups = varargin{2};
    cfg = varargin{3};
    machine = varargin{4};
    return;
end

error('run_post_filter_baseline:invalidInput', ...
    'Expected either cfg or (map, groups, cfg, machine) after wInit.');
end

function validate_inputs(dij, cst, pln, wInit, cfg)
if ~isstruct(dij) || ~isfield(dij, 'totalNumOfBixels')
    error('run_post_filter_baseline:invalidInput', ...
        'dij.totalNumOfBixels is required.');
end
if ~iscell(cst) || size(cst, 2) < 6
    error('run_post_filter_baseline:invalidInput', ...
        'cst must be a matRad cst cell array.');
end
if ~isstruct(pln)
    error('run_post_filter_baseline:invalidInput', ...
        'pln must be a struct.');
end
if ~isnumeric(wInit) || ~isvector(wInit) || numel(wInit) ~= double(dij.totalNumOfBixels)
    error('run_post_filter_baseline:invalidInput', ...
        'wInit must match dij.totalNumOfBixels.');
end
if ~isstruct(cfg)
    error('run_post_filter_baseline:invalidInput', ...
        'cfg must be a struct.');
end
end

function threshold = resolve_post_filter_threshold(cfg)
thresholdCandidates = [];

if isfield(cfg, 'postFilter') && isstruct(cfg.postFilter) && ...
        isfield(cfg.postFilter, 'threshold') && isnumeric(cfg.postFilter.threshold) && ...
        isscalar(cfg.postFilter.threshold) && isfinite(cfg.postFilter.threshold)
    thresholdCandidates(end + 1, 1) = double(cfg.postFilter.threshold); %#ok<AGROW>
end

if isfield(cfg, 'mmuThreshold') && isnumeric(cfg.mmuThreshold) && isscalar(cfg.mmuThreshold) && isfinite(cfg.mmuThreshold)
    thresholdCandidates(end + 1, 1) = double(cfg.mmuThreshold); %#ok<AGROW>
end

if isfield(cfg, 'spotThreshold') && isnumeric(cfg.spotThreshold) && isscalar(cfg.spotThreshold) && isfinite(cfg.spotThreshold)
    thresholdCandidates(end + 1, 1) = double(cfg.spotThreshold); %#ok<AGROW>
end

thresholdCandidates = thresholdCandidates(thresholdCandidates >= 0);
if isempty(thresholdCandidates)
    threshold = 0;
else
    threshold = max(thresholdCandidates);
end
end

function wResult = extract_result_weights(resultGui, expectedLength)
if ~isstruct(resultGui) || ~isfield(resultGui, 'w') || isempty(resultGui.w)
    error('run_post_filter_baseline:missingWeights', ...
        'Reduced optimization result did not return result.w.');
end
wResult = reshape(double(resultGui.w), [], 1);
if numel(wResult) ~= expectedLength
    error('run_post_filter_baseline:dimensionMismatch', ...
        'Reduced result.w length %d does not match expected length %d.', ...
        numel(wResult), expectedLength);
end
end
