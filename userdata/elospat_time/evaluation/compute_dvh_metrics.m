function doseMetrics = compute_dvh_metrics(doseCube, ct, cst, cfg)
% compute_dvh_metrics  Compute DVH-derived dose metrics for all structures.
%
% INPUTS
%   doseCube - dose cube or dose vector on the evaluation grid
%   ct       - optional ct struct (used only to infer dimensions when needed)
%   cst      - matRad cst cell array
%   cfg      - configuration struct (optional, used for prescription fallback)
%
% OUTPUT
%   doseMetrics - struct with fields:
%                 .doseGrid
%                 .structures
%                 .summary
%
% Each entry in doseMetrics.structures contains:
%   .name
%   .type
%   .isTarget
%   .nVoxels
%   .prescribedDose
%   .D98
%   .D95
%   .D5
%   .V95
%   .Dmean
%   .Dmax
%   .homogeneityIndex
%   .dvhCurve

narginchk(3, 4);

if nargin < 4
    cfg = struct();
end

validate_inputs(doseCube, cst, cfg);

[doseLinear, doseForDvh] = normalize_dose_input(doseCube, ct, cst);
dvhDoseGrid = build_dvh_grid(doseLinear);
dvhStruct = compute_dvh_curves(cst, doseForDvh, dvhDoseGrid);

nStructures = size(cst, 1);
structureMetrics = repmat(struct( ...
    'name', '', ...
    'type', '', ...
    'isTarget', false, ...
    'nVoxels', 0, ...
    'prescribedDose', NaN, ...
    'D98', NaN, ...
    'D95', NaN, ...
    'D5', NaN, ...
    'V95', NaN, ...
    'Dmean', NaN, ...
    'Dmax', NaN, ...
    'homogeneityIndex', NaN, ...
    'dvhCurve', []), nStructures, 1);

for s = 1:nStructures
    indices = cst{s, 4}{1};
    name = char(string(cst{s, 2}));
    type = char(string(cst{s, 3}));
    isTarget = strcmpi(type, 'TARGET');

    structureMetrics(s).name = name;
    structureMetrics(s).type = type;
    structureMetrics(s).isTarget = isTarget;
    structureMetrics(s).dvhCurve = reshape(double(dvhStruct(s).volumePoints), 1, []);

    if isempty(indices)
        continue;
    end

    doseValues = reshape(double(doseLinear(indices)), [], 1);
    doseValues = doseValues(isfinite(doseValues));

    structureMetrics(s).nVoxels = numel(doseValues);
    if isempty(doseValues)
        continue;
    end

    prescribedDose = infer_prescribed_dose(cst(s, :), cfg);
    d98 = percentile_linear(doseValues, 2);
    d95 = percentile_linear(doseValues, 5);
    d5 = percentile_linear(doseValues, 95);
    dmean = mean(doseValues);
    dmax = max(doseValues);

    if isfinite(prescribedDose) && prescribedDose > 0
        v95 = 100 * mean(doseValues >= 0.95 * prescribedDose);
        hi = (d5 - d95) / prescribedDose;
    else
        v95 = NaN;
        hi = NaN;
    end

    structureMetrics(s).prescribedDose = prescribedDose;
    structureMetrics(s).D98 = d98;
    structureMetrics(s).D95 = d95;
    structureMetrics(s).D5 = d5;
    structureMetrics(s).V95 = v95;
    structureMetrics(s).Dmean = dmean;
    structureMetrics(s).Dmax = dmax;
    structureMetrics(s).homogeneityIndex = hi;
end

doseMetrics = struct();
doseMetrics.doseGrid = dvhDoseGrid;
doseMetrics.structures = structureMetrics;
doseMetrics.summary = build_summary(structureMetrics);

end

function validate_inputs(doseCube, cst, cfg)
if ~isnumeric(doseCube) || isempty(doseCube)
    error('compute_dvh_metrics:invalidInput', ...
        'doseCube must be a non-empty numeric array or vector.');
end

if ~iscell(cst) || size(cst, 2) < 6
    error('compute_dvh_metrics:invalidInput', ...
        'cst must be a matRad cst cell array with at least 6 columns.');
end

if ~isstruct(cfg)
    error('compute_dvh_metrics:invalidInput', ...
        'cfg must be a struct.');
end
end

function [doseLinear, doseForDvh] = normalize_dose_input(doseCube, ct, cst)
doseLinear = reshape(double(doseCube), [], 1);
maxIndex = find_max_cst_index(cst);

if maxIndex > numel(doseLinear)
    dims = infer_dimensions_from_ct(ct);
    if isempty(dims) || prod(double(dims)) ~= numel(doseLinear)
        error('compute_dvh_metrics:dimensionMismatch', ...
            'Dose array has %d voxels, but cst references index %d.', ...
            numel(doseLinear), maxIndex);
    end
    doseForDvh = reshape(doseLinear, dims);
else
    dims = infer_dimensions_from_ct(ct);
    if isempty(dims) || prod(double(dims)) ~= numel(doseLinear)
        doseForDvh = doseLinear;
    else
        doseForDvh = reshape(doseLinear, dims);
    end
end
end

function maxIndex = find_max_cst_index(cst)
maxIndex = 0;
for i = 1:size(cst, 1)
    if isempty(cst{i, 4}) || isempty(cst{i, 4}{1})
        continue;
    end
    maxIndex = max(maxIndex, max(double(cst{i, 4}{1}(:))));
end
end

function dims = infer_dimensions_from_ct(ct)
dims = [];

if ~isstruct(ct) || isempty(fieldnames(ct))
    return;
end

if isfield(ct, 'cubeDim') && isnumeric(ct.cubeDim) && numel(ct.cubeDim) == 3
    dims = reshape(double(ct.cubeDim), 1, []);
    return;
end

if isfield(ct, 'cubeHU') && ~isempty(ct.cubeHU) && iscell(ct.cubeHU) && ~isempty(ct.cubeHU{1})
    dims = size(ct.cubeHU{1});
    return;
end

if isfield(ct, 'cube') && isnumeric(ct.cube) && ndims(ct.cube) == 3
    dims = size(ct.cube);
end
end

function dvhDoseGrid = build_dvh_grid(doseLinear)
maxDose = max(double(doseLinear(:)));
if ~isfinite(maxDose) || maxDose < 0
    maxDose = 0;
end

dvhDoseGrid = 0:0.1:(ceil(maxDose / 0.1) * 0.1 + 0.1);
if isempty(dvhDoseGrid)
    dvhDoseGrid = 0;
end
end

function dvhStruct = compute_dvh_curves(cst, doseForDvh, dvhDoseGrid)
if exist('matRad_calcDVH', 'file') == 2
    dvhStruct = matRad_calcDVH(cst, doseForDvh, 'cum', dvhDoseGrid);
    return;
end

nStructures = size(cst, 1);
dvhStruct = repmat(struct('doseGrid', dvhDoseGrid, 'volumePoints', [], 'name', ''), nStructures, 1);
for s = 1:nStructures
    indices = cst{s, 4}{1};
    doses = reshape(double(doseForDvh(indices)), [], 1);
    curve = zeros(1, numel(dvhDoseGrid));
    for j = 1:numel(dvhDoseGrid)
        curve(j) = 100 * mean(doses >= dvhDoseGrid(j));
    end
    dvhStruct(s).doseGrid = dvhDoseGrid;
    dvhStruct(s).volumePoints = curve;
    dvhStruct(s).name = cst{s, 2};
end
end

function prescribedDose = infer_prescribed_dose(cstRow, cfg)
prescribedDose = NaN;

if isstruct(cfg) && isfield(cfg, 'prescriptionDose') && ...
        isnumeric(cfg.prescriptionDose) && isscalar(cfg.prescriptionDose) && ...
        isfinite(cfg.prescriptionDose) && cfg.prescriptionDose > 0
    prescribedDose = double(cfg.prescriptionDose);
end

objectives = cstRow{6};
if isempty(objectives)
    return;
end

candidateValues = [];
for j = 1:numel(objectives)
    candidateValues = [candidateValues; extract_objective_dose_values(objectives{j})]; %#ok<AGROW>
end

candidateValues = candidateValues(isfinite(candidateValues) & candidateValues > 0);
if ~isempty(candidateValues)
    prescribedDose = min(candidateValues);
end
end

function values = extract_objective_dose_values(objective)
values = [];

fieldCandidates = {'dose', 'doseValue', 'referenceDose', 'refDose', 'dRef'};
for i = 1:numel(fieldCandidates)
    fieldName = fieldCandidates{i};
    if has_field_or_prop(objective, fieldName)
        values = [values; flatten_positive_numeric(read_field_or_prop(objective, fieldName))]; %#ok<AGROW>
    end
end

if has_field_or_prop(objective, 'parameters')
    parameters = read_field_or_prop(objective, 'parameters');
    if iscell(parameters)
        for i = 1:numel(parameters)
            values = [values; flatten_positive_numeric(parameters{i})]; %#ok<AGROW>
        end
    else
        values = [values; flatten_positive_numeric(parameters)]; %#ok<AGROW>
    end
end
end

function tf = has_field_or_prop(value, name)
tf = (isstruct(value) && isfield(value, name)) || ...
    (isobject(value) && isprop(value, name));
end

function out = read_field_or_prop(value, name)
if isstruct(value)
    out = value.(name);
else
    out = value.(name);
end
end

function values = flatten_positive_numeric(value)
if isnumeric(value) && ~isempty(value)
    values = reshape(double(value), [], 1);
else
    values = [];
end
end

function p = percentile_linear(values, percentileValue)
values = reshape(double(values), [], 1);
values = values(isfinite(values));
if isempty(values)
    p = NaN;
    return;
end

values = sort(values, 'ascend');
n = numel(values);
if n == 1
    p = values;
    return;
end

position = 1 + (n - 1) * (percentileValue / 100);
lowerIx = floor(position);
upperIx = ceil(position);

if lowerIx == upperIx
    p = values(lowerIx);
    return;
end

alpha = position - lowerIx;
p = (1 - alpha) * values(lowerIx) + alpha * values(upperIx);
end

function summary = build_summary(structureMetrics)
targetMask = [structureMetrics.isTarget]';
oarMask = ~targetMask;

targetD95 = [structureMetrics(targetMask).D95]';
targetD98 = [structureMetrics(targetMask).D98]';
targetV95 = [structureMetrics(targetMask).V95]';
oarDmax = [structureMetrics(oarMask).Dmax]';

summary = struct();
summary.nStructures = numel(structureMetrics);
summary.nTargets = sum(targetMask);
summary.nOars = sum(oarMask);
summary.targetD95Min = safe_reduce(targetD95, @min);
summary.targetD98Min = safe_reduce(targetD98, @min);
summary.targetV95Min = safe_reduce(targetV95, @min);
summary.oarDmaxMax = safe_reduce(oarDmax, @max);
end

function value = safe_reduce(values, reducer)
values = values(isfinite(values));
if isempty(values)
    value = NaN;
else
    value = reducer(values);
end
end
