function stfReduced = cut_stf_by_mask(stf, mask, map)
% cut_stf_by_mask  Remove inactive rays and bixels from stf based on a bixel mask.
%
% INPUTS
%   stf  - original steering information struct array
%   mask - logical bixel-level keep mask (N x 1), true = keep
%   map  - original bixel mapping with beamNum, rayNum, bixelNum, globalIx
%
% OUTPUT
%   stfReduced - reduced stf consistent with cut_dij_by_mask output

narginchk(3, 3);

validate_stf(stf);
nBixels = validate_map_and_mask(map, mask);
mask = reshape(logical(mask), [], 1);

stfReduced = repmat(stf(1), 0, 1);

for b = 1:numel(stf)
    beam = stf(b);
    if ~isfield(beam, 'ray') || isempty(beam.ray)
        continue;
    end

    reducedRays = repmat(beam.ray(1), 0, 1);

    for r = 1:numel(beam.ray)
        ray = beam.ray(r);
        nRayBixels = infer_ray_bixel_count(ray, b, r);
        keepLocal = build_local_keep_mask(map, mask, b, r, nRayBixels, nBixels);

        if ~any(keepLocal)
            continue;
        end

        reducedRays(end + 1, 1) = slice_ray_fields(ray, keepLocal, nRayBixels); %#ok<AGROW>
    end

    if isempty(reducedRays)
        continue;
    end

    beam.ray = reducedRays;
    if isfield(beam, 'numOfRays')
        beam.numOfRays = numel(reducedRays);
    end
    if isfield(beam, 'totalNumOfBixels')
        beam.totalNumOfBixels = sum(arrayfun(@(ray) infer_ray_bixel_count(ray, b, 1), reducedRays));
    end

    stfReduced(end + 1, 1) = beam; %#ok<AGROW>
end

end

function validate_stf(stf)
if ~isstruct(stf) || isempty(stf)
    error('cut_stf_by_mask:invalidInput', ...
        'stf must be a non-empty struct array.');
end

for b = 1:numel(stf)
    if ~isfield(stf(b), 'ray') || ~isstruct(stf(b).ray)
        error('cut_stf_by_mask:invalidInput', ...
            'stf(%d).ray must be a struct array.', b);
    end
end
end

function nBixels = validate_map_and_mask(map, mask)
requiredFields = {'globalIx', 'beamNum', 'rayNum', 'bixelNum'};

if ~isstruct(map) || isempty(map)
    error('cut_stf_by_mask:invalidInput', ...
        'map must be a non-empty struct array.');
end

for i = 1:numel(requiredFields)
    if ~isfield(map, requiredFields{i})
        error('cut_stf_by_mask:missingField', ...
            'map.%s is required.', requiredFields{i});
    end
end

nBixels = numel(map);

if ~(islogical(mask) || isnumeric(mask)) || ~isvector(mask) || numel(mask) ~= nBixels
    error('cut_stf_by_mask:dimensionMismatch', ...
        'mask length %d must match numel(map)=%d.', numel(mask), nBixels);
end

if isnumeric(mask) && (any(~isfinite(mask(:))) || any(mask(:) ~= 0 & mask(:) ~= 1))
    error('cut_stf_by_mask:invalidInput', ...
        'Numeric mask must contain only 0/1 values.');
end
end

function nRayBixels = infer_ray_bixel_count(ray, beamIx, rayIx)
if isfield(ray, 'energy') && ~isempty(ray.energy)
    nRayBixels = numel(ray.energy);
    return;
end

fieldNames = fieldnames(ray);
for i = 1:numel(fieldNames)
    value = ray.(fieldNames{i});
    if isnumeric(value) && isvector(value) && ~isempty(value)
        nRayBixels = numel(value);
        return;
    end
end

error('cut_stf_by_mask:cannotInferRaySize', ...
    'Could not infer number of bixels for stf(%d).ray(%d).', beamIx, rayIx);
end

function keepLocal = build_local_keep_mask(map, mask, beamIx, rayIx, nRayBixels, nBixels)
beamNum = reshape(double([map.beamNum]), [], 1);
rayNum = reshape(double([map.rayNum]), [], 1);
bixelNum = reshape(double([map.bixelNum]), [], 1);
globalIx = reshape(double([map.globalIx]), [], 1);

if numel(globalIx) ~= nBixels
    error('cut_stf_by_mask:invalidMap', ...
        'map.globalIx length mismatch.');
end

entryMask = beamNum == beamIx & rayNum == rayIx;
localBixels = bixelNum(entryMask);
localGlobalIx = globalIx(entryMask);

keepLocal = false(nRayBixels, 1);

if isempty(localBixels)
    return;
end

if any(localBixels < 1) || any(localBixels > nRayBixels) || any(floor(localBixels) ~= localBixels)
    error('cut_stf_by_mask:invalidMap', ...
        'map contains invalid local bixel indices for beam %d ray %d.', beamIx, rayIx);
end

keepLocal(localBixels) = mask(localGlobalIx);
end

function rayReduced = slice_ray_fields(ray, keepLocal, nRayBixels)
rayReduced = ray;
fieldNames = fieldnames(ray);

for i = 1:numel(fieldNames)
    fieldName = fieldNames{i};
    value = ray.(fieldName);
    rayReduced.(fieldName) = slice_field_value(value, keepLocal, nRayBixels, fieldName);
end
end

function valueReduced = slice_field_value(value, keepLocal, nRayBixels, fieldName)
if isempty(value)
    valueReduced = value;
    return;
end

if isnumeric(value) || islogical(value)
    if isvector(value) && numel(value) == nRayBixels
        if isrow(value)
            valueReduced = value(keepLocal.');
        else
            valueReduced = value(keepLocal);
        end
        return;
    end

    if ismatrix(value) && size(value, 1) == nRayBixels && size(value, 2) ~= nRayBixels
        valueReduced = value(keepLocal, :);
        return;
    end

    if ismatrix(value) && size(value, 2) == nRayBixels && size(value, 1) ~= nRayBixels
        valueReduced = value(:, keepLocal);
        return;
    end
end

if iscell(value)
    if isvector(value) && numel(value) == nRayBixels
        if isrow(value)
            valueReduced = value(keepLocal.');
        else
            valueReduced = value(keepLocal);
        end
        return;
    end
end

if isstruct(value)
    valueReduced = value;
    return;
end

if ischar(value) || isstring(value)
    valueReduced = value;
    return;
end

% Leave unrelated metadata unchanged.
valueReduced = value;

if isnumeric(value) && ismatrix(value) && size(value, 1) == nRayBixels && size(value, 2) == nRayBixels
    error('cut_stf_by_mask:ambiguousField', ...
        'Field %s has ambiguous %dx%d size; explicit handling is required.', ...
        fieldName, nRayBixels, nRayBixels);
end
end
