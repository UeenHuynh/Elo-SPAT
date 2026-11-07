function [dijReduced, mapReduced, wReduced] = cut_dij_by_mask(dij, mask, map, w)
% cut_dij_by_mask  Remove inactive bixels from a dij structure.
%
% INPUTS
%   dij  - original dose influence matrix struct
%   mask - logical bixel-level keep mask (N x 1), true = keep
%   map  - optional original bixel mapping
%   w    - optional original weight vector
%
% OUTPUTS
%   dijReduced - reduced dij with only active bixels retained
%   mapReduced - reduced mapping, preserving original globalIx
%   wReduced   - reduced weight vector for warm-start use

narginchk(2, 4);

if nargin < 3
    map = [];
end

if nargin < 4
    w = [];
end

validate_dij(dij);
nBixels = double(dij.totalNumOfBixels);
mask = validate_mask(mask, nBixels);
keepIdx = find(mask);

if ~isempty(map)
    validate_map(map, nBixels);
end

if ~isempty(w)
    validate_weight_vector(w, nBixels);
end

dijReduced = dij;

doseFields = {'physicalDose', 'mAlphaDose', 'mSqrtBetaDose', 'mLETDose'};
for i = 1:numel(doseFields)
    fieldName = doseFields{i};
    if isfield(dijReduced, fieldName) && ~isempty(dijReduced.(fieldName))
        dijReduced.(fieldName) = slice_dose_field(dijReduced.(fieldName), keepIdx, nBixels, fieldName);
    end
end

vectorFields = {'beamNum', 'rayNum', 'bixelNum', 'minMU', 'maxMU', 'numParticlesPerMU'};
for i = 1:numel(vectorFields)
    fieldName = vectorFields{i};
    if isfield(dijReduced, fieldName) && ~isempty(dijReduced.(fieldName))
        dijReduced.(fieldName) = slice_vector_field(dijReduced.(fieldName), keepIdx, nBixels, fieldName);
    end
end

dijReduced.totalNumOfBixels = numel(keepIdx);
dijReduced.keptBixelIndices = reshape(keepIdx, [], 1);

if ~isempty(map)
    mapReduced = map(mask);
else
    mapReduced = [];
end

if ~isempty(w)
    wReduced = w(mask);
else
    wReduced = [];
end

end

function validate_dij(dij)
if ~isstruct(dij) || ~isfield(dij, 'totalNumOfBixels')
    error('cut_dij_by_mask:invalidInput', ...
        'dij.totalNumOfBixels is required.');
end

if ~isscalar(dij.totalNumOfBixels) || ~isnumeric(dij.totalNumOfBixels) || ...
        ~isfinite(dij.totalNumOfBixels) || dij.totalNumOfBixels < 1 || ...
        floor(dij.totalNumOfBixels) ~= dij.totalNumOfBixels
    error('cut_dij_by_mask:invalidInput', ...
        'dij.totalNumOfBixels must be a positive integer scalar.');
end
end

function mask = validate_mask(mask, nBixels)
if ~(islogical(mask) || isnumeric(mask)) || isempty(mask) || ~isvector(mask)
    error('cut_dij_by_mask:invalidMask', ...
        'mask must be a non-empty logical or numeric vector.');
end

mask = reshape(mask, [], 1);
if numel(mask) ~= nBixels
    error('cut_dij_by_mask:dimensionMismatch', ...
        'mask length %d does not match dij.totalNumOfBixels=%d.', ...
        numel(mask), nBixels);
end

if isnumeric(mask)
    if any(~isfinite(mask)) || any(mask ~= 0 & mask ~= 1)
        error('cut_dij_by_mask:invalidMask', ...
            'Numeric mask must contain only 0/1 values.');
    end
    mask = logical(mask);
end
end

function validate_map(map, nBixels)
if ~isstruct(map) || numel(map) ~= nBixels
    error('cut_dij_by_mask:invalidMap', ...
        'map must be a struct array with length %d.', nBixels);
end
end

function validate_weight_vector(w, nBixels)
if ~isnumeric(w) || isempty(w) || ~isvector(w) || any(~isfinite(w(:)))
    error('cut_dij_by_mask:invalidWeight', ...
        'w must be a non-empty finite numeric vector.');
end

if numel(w) ~= nBixels
    error('cut_dij_by_mask:invalidWeight', ...
        'w length %d does not match dij.totalNumOfBixels=%d.', ...
        numel(w), nBixels);
end
end

function fieldReduced = slice_dose_field(fieldValue, keepIdx, nBixels, fieldName)
if iscell(fieldValue)
    fieldReduced = cell(size(fieldValue));
    for s = 1:numel(fieldValue)
        fieldReduced{s} = slice_matrix_columns(fieldValue{s}, keepIdx, nBixels, fieldName, s);
    end
else
    fieldReduced = slice_matrix_columns(fieldValue, keepIdx, nBixels, fieldName, []);
end
end

function matrixReduced = slice_matrix_columns(matrixValue, keepIdx, nBixels, fieldName, scenarioIx)
if isempty(matrixValue)
    matrixReduced = matrixValue;
    return;
end

if ~ismatrix(matrixValue) || size(matrixValue, 2) ~= nBixels
    if isempty(scenarioIx)
        error('cut_dij_by_mask:invalidDoseField', ...
            'dij.%s must have %d columns.', fieldName, nBixels);
    else
        error('cut_dij_by_mask:invalidDoseField', ...
            'dij.%s{%d} must have %d columns.', fieldName, scenarioIx, nBixels);
    end
end

matrixReduced = matrixValue(:, keepIdx);
end

function vectorReduced = slice_vector_field(vectorValue, keepIdx, nBixels, fieldName)
if ~isvector(vectorValue) || numel(vectorValue) ~= nBixels
    error('cut_dij_by_mask:invalidVectorField', ...
        'dij.%s must be a vector with %d entries.', fieldName, nBixels);
end

vectorReduced = vectorValue(keepIdx);
end
