function [resultGUI, optimizer] = reoptimize_reduced_problem(dijReduced, cst, pln, wInit)
% reoptimize_reduced_problem  Re-optimize a reduced dij with optional warm-start.
%
% INPUTS
%   dijReduced - reduced matRad dose influence matrix
%   cst        - matRad clinical structure table
%   pln        - matRad plan structure
%   wInit      - optional warm-start weight vector matching reduced dij
%
% OUTPUTS
%   resultGUI  - matRad optimization result structure
%   optimizer  - optimizer object/info returned by matRad

narginchk(3, 4);

validate_reduced_dij(dijReduced);
validate_cst(cst);
validate_pln(pln);

if nargin < 4
    wInit = [];
end

if ~isempty(wInit)
    wInit = validate_warm_start(wInit, double(dijReduced.totalNumOfBixels));
    [resultGUI, optimizer] = matRad_fluenceOptimization(dijReduced, cst, pln, wInit);
else
    [resultGUI, optimizer] = matRad_fluenceOptimization(dijReduced, cst, pln);
end

end

function validate_reduced_dij(dijReduced)
if ~isstruct(dijReduced) || ~isfield(dijReduced, 'totalNumOfBixels')
    error('reoptimize_reduced_problem:invalidInput', ...
        'dijReduced.totalNumOfBixels is required.');
end

if ~isscalar(dijReduced.totalNumOfBixels) || ~isnumeric(dijReduced.totalNumOfBixels) || ...
        ~isfinite(dijReduced.totalNumOfBixels) || dijReduced.totalNumOfBixels < 1 || ...
        floor(dijReduced.totalNumOfBixels) ~= dijReduced.totalNumOfBixels
    error('reoptimize_reduced_problem:invalidInput', ...
        'dijReduced.totalNumOfBixels must be a positive integer scalar.');
end
end

function validate_cst(cst)
if ~iscell(cst) || size(cst, 2) < 6
    error('reoptimize_reduced_problem:invalidInput', ...
        'cst must be a matRad cell array with at least 6 columns.');
end
end

function validate_pln(pln)
if ~isstruct(pln)
    error('reoptimize_reduced_problem:invalidInput', ...
        'pln must be a matRad plan struct.');
end
end

function wInit = validate_warm_start(wInit, nBixels)
if ~isnumeric(wInit) || isempty(wInit) || ~isvector(wInit) || any(~isfinite(wInit(:)))
    error('reoptimize_reduced_problem:invalidWarmStart', ...
        'wInit must be a non-empty finite numeric vector.');
end

wInit = reshape(double(wInit), [], 1);
if numel(wInit) ~= nBixels
    error('reoptimize_reduced_problem:dimensionMismatch', ...
        'wInit length %d does not match dijReduced.totalNumOfBixels=%d.', ...
        numel(wInit), nBixels);
end
end
