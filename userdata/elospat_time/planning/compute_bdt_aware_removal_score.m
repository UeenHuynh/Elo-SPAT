function removalScore = compute_bdt_aware_removal_score(doseImportance, bdtSaving, epsValue)
% compute_bdt_aware_removal_score  Wrapper removal score for layer/spot pruning.
%
% removalScore = doseImportance ./ (bdtSaving + epsValue)

narginchk(2, 3);

if nargin < 3 || isempty(epsValue)
    epsValue = 1e-12;
end

if ~isnumeric(doseImportance) || ~isnumeric(bdtSaving) || ...
        ~isequal(size(doseImportance), size(bdtSaving))
    error('compute_bdt_aware_removal_score:invalidInput', ...
        'doseImportance and bdtSaving must be numeric arrays of equal size.');
end

if ~isscalar(epsValue) || ~isnumeric(epsValue) || ~isfinite(epsValue) || epsValue <= 0
    error('compute_bdt_aware_removal_score:invalidInput', ...
        'epsValue must be a positive finite scalar.');
end

removalScore = double(doseImportance) ./ (double(bdtSaving) + double(epsValue));
removalScore(~isfinite(removalScore)) = inf;
end
