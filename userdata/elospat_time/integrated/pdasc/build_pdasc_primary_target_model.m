function model = build_pdasc_primary_target_model(dij, cst, pln)
%BUILD_PDASC_PRIMARY_TARGET_MODEL Build Zhao-style U and uniform p0.
%
% Zhao et al. optimized only the primary CTV. For a multi-target matRad cst,
% this comparator defines "primary" deterministically as the TARGET with the
% largest positive prescribed dose inferred from its objective parameters.

narginchk(3, 3);
if iscell(dij.physicalDose), D = dij.physicalDose{1}; else, D = dij.physicalDose; end
targetRows = find(cellfun(@(x) strcmpi(char(string(x)), 'TARGET'), cst(:,3)));
if isempty(targetRows)
    error('build_pdasc_primary_target_model:noTarget', ...
        'cst contains no TARGET structure.');
end
prescriptions = NaN(numel(targetRows), 1);
for j = 1:numel(targetRows)
    prescriptions(j) = infer_prescription(cst{targetRows(j), 6});
end
if any(isfinite(prescriptions) & prescriptions > 0)
    valid = prescriptions; valid(~isfinite(valid) | valid <= 0) = -Inf;
    [prescription, local] = max(valid);
else
    error('build_pdasc_primary_target_model:noPrescription', ...
        'No positive target prescription could be inferred from cst.');
end
row = targetRows(local);
voxels = unique(reshape(double(cst{row,4}{1}), [], 1));
if isempty(voxels), error('build_pdasc_primary_target_model:emptyTarget', ...
        'The selected primary target has no voxels.'); end
fractions = 1;
if isfield(pln, 'numOfFractions') && ~isempty(pln.numOfFractions)
    fractions = double(pln.numOfFractions);
end
prescriptionPerFraction = prescription / fractions;
model = struct('U', sparse(double(D(voxels,:))), ...
    'p0', prescriptionPerFraction * ones(numel(voxels),1), ...
    'targetRow', row, 'targetName', char(string(cst{row,2})), ...
    'targetVoxels', voxels, 'prescriptionTotal', prescription, ...
    'prescriptionPerFraction', prescriptionPerFraction, ...
    'numFractions', fractions, ...
    'selectionRule', 'highest prescribed-dose TARGET (primary target proxy)');
end

function value = infer_prescription(objectives)
value = NaN; candidates = zeros(0,1);
for i = 1:numel(objectives)
    objective = objectives{i}; parameters = [];
    if isstruct(objective) && isfield(objective, 'parameters')
        parameters = objective.parameters;
    elseif isobject(objective)
        try, parameters = objective.getDoseParameters(); catch, end
    end
    if iscell(parameters)
        for j = 1:numel(parameters)
            if isnumeric(parameters{j})
                candidates = [candidates; double(parameters{j}(:))]; %#ok<AGROW>
            end
        end
    elseif isnumeric(parameters)
        candidates = [candidates; double(parameters(:))]; %#ok<AGROW>
    end
end
candidates = candidates(isfinite(candidates) & candidates > 0);
if ~isempty(candidates), value = min(candidates); end
end
