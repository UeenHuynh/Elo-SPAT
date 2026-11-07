% test_integrated_native_matrad_objective  Audit native matRad callback use.

testDir = fileparts(mfilename('fullpath'));
workspaceDir = fileparts(testDir);
eloRoot = fileparts(fileparts(workspaceDir));
matRadDir = fullfile(fileparts(eloRoot), 'matRad');
addpath(matRadDir);
matRad_rc;
addpath(fullfile(workspaceDir, 'integrated'));

D = sparse([1.0, 0.2, 0.0; 0.1, 0.8, 0.3; 0.0, 0.2, 1.1]);
dij = struct();
dij.physicalDose = {D};
dij.totalNumOfBixels = 3;
dij.doseGrid.numOfVoxels = 3;

cst = cell(1, 6);
cst{1, 3} = 'TARGET';
cst{1, 4} = {(1:3)'};
cst{1, 5} = struct('alphaX', 0.1, 'betaX', 0.05);
cst{1, 6} = {DoseObjectives.matRad_SquaredDeviation(2.5, 1.0), ...
    DoseConstraints.matRad_MinMaxDose(0, 0.5, 'voxelwise')};

bp = matRad_DoseProjection;
bp.scenarios = 1;
bp.scenarioProb = 1;
bp.nominalCtScenarios = 1;
optiProb = matRad_OptimizationProblem(bp);
model = struct();
model.optiProb = optiProb;
model.dij = dij;
model.cst = cst;
[model.constraintLower, model.constraintUpper] = ...
    optiProb.matRad_getConstraintBounds(cst);
model.constraintPenalty = 4;
model.nWeights = 3;
model.hasConstraints = true;

w = [0.7; 0.4; 0.9];
[value, gradient] = evaluate_matrad_dose_model(model, w, true);
constraintResidual = max(D*w - 0.5, 0);
directValue = 2.5 / 3 * sum((D*w - 1).^2) + ...
    0.5 * model.constraintPenalty * sum(constraintResidual.^2);
directGradient = 2 * 2.5 / 3 * D' * (D*w - 1) + ...
    model.constraintPenalty * D' * constraintResidual;
assert(abs(value - directValue) < 1e-12, ...
    'Native matRad objective value changed in the adapter.');
assert(norm(gradient - directGradient) < 1e-11, ...
    'Native matRad objective gradient changed in the adapter.');
fd = validate_matrad_objective_gradient(model, w, (1:3)', 1e-6);
assert(fd.maximumRelativeError < 1e-6, ...
    'Native matRad objective gradient failed finite differences.');
fprintf('[OK] Native matRad dose objective/gradient adapter validated.\n');
