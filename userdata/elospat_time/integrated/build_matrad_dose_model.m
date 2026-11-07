function model = build_matrad_dose_model(dij, cst, pln, constraintPenalty)
% build_matrad_dose_model  Prepare the native matRad objective callbacks.
%
% This mirrors the objective-related preparation in matRad_fluenceOptimization
% without invoking an optimizer. The returned model evaluates the same matRad
% dose objectives. Hard dose constraints are retained through a differentiable
% squared-violation merit term; their raw violations are always reported.

narginchk(3, 4);
if nargin < 4 || isempty(constraintPenalty)
    constraintPenalty = 0;
end
if ~isscalar(constraintPenalty) || ~isfinite(constraintPenalty) || constraintPenalty < 0
    error('build_matrad_dose_model:invalidPenalty', ...
        'constraintPenalty must be a non-negative finite scalar.');
end
if ~isstruct(dij) || ~isfield(dij, 'physicalDose') || ...
        ~isfield(dij, 'totalNumOfBixels')
    error('build_matrad_dose_model:invalidDij', ...
        'dij must be a matRad dose influence structure.');
end
if ~iscell(cst) || size(cst, 2) < 6 || ~isstruct(pln)
    error('build_matrad_dose_model:invalidInput', ...
        'cst and pln must be valid matRad inputs.');
end

cstPrepared = matRad_setOverlapPriorities(cst);
if ~isfield(pln, 'numOfFractions') || isempty(pln.numOfFractions)
    error('build_matrad_dose_model:missingFractions', ...
        'pln.numOfFractions is required.');
end

for i = 1:size(cstPrepared, 1)
    if isstruct(cstPrepared{i, 6})
        cstPrepared{i, 6} = arrayfun( ...
            @matRad_DoseOptimizationFunction.convertOldOptimizationStruct, ...
            cstPrepared{i, 6}, 'UniformOutput', false);
    end
    for j = 1:numel(cstPrepared{i, 6})
        f = cstPrepared{i, 6}{j};
        if ~isa(f, 'matRad_DoseOptimizationFunction')
            f = matRad_DoseOptimizationFunction.createInstanceFromStruct(f);
        end
        f = f.setDoseParameters(f.getDoseParameters() / pln.numOfFractions);
        cstPrepared{i, 6}{j} = f;
    end
end

cstPrepared = matRad_resizeCstToGrid(cstPrepared, ...
    dij.ctGrid.x, dij.ctGrid.y, dij.ctGrid.z, ...
    dij.doseGrid.x, dij.doseGrid.y, dij.doseGrid.z);

if ~isfield(pln, 'propOpt') || isempty(pln.propOpt)
    pln.propOpt = struct();
end
cfg = MatRad_Config.instance();
if ~isfield(pln.propOpt, 'clearUnusedVoxels')
    pln.propOpt.clearUnusedVoxels = cfg.defaults.propOpt.clearUnusedVoxels;
end
dijPrepared = dij;
if pln.propOpt.clearUnusedVoxels
    dijPrepared = matRad_clearUnusedVoxelsFromDij(cstPrepared, dijPrepared);
end

if ~isfield(pln, 'bioModel')
    pln.bioModel = 'none';
end
if ~isa(pln.bioModel, 'matRad_BiologicalModel')
    pln.bioModel = matRad_BiologicalModel.validate(pln.bioModel, pln.radiationMode);
end
if ~isfield(pln.propOpt, 'quantityOpt') || isempty(pln.propOpt.quantityOpt)
    pln.propOpt.quantityOpt = pln.bioModel.defaultReportQuantity;
end

switch pln.propOpt.quantityOpt
    case 'effect'
        backProjection = matRad_EffectProjection;
    case 'RBExDose'
        if isa(pln.bioModel, 'matRad_ConstantRBE') || ...
                (isstruct(pln.bioModel) && strcmp(pln.bioModel.model, 'constRBE'))
            backProjection = matRad_ConstantRBEProjection;
            if ~isfield(dijPrepared, 'RBE')
                dijPrepared.RBE = 1.1;
            end
        else
            backProjection = matRad_VariableRBEProjection;
        end
    case 'BED'
        backProjection = matRad_BEDProjection;
    otherwise
        backProjection = matRad_DoseProjection;
end

totNumCtScen = size(dijPrepared.physicalDose, 1);
if ~isfield(pln, 'multScen')
    pln.multScen = 'nomScen';
end
if ~isa(pln.multScen, 'matRad_ScenarioModel')
    pln.multScen = matRad_ScenarioModel.create(pln.multScen, ...
        struct('numOfCtScen', totNumCtScen));
end
if isfield(pln.propOpt, 'scen4D')
    scen4D = pln.propOpt.scen4D;
else
    scen4D = 1;
end
if isequal(scen4D, 'all')
    scen4D = 1:totNumCtScen;
end

linIx = find(~cellfun(@isempty, dijPrepared.physicalDose(scen4D, :, :)))';
linNominal = find(~cellfun(@isempty, dijPrepared.physicalDose(scen4D, 1, 1)))';
hasProb = false;
hasRobust = false;
for i = 1:size(cstPrepared, 1)
    for j = 1:numel(cstPrepared{i, 6})
        robustness = cstPrepared{i, 6}{j}.robustness;
        hasProb = hasProb || (strcmp(robustness, 'PROB') && numel(linIx) > 1);
        hasRobust = hasRobust || (~strcmp(robustness, 'none') && numel(linIx) > 1);
    end
end
if hasProb
    dijPrepared = matRad_calculateProbabilisticQuantities( ...
        dijPrepared, cstPrepared, pln);
end
if ~hasRobust || hasProb
    ixForOpt = scen4D;
else
    ixForOpt = linIx;
end
backProjection.scenarios = ixForOpt;
backProjection.scenarioProb = pln.multScen.scenProb;
backProjection.nominalCtScenarios = linNominal;

optiProb = matRad_OptimizationProblem(backProjection);
if isfield(pln.propOpt, 'useLogSumExpForRobOpt')
    if pln.propOpt.useLogSumExpForRobOpt
        optiProb.useMaxApprox = 'logsumexp';
    else
        optiProb.useMaxApprox = 'none';
    end
end

[cl, cu] = optiProb.matRad_getConstraintBounds(cstPrepared);
model = struct();
model.optiProb = optiProb;
model.dij = dijPrepared;
model.cst = cstPrepared;
model.pln = pln;
model.constraintLower = reshape(double(cl), [], 1);
model.constraintUpper = reshape(double(cu), [], 1);
model.constraintPenalty = double(constraintPenalty);
model.nWeights = double(dijPrepared.totalNumOfBixels);
model.hasConstraints = ~isempty(cl);
model.source = 'native matRad objective/gradient callbacks';
end
