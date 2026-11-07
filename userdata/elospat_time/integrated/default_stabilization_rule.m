function rule = default_stabilization_rule()
% default_stabilization_rule  Scale-free practical stabilization criteria.
%
% These criteria are deliberately independent of the proximal-gradient
% stationarity diagnostic. They are fixed before the full 20-degree run and
% are not calibrated to make any particular alpha value pass.

rule = struct();
rule.minimumIterations = 100;
rule.trailingWindow = 50;
rule.relativeObjectiveWindowTolerance = 1e-4;
rule.maxPerIterationRelativeObjectiveTolerance = 1e-5;
rule.objectiveIncreaseTolerance = 1e-12;
rule.maximumStepRatio = 1.25;
rule.minimumBDTCheckpoints = 3;
rule.relativeBDTRangeTolerance = 1e-3; % 0.1 percent
rule.constraintTolerance = 1e-6;
rule.supportLeakageTolerance = 1e-10;

% A stationarity threshold is intentionally not calibrated yet. Setting a
% finite value in a future calibration study enables stationary_converged;
% it never changes numerically_stabilized.
rule.calibratedStationarityTolerance = NaN;
end
