function output = run_integrated_rms_elospat(dij, cst, pln, denseW, map, groups, cpGroups, cfg, machine)
% run_integrated_rms_elospat  Integrated MATLAB-only RMS-ELO-SPAT pipeline.
%
% Dense Arc -> integrated optimizer -> fixed support -> reduced-dij IPOPT.
% The legacy layer/spot pruning wrapper is not called by this function.

narginchk(9, 9);
effectiveCfg = configure_integrated_stage(cfg);
doseModel = build_matrad_dose_model( ...
    dij, cst, pln, effectiveCfg.integrated.constraintPenalty);
context = build_integrated_context(groups, cpGroups, map, machine, effectiveCfg);
optimizerResult = run_integrated_optimizer( ...
    doseModel, context, denseW, effectiveCfg);
support = extract_integrated_support( ...
    optimizerResult.w, denseW, context, effectiveCfg);

finalCfg = effectiveCfg;
finalCfg.enforceSelectedMinimumWeight = ...
    effectiveCfg.integrated.enforceFinalMinimumWeight;
finalResult = run_final_reoptimization(dij, cst, pln, support.w, ...
    finalCfg, support.mask, support.layerMask);
if any(abs(finalResult.w(~support.mask)) > 0)
    error('run_integrated_rms_elospat:supportLeakage', ...
        'Final IPOPT reoptimization reactivated an excluded spot.');
end

optimizerSeq = extract_delivery_sequence( ...
    support.w, map, groups, effectiveCfg);
finalSeq = extract_delivery_sequence( ...
    finalResult.w, map, groups, effectiveCfg);
output = struct();
output.optimizer = optimizerResult;
output.support = support;
output.final = finalResult;
output.effectiveConfig = effectiveCfg;
output.optimizationTimeSurrogate = optimizerResult.parts.timeSurrogate;
output.selectedSupportBDT = estimate_bdt(optimizerSeq, machine);
output.finalBDT = estimate_bdt(finalSeq, machine);
output.formulationStage = effectiveCfg.integrated.componentStage;
output.runtime = 'MATLAB-only';
end
