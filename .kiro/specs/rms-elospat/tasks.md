# Implementation Plan: RMS-ELO-SPAT

## Overview

This plan implements the RMS-ELO-SPAT research framework as an external MATLAB module on top of matRad. The implementation follows 8 phases: infrastructure, arc geometry + mapping, objective functions, delivery model, reduction, wrapper loop integration, robustness, and evaluation/export. All code resides in `userdata/elospat_time/` without modifying matRad core.

## Tasks

- [ ] 1. Create project folder structure (`userdata/elospat_time/` with subdirectories config/, planning/, mapping/, objective/, delivery/, reduction/, robustness/, evaluation/, tests/) and README.md with research-only disclaimer. Verify matRad is available on path. (Requirements: 17, 18)
- [ ] 2. Implement `config/default_config.m` returning cfg struct with all tunable parameters (mode, patientFile, radiationMode, machine, bioModel, geometry, penalty weights, smoothing params, thresholds, robust settings, wrapper settings, evaluation tolerances, output settings, randomSeed). (Requirements: 17)
- [ ] 3. Implement `config/default_machine_model.m` returning machine struct with energy switching (upBase=1.5s, downBase=0.5s, perMeVUp=0.02, perMeVDown=0.005), spot switching (0.002s), irradiation (muRate=1000), gantry (maxSpeed=6, maxAccel=1, toleranceWindow=0.5), buffer (perLayer=0.05, perCP=0.02), and validation.status='proxy_unvalidated'. (Requirements: 8, 17)
- [ ] 4. Implement `config/validate_config.m` checking all required fields exist, numeric values in valid ranges, patient file exists on disk, arc spacing divides 360. Throws descriptive error on invalid input. (Requirements: 17)
- [ ] 5. Implement `planning/build_impt_plan.m` that configures pln for protons with cfg.imptAngles (default [90,270]), calls matRad_generateStf, matRad_calcDoseInfluence, matRad_fluenceOptimization, checks convergence, returns pln/stf/dij/result. Verify on PROSTATE.mat. (Requirements: 1)
- [ ] 6. Implement `planning/build_discrete_arc_plan.m` generating gantry angles 0:arcSpacing:(360-arcSpacing) with couch=0, calling matRad_generateStf and matRad_calcDoseInfluence. Reject invalid spacing. Verify 20° produces 18 CPs with valid dij. (Requirements: 2)
- [ ] 7. Implement `planning/run_dense_arc_baseline.m` calling matRad_fluenceOptimization on full arc dij. This is the Dense Arc baseline method. Verify convergence on 20° arc. (Requirements: 2, 13)
- [ ] 8. Implement `mapping/map_bixel_to_cp_energy_spot.m` iterating over all bixels using dij.beamNum/rayNum/bixelNum and stf(b).ray(r).energy(s) to build map struct array with globalIx, beamNum, rayNum, bixelNum, gantryAngle, energy, cpIndex, energyIndex. Assign energyIndex from sorted unique energies per beam. (Requirements: 3)
- [ ] 9. Implement `mapping/build_group_index.m` extracting unique (cpIndex, energyIndex) pairs into groups struct array (.ix, .nSpots, .energy, .cpIndex, .energyIndex) and unique cpIndex into cpGroups struct array (.ix, .gantryAngle, .cpIndex). (Requirements: 3)
- [ ] 10. Implement `mapping/validate_mapping.m` asserting numel(map)==dij.totalNumOfBixels, all globalIx are 1:N consecutive, every bixel in exactly one group. Throws error with affected indices on failure. Also implement `mapping/get_group_indices.m` and `mapping/get_control_point_indices.m` accessors. (Requirements: 3)
- [ ] 11. Implement `objective/l21half_group_penalty.m` computing val = sum_g (sum_{i in g} w_i^2 + epsEL)^(1/4) and grad_i = 0.5*w_i*(sum_{j in g} w_j^2 + epsEL)^(-3/4). Create `tests/test_l21half_gradient.m` verifying finite-difference agreement (95% pass at relErr < 1e-4). (Requirements: 4)
- [ ] 12. Implement `objective/control_point_barrier.m` computing val = -sum_b log(sqrt(sum_{i in C_b} w_i^2 + epsB)) and grad_i = -w_i/(sum_{j in C_b} w_j^2 + epsB). Create `tests/test_barrier_gradient.m` with finite-difference validation. (Requirements: 5)
- [ ] 13. Implement `objective/energy_sequencing_penalty.m` computing Phi_ES using soft activity z_tilde = a/(a+epsZ) and switch-up penalty sum over consecutive CP pairs. Return both smooth penalty value and hard switch-up count. Create `tests/test_sequencing_penalty.m` verifying zero penalty for monotonic descending sequences. (Requirements: 6)
- [ ] 14. Implement `objective/spot_sparsity_penalty.m` computing val = sum_i w_i/(w_i+epsS) and grad_i = epsS/(w_i+epsS)^2. Create `tests/test_spot_sparsity_gradient.m` with finite-difference validation. (Requirements: 10)
- [ ] 15. Implement `objective/composite_score.m` computing total score = lambdaEL*R_EL + lambdaB*B + lambdaES*Phi_ES + lambdaS*R_spot + lambdaT*T_estimated, returning score and components struct. Verify all-zero lambdas give score=0. (Requirements: 14, 15)
- [ ] 16. Implement `delivery/extract_delivery_sequence.m` identifying active layers/spots by threshold, ordering CPs by gantry angle ascending, ordering layers within CP by energy descending (default). Build seq struct array with cpIndex, gantryAngle, energy, spotIndices, spotWeights, numActiveSpots, layerWeightSum, prevEnergy, nextEnergy. (Requirements: 7)
- [ ] 17. Implement `delivery/order_layers_within_control_point.m` supporting 'descending' (default), 'ascending', 'custom' modes. (Requirements: 7)
- [ ] 18. Implement `delivery/estimate_bdt.m` computing T_energy (direction-dependent: upBase+perMeVUp*dE for up, downBase+perMeVDown*|dE| for down), T_spot (max(nSpots-1,0)*switchTime per layer), T_irradiation (totalWeight/muRate), T_gantry (trapezoidal profile), T_buffer (nLayers*perLayer + nCP*perCP). Verify T_total == sum of components. (Requirements: 8)
- [ ] 19. Implement `delivery/estimate_gantry_time.m` with trapezoidal motion profile (if theta<=vmax^2/amax: t=2*sqrt(theta/amax), else: t=theta/vmax+vmax/amax) and `delivery/validate_gantry_constraints.m`. Create `tests/test_bdt_estimator.m`. (Requirements: 8)
- [ ] 20. Implement `reduction/build_layer_mask.m` (logical nGroups x 1 based on group L2 norm > threshold) and `reduction/build_spot_mask.m` (logical N x 1 based on w > threshold). (Requirements: 11)
- [ ] 21. Implement `reduction/expand_layer_mask_to_bixel_mask.m` converting group-level layerMask (nGroups) to bixel-level bixelMask (totalNumOfBixels) by setting bixelMask(groups(g).ix)=true for active groups. Create `tests/test_expand_layer_mask.m` and `tests/test_mask_dimension_safety.m` verifying layerMask cannot be accidentally used as bixelMask. (Requirements: 11)
- [ ] 22. Implement `reduction/cut_dij_by_mask.m` handling all dose-influence fields (physicalDose, mAlphaDose, mSqrtBetaDose, mLETDose as cell arrays per scenario; beamNum, rayNum, bixelNum, minMU, maxMU, numParticlesPerMU as vectors). Assert mask length == dij.totalNumOfBixels. Skip missing fields safely. Preserve reverse mapping. Create `tests/test_dij_reduction.m` including multi-scenario test. (Requirements: 11)
- [ ] 23. Implement `reduction/cut_stf_by_mask.m` removing inactive rays/bixels from stf, adjusting ray/bixel counts per beam, removing empty rays. Must be consistent with cut_dij_by_mask output. (Requirements: 11)
- [ ] 24. Implement `reduction/reoptimize_reduced_problem.m` calling matRad_fluenceOptimization with reduced dij and warm-start wInit. Verify dimensions match before calling. (Requirements: 11, 15)
- [ ] 25. Implement `planning/run_elospat_like_stage.m` with iterative layer filtering: compute doseContribution(g)=norm(D(:,ix_g)*w(ix_g),2), compute bdtSaving(g)=T_current-T_without_g, compute removalScore=doseImportance/(bdtSaving+eps), remove bottom-K groups, expand to bixel mask, cut dij, re-optimize, check convergence. Document as ELO-SPAT-inspired wrapper approximation. Create `tests/test_bdt_aware_removal_score.m`. (Requirements: 4, 6, 8, 15)
- [ ] 26. Implement `planning/run_spot_sparse_stage.m` with iterative spot filtering: compute per-spot dose importance and BDT saving, removal score, identify sub-MMU + low-score candidates, cut dij, re-optimize, check CTV D95 degradation (revert if > tolerance). (Requirements: 10, 15)
- [ ] 27. Implement `planning/run_final_reoptimization.m` optionally configuring robust scenarios (COWC/OWC) on pln/cst before calling matRad_fluenceOptimization. Support both nominal and robust modes. (Requirements: 12)
- [ ] 28. Implement `robustness/configure_robust_scenarios.m` (set pln.multScen with configured uncertainty) and `robustness/configure_robust_objectives.m` (set cst robustness mode). (Requirements: 12)
- [ ] 29. Implement `robustness/evaluate_robust_scenarios.m` (accepts ct, cst, stf, dij, pln, w, cfg) computing per-scenario dose, DVH metrics, worst-case identification, scenario pass rate, std dose cube. Distinguish robust optimization vs evaluation-only. Implement `robustness/compute_worst_case_metrics.m`. (Requirements: 12)
- [ ] 30. Implement `evaluation/compute_dvh_metrics.m` computing D98, D95, D5, V95, Dmean, Dmax, homogeneity index (D5-D95)/Dprescribed, DVH curves per structure at 0.1 Gy resolution. (Requirements: 14)
- [ ] 31. Implement `evaluation/compute_delivery_metrics.m` computing nControlPoints, nActiveLayers, nActiveSpots, nSwitchUp, totalSwitchMagnitude, full BDT decomposition, spot weight statistics (min, median, pctSubMMU). (Requirements: 14)
- [ ] 32. Implement `evaluation/compute_solver_metrics.m` extracting final objective, dose component, regularization component, iterations, runtime, convergence status from matRad result. (Requirements: 14)
- [ ] 33. Implement `evaluation/compare_methods_table.m` building MATLAB table with one row per method (IMPT, Dense Arc, ELO-SPAT-like, Post-filter, RMS-ELO-SPAT) and columns for all dose/delivery/solver metrics. (Requirements: 13, 14)
- [ ] 34. Implement `evaluation/export_latex_table.m` generating LaTeX tabular with columns: Method, CTV D95, CTV D98, OAR Dmax, Active EL, Active Spots, Switch-up, T_energy, T_spot, T_gantry, T_total. (Requirements: 14)
- [ ] 35. Implement `evaluation/export_tradeoff_data.m` saving dose-vs-estimated-BDT pairs per method for Pareto front visualization. (Requirements: 14)
- [ ] 36. Implement `planning/run_post_filter_baseline.m` applying threshold-based spot removal on dense arc result, cutting dij, re-optimizing. This is Method 4 comparison baseline. (Requirements: 13)
- [ ] 37. Implement `run_rms_elospat_main.m` orchestrator integrating all modules: setup, IMPT baseline, arc setup, dense arc baseline, mapping, ELO-SPAT-like stage, layer reduction (with expand_layer_mask_to_bixel_mask), re-opt, spot sparse stage, spot reduction, final re-opt, post-filter baseline, robust eval, metrics, export. Support cfgOverride for sweeps. Include research-only disclaimer. (Requirements: 1, 2, 3, 13, 15, 17, 18)
- [ ] 38. Run end-to-end integration test on PROSTATE.mat with 20° arc (18 CPs). Verify all 5 methods produce valid results, RMS-ELO-SPAT estimated BDT < Dense Arc estimated BDT, CTV D95 within 2% tolerance, comparison table and LaTeX export generated. Document runtime. (Requirements: 1, 2, 13, 14, 15)
- [ ] 39. Run robustness integration test with cfg.robust.enabled=true, mode='COWC', setupMM=[3,3,3], rangePercent=3.5. Verify worst-case CTV D95 computed and within tolerance, scenario pass rate reported. (Requirements: 12)
- [ ] 40. Add parameter sweep logic to orchestrator: iterate over Cartesian product of specified grid arrays, save each result to separate .mat file, log failures and continue remaining combinations. Verify on small 2x2 grid. (Requirements: 17)
- [ ] 41. Run validation with finer arc spacing: 10° (36 CPs) and 5° (72 CPs). Compare estimated BDT reduction across spacing levels. Document memory requirements. Report if memory insufficient for user to address. (Requirements: 2, 14)
- [ ] 42. Final documentation pass: ensure all outputs include research-only disclaimer, machine model labeled as estimated proxy, README.md updated with full usage guide and citation placeholder. (Requirements: 18)

## Task Dependency Graph

```json
{
  "waves": [
    [1],
    [2, 3, 4],
    [5],
    [6, 7],
    [8, 9, 10],
    [11, 12, 13, 14, 15],
    [16, 17, 18, 19],
    [20, 21],
    [22, 23, 24],
    [25, 26],
    [27, 36],
    [28, 29],
    [30, 31, 32],
    [33, 34, 35],
    [37],
    [38, 39, 40, 41, 42]
  ]
}
```

## Notes

- matRad must be cloned/available before Task 5. If not on path, Task 1 will report to user for manual download.
- PROSTATE.mat must be available in matRad's data path. If missing, report to user.
- Tasks 1-4 (infrastructure) can be implemented without MATLAB execution.
- Tasks 5-7 require MATLAB + matRad to verify.
- Memory may be insufficient for 5° or finer arc spacing (Task 41); document and report if so.
- Integrated Mode (custom matRad_OptimizationProblemELOSPAT class) is future work, not included in this task list.
