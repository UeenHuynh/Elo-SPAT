# RMS-ELO-SPAT v2.0 Technical Specification

**Project name:** RMS-ELO-SPAT  
**Full name:** Robust Machine-specific Sparse ELO-SPAT for Discrete Proton Arc Therapy  
**Document type:** Research Engineering Requirements Specification  
**Target platform:** MATLAB + matRad  
**Implementation location:** `userdata/elospat_time/`  
**Status:** Draft v2.0, mathematically and implementation-aligned  
**Primary research objective:** Reduce estimated machine-specific beam delivery time while preserving nominal and robust target/OAR dose quality.

---

## 1. Purpose and Scope

RMS-ELO-SPAT is a matRad-based research sandbox for discrete proton arc therapy (PAT) optimization. It reproduces the main mathematical idea of ELO-SPAT, namely \(L_{2,1/2}\) energy-layer group sparsification, and extends it toward a robust, machine-specific, spot-sparse and deliverability-aware optimization framework.

The system is designed to test the following research hypothesis:

> Extending ELO-SPAT from energy-layer sparsification and energy-sequencing regularization to robust machine-specific sparse optimization can reduce estimated total beam delivery time while preserving nominal and worst-case target coverage, OAR sparing, and deliverability constraints.

This project is not intended for clinical treatment planning. matRad is used as a dose-calculation and inverse-planning research sandbox. The machine-specific delivery model is initially an estimated configurable proxy unless replaced by validated vendor timing data or measured delivery logs.

---

## 2. Literature Alignment

The specification is aligned with the following literature-derived observations:

1. ELO-SPAT already addresses energy-layer group sparsity using an \(L_{2,1/2}\) penalty and includes an energy-sequencing penalty.
2. ELO-SPAT does not directly optimize full machine-specific beam delivery time.
3. Energy-layer sparsity is not equivalent to spot sparsity.
4. Minimum monitor-unit constraints require zero-or-above-threshold logic, not a simple positive lower bound on all spots.
5. Robustness under setup and range uncertainty should be embedded in optimization or at least evaluated across uncertainty scenarios.
6. Machine deliverability depends on energy switching, spot switching, irradiation time, gantry motion, and buffer/tolerance behavior.

The proposed research gap is therefore:

> Develop and evaluate a robust machine-specific sparse extension of ELO-SPAT that jointly considers energy-layer sparsity, spot sparsity, estimated machine-specific delivery time, MMU-aware re-optimization, and robust scenario performance.

---

## 3. Implementation Philosophy

The system SHALL support two implementation modes.

### 3.1 Mode A: Wrapper Prototype Mode

Mode A is the required minimum viable implementation.

In this mode, matRad core code SHALL NOT be modified. The RMS-ELO-SPAT module SHALL operate as an external wrapper around matRad by:

1. generating discrete proton arc geometries;
2. running standard matRad fluence optimization;
3. mapping bixels to control point, energy layer, and spot identities;
4. computing ELO-SPAT-like sparsity, sequencing, and BDT scores;
5. filtering inactive or inefficient energy layers/spots;
6. reducing the dose-influence matrix;
7. re-optimizing with matRad;
8. performing robust re-optimization or robust evaluation;
9. exporting dose, robustness, complexity, and delivery-time metrics.

In Mode A, delivery-time and sequencing terms MAY be used as outer-loop selection criteria rather than fully differentiable IPOPT objective terms.

### 3.2 Mode B: Integrated Optimization Mode

Mode B is the advanced implementation target.

In this mode, the system SHALL implement a custom matRad optimization problem class, for example:

```text
matRad/optimization/@matRad_OptimizationProblemELOSPAT/
```

The custom problem SHALL include differentiable regularization terms directly in the IPOPT objective:

\[
F_{\mathrm{total}}(\mathbf{x})
=
F_{\mathrm{robust}}(\mathbf{x})
+
\lambda_{EL}R_{2,1/2}(\mathbf{x})
+
\lambda_B B(\mathbf{x})
+
\lambda_{ES}\Phi_{ES}(\mathbf{x})
+
\lambda_S R_{\mathrm{spot}}(\mathbf{x})
+
\lambda_T T_{\mathrm{machine}}^{\mathrm{smooth}}(\mathbf{x})
\]

Mode B SHALL only be enabled after analytical or finite-difference-validated gradients are available for all active objective components.

---

## 4. Definitions and Glossary

| Term | Definition |
|---|---|
| matRad | Open-source MATLAB-based radiation treatment planning toolkit used here as a research sandbox |
| PAT | Proton Arc Therapy, approximated here as many discrete proton beam/control-point angles |
| ELO-SPAT | Energy Layer Optimization for Spot-scanning Proton Arc Therapy |
| RMS-ELO-SPAT | Robust Machine-specific Sparse ELO-SPAT, the proposed extension |
| BDT | Beam Delivery Time |
| ELST | Energy Layer Switching Time |
| Bixel | A single optimized beamlet/spot weight in matRad |
| Control point | A discrete gantry-angle position in the arc |
| Energy layer | A proton energy level at a control point |
| Spot | A pencil-beam position within an energy layer |
| Group | A set of bixels sharing the same control point and energy layer |
| \(L_{2,1/2}\) | Nonconvex group-sparsity penalty \(\sum_g \|\mathbf{x}_g\|_2^{1/2}\) |
| MMU | Minimum monitor unit or minimum deliverable spot-weight threshold |
| COWC | Composite Objective Worst Case robust optimization paradigm |
| OWC | Objective-Wise Worst Case robust optimization paradigm |
| DVH | Dose-volume histogram |
| `dij` | matRad dose-influence matrix structure |
| `stf` | matRad steering information structure |
| `pln` | matRad treatment plan structure |
| `cst` | matRad structure/objective table |
| Machine model | Configurable delivery-time parameter structure |
| Wrapper mode | External iterative filtering and re-optimization mode |
| Integrated mode | Custom objective/gradient mode inside matRad optimizer architecture |

---

## 5. Mathematical Formulation

### 5.1 Decision Variables

Let

\[
x_{b,e,s} \geq 0
\]

denote the spot weight for spot \(s\) in energy layer \(e\) at control point \(b\).

Let

\[
\mathbf{x}_{b,e}
=
(x_{b,e,1},x_{b,e,2},\ldots,x_{b,e,|\mathcal{S}_{b,e}|})^\top
\]

be the group of all spot weights in energy layer \(e\) at control point \(b\).

The full spot-weight vector is

\[
\mathbf{x}
=
[\mathbf{x}_{1,1}^{\top},\ldots,\mathbf{x}_{b,e}^{\top},\ldots]^\top.
\]

### 5.2 Dose Model

Let \(\mathbf{A}\) denote the dose-influence matrix. The nominal dose is

\[
\mathbf{d}(\mathbf{x})=\mathbf{A}\mathbf{x}.
\]

For robust optimization, each uncertainty scenario \(\omega \in \Omega\) may have a scenario-specific influence matrix \(\mathbf{A}^{(\omega)}\), giving

\[
\mathbf{d}^{(\omega)}(\mathbf{x})=
\mathbf{A}^{(\omega)}\mathbf{x}.
\]

### 5.3 ELO-SPAT-like Energy-Layer Group Sparsity

The energy-layer group-sparsity penalty SHALL be

\[
R_{2,1/2}(\mathbf{x})
=
\sum_{b=1}^{B}
\sum_{e\in\mathcal{E}_b}
\|\mathbf{x}_{b,e}\|_2^{1/2}.
\]

For numerical smoothing, the implemented group penalty SHALL be

\[
r_g(\mathbf{x})
=
\left(
\sum_{i\in g}x_i^2+\epsilon_{EL}
\right)^{1/4}.
\]

The total implemented penalty SHALL be

\[
R_{EL}(\mathbf{x})
=
\sum_g r_g(\mathbf{x}).
\]

The gradient for \(i\in g\) SHALL be

\[
\frac{\partial r_g}{\partial x_i}
=
\frac{1}{2}
x_i
\left(
\sum_{j\in g}x_j^2+\epsilon_{EL}
\right)^{-3/4}.
\]

### 5.4 Control-Point Activity Preservation Barrier

The log-barrier term is NOT a sparsity penalty. Its role is to discourage empty control points after energy-layer sparsification.

For control point \(b\), define

\[
u_b(\mathbf{x})
=
\left(
\sum_{i\in C_b}x_i^2+\epsilon_B
\right)^{1/2},
\]

where \(C_b\) is the set of all bixels at control point \(b\).

The barrier term SHALL be

\[
B(\mathbf{x})
=
-\sum_{b=1}^{B}\log(u_b(\mathbf{x})).
\]

Its gradient for \(i\in C_b\) SHALL be

\[
\frac{\partial B}{\partial x_i}
=
-\frac{x_i}
{\sum_{j\in C_b}x_j^2+\epsilon_B}.
\]

### 5.5 Energy-Sequencing Penalty

For energy layer \(e\) at control point \(b\), define activity

\[
a_{b,e}(\mathbf{x})
=
\|\mathbf{x}_{b,e}\|_2.
\]

A smooth activity approximation SHALL be

\[
\tilde{z}_{b,e}
=
\frac{a_{b,e}}
{a_{b,e}+\epsilon_Z}.
\]

The switch-up penalty SHALL be

\[
\Phi_{ES}(\mathbf{x})
=
\sum_{b=1}^{B-1}
\sum_{e\in\mathcal{E}_b}
\sum_{e'\in\mathcal{E}_{b+1}}
\tilde{z}_{b,e}\tilde{z}_{b+1,e'}
\left[
E_{b+1,e'}-E_{b,e}
\right]_+.
\]

Here,

\[
[y]_+=\max(y,0).
\]

In Wrapper Mode, \(\Phi_{ES}\) MAY be used as an outer-loop score. In Integrated Mode, \(\Phi_{ES}\) SHALL require analytical or validated numerical gradients.

### 5.6 Spot Sparsity

The smooth spot-sparsity penalty SHALL be

\[
R_{\mathrm{spot}}(\mathbf{x})
=
\sum_i
\frac{x_i}{x_i+\epsilon_S}.
\]

The gradient SHALL be

\[
\frac{\partial R_{\mathrm{spot}}}{\partial x_i}
=
\frac{\epsilon_S}{(x_i+\epsilon_S)^2}.
\]

This term is a smooth approximation encouraging individual spot elimination. It is not sufficient alone to enforce MMU deliverability; thresholding and re-optimization are still required.

### 5.7 MMU Deliverability Condition

The intended deliverability condition is

\[
x_i = 0
\quad\text{or}\quad
x_i \geq x_{\min}.
\]

The system SHALL NOT enforce a simple lower bound \(x_i \geq x_{\min}\) on all candidate spots, because that would incorrectly force every spot to be active.

Instead, the system SHALL apply:

1. smooth spot sparsity during optimization or outer-loop scoring;
2. thresholding of sub-MMU spots;
3. removal of sub-MMU spots from the candidate set;
4. dose re-optimization on the reduced candidate set;
5. fallback if dose degradation exceeds tolerance.

### 5.8 Estimated Machine-Specific BDT Model

The estimated delivery time SHALL be decomposed as

\[
T_{\mathrm{machine}}
=
T_{\mathrm{energy}}
+
T_{\mathrm{spot}}
+
T_{\mathrm{irradiation}}
+
T_{\mathrm{gantry}}
+
T_{\mathrm{buffer}}.
\]

This is a configurable estimated BDT model. It SHALL be called "true BDT" only if machine timing parameters are validated against machine logs, vendor timing specifications, or measured delivery.

The default implementation SHALL use a conservative serial-time approximation. Future implementations MAY replace this with a synchronized delivery timeline where gantry motion, irradiation, and switching partially overlap subject to machine constraints.

### 5.9 Robust Dose Objective

The robust dose objective MAY be implemented as worst-case dose cost:

\[
F_{\mathrm{robust}}(\mathbf{x})
=
\max_{\omega\in\Omega}
F_{\mathrm{dose}}^{(\omega)}(\mathbf{x}).
\]

Alternatively, it MAY use a matRad-supported robust paradigm such as COWC or OWC.

---

## 6. System Architecture

### 6.1 Folder Structure

The implementation SHALL use the following project structure:

```text
userdata/elospat_time/
  run_rms_elospat_main.m
  README.md

  config/
    default_config.m
    default_machine_model.m
    validate_config.m

  planning/
    build_impt_plan.m
    build_discrete_arc_plan.m
    run_dense_arc_baseline.m
    run_elospat_like_stage.m
    run_spot_sparse_stage.m
    run_final_reoptimization.m

  mapping/
    map_bixel_to_cp_energy_spot.m
    build_group_index.m
    validate_mapping.m
    get_group_indices.m
    get_control_point_indices.m

  objective/
    l21half_group_penalty.m
    l21half_group_gradient.m
    control_point_barrier.m
    control_point_barrier_gradient.m
    energy_sequencing_penalty.m
    spot_sparsity_penalty.m
    spot_sparsity_gradient.m
    composite_score.m

  delivery/
    extract_delivery_sequence.m
    order_layers_within_control_point.m
    estimate_bdt.m
    estimate_energy_switch_time.m
    estimate_spot_switch_time.m
    estimate_irradiation_time.m
    estimate_gantry_time.m
    validate_gantry_constraints.m

  reduction/
    build_layer_mask.m
    build_spot_mask.m
    cut_dij_by_mask.m
    cut_stf_by_mask.m
    reoptimize_reduced_problem.m

  robustness/
    configure_robust_scenarios.m
    configure_robust_objectives.m
    evaluate_robust_scenarios.m
    compute_worst_case_metrics.m

  evaluation/
    compute_dvh_metrics.m
    compute_delivery_metrics.m
    compute_solver_metrics.m
    compare_methods_table.m
    export_latex_table.m
    export_tradeoff_data.m

  tests/
    test_mapping_integrity.m
    test_l21half_gradient.m
    test_spot_sparsity_gradient.m
    test_bdt_estimator.m
    test_dij_reduction.m
```

### 6.2 Core Data Flow

```text
Patient data
  -> matRad plan setup
  -> stf generation
  -> dij dose influence calculation
  -> dense arc baseline optimization
  -> bixel-to-group mapping
  -> layer sparsification
  -> layer reduction
  -> re-optimization
  -> spot sparsity / MMU filtering
  -> spot reduction
  -> final nominal or robust re-optimization
  -> robust evaluation
  -> dose/delivery/robustness/solver metric export
```

---

## 7. Requirements

### Requirement 1: Baseline IMPT Planning

**User story:** As a researcher, I want to run a standard multi-beam IMPT plan in matRad so that I have a reference baseline for comparison with arc-based methods.

#### Acceptance Criteria

1. The system SHALL support a baseline IMPT configuration using 2 to 5 user-defined gantry angles.
2. The default IMPT gantry angles SHALL be \([90^\circ,270^\circ]\).
3. The system SHALL support `PROSTATE.mat` as the primary initial test dataset.
4. The system SHALL call `matRad_generateStf`, `matRad_calcDoseInfluence` or equivalent particle dose calculation, and `matRad_fluenceOptimization`.
5. A plan SHALL be considered converged only when matRad/IPOPT reports solved or solved to acceptable level.
6. If the optimizer reaches maximum iterations, the plan SHALL be stored but marked as non-converged.
7. If optimization fails due to infeasibility or numerical error, the system SHALL return the failure status, final available objective value if available, iteration count if available, and error message.
8. The output SHALL include optimized weights, dose cube, DVH metrics, solver status, runtime, and configuration.

---

### Requirement 2: Discrete Proton Arc Surrogate

**User story:** As a researcher, I want to create a discrete proton arc geometry using many gantry angles so that arc delivery can be approximated with matRad's multi-beam infrastructure.

#### Acceptance Criteria

1. The system SHALL generate coplanar control points using

\[
\theta_b = 0:\Delta\theta:(360-\Delta\theta).
\]

2. Supported default angular spacings SHALL include \(20^\circ,10^\circ,5^\circ,4^\circ,3^\circ,2^\circ\).
3. The number of control points SHALL be

\[
B = 360/\Delta\theta
\]

when \(\Delta\theta\) divides \(360^\circ\).
4. All couch angles SHALL default to \(0^\circ\).
5. The system SHALL reject angular spacing values that do not divide \(360^\circ\), unless custom nonuniform gantry angles are explicitly provided.
6. The system SHALL compute the dose-influence matrix for all control points using a single matRad dose calculation workflow where memory permits.
7. If memory limits are exceeded, the system SHALL report the failure and recommend a coarser angle spacing or dose grid.

---

### Requirement 3: Bixel-to-Control-Point/Energy/Spot Mapping

**User story:** As a researcher, I want every bixel mapped to its control point, energy layer, and spot identity so that group sparsity can be applied at the energy-layer level.

#### Acceptance Criteria

1. The mapping length SHALL equal `dij.totalNumOfBixels`.
2. Each global bixel index SHALL appear exactly once.
3. Each bixel SHALL be mapped to:
   - global index;
   - beam/control-point index;
   - gantry angle;
   - ray index;
   - bixel/spot index;
   - physical energy;
   - energy-layer index within the control point.
4. The system SHALL build group indices for every unique \((b,e)\) pair.
5. The system SHALL provide accessor functions for:
   - all bixels in a given energy-layer group;
   - all bixels in a given control point;
   - all active groups above threshold.
6. If a requested group does not exist, the accessor SHALL return an empty index set without error.
7. If any bixel is unmapped or multiply mapped, validation SHALL fail and report the affected indices.
8. The mapping SHALL be recomputed whenever `stf`, `dij`, or gantry geometry changes.

---

### Requirement 4: ELO-SPAT-like \(L_{2,1/2}\) Group Sparsity

**User story:** As a researcher, I want to apply an \(L_{2,1/2}\) group sparsity penalty over energy-layer groups so that inactive layers are driven toward zero.

#### Acceptance Criteria

1. The group penalty SHALL be computed as

\[
R_{EL}(\mathbf{x})
=
\sum_g
\left(
\sum_{i\in g}x_i^2+\epsilon_{EL}
\right)^{1/4}.
\]

2. The gradient for \(i\in g\) SHALL be

\[
\frac{\partial R_{EL}}{\partial x_i}
=
\frac{1}{2}x_i
\left(
\sum_{j\in g}x_j^2+\epsilon_{EL}
\right)^{-3/4}.
\]

3. The default \(\epsilon_{EL}\) SHALL be \(10^{-8}\), with valid range \(10^{-10}\) to \(10^{-6}\).
4. The default \(\lambda_{EL}\) SHALL be 0, with recommended sweep range \(0\) to \(10^{-2}\).
5. A group SHALL be classified as inactive for delivery evaluation if

\[
\|\mathbf{x}_{b,e}\|_2 < \tau_{EL}.
\]

6. The default layer activity threshold \(\tau_{EL}\) SHALL be \(10^{-5}\), with recommended range \(10^{-6}\) to \(10^{-4}\).
7. Integrated Mode SHALL require finite-difference validation of the group-sparsity gradient before use in IPOPT.

---

### Requirement 5: Control-Point Activity Preservation Barrier

**User story:** As a researcher, I want to discourage empty control points so that layer sparsification does not remove all dose-delivering degrees of freedom at a gantry angle.

#### Acceptance Criteria

1. The barrier term SHALL be

\[
B(\mathbf{x})
=
-\sum_b
\log
\left(
\sqrt{
\sum_{i\in C_b}x_i^2+\epsilon_B
}
\right).
\]

2. The gradient for \(i\in C_b\) SHALL be

\[
\frac{\partial B}{\partial x_i}
=
-\frac{x_i}
{\sum_{j\in C_b}x_j^2+\epsilon_B}.
\]

3. The barrier SHALL NOT be described as a sparsity-inducing term.
4. The barrier SHALL be optional.
5. The default \(\lambda_B\) SHALL be 0 in RMS-ELO-SPAT experiments unless reproducing ELO-SPAT-like behavior requires it.
6. If enabled, \(\epsilon_B\) SHALL default to \(10^{-8}\).
7. The system SHALL report the number of empty or near-empty control points before and after applying the barrier.

---

### Requirement 6: Energy-Sequencing Penalty

**User story:** As a researcher, I want to penalize unfavorable low-to-high energy transitions so that plans avoid excessive energy switch-up behavior.

#### Acceptance Criteria

1. The system SHALL compute layer activity as

\[
a_{b,e}=\|\mathbf{x}_{b,e}\|_2.
\]

2. The system SHALL compute soft activity as

\[
\tilde{z}_{b,e}=
\frac{a_{b,e}}
{a_{b,e}+\epsilon_Z}.
\]

3. The energy-sequencing penalty SHALL be

\[
\Phi_{ES}
=
\sum_{b=1}^{B-1}
\sum_{e,e'}
\tilde{z}_{b,e}
\tilde{z}_{b+1,e'}
[E_{b+1,e'}-E_{b,e}]_+.
\]

4. The default \(\epsilon_Z\) SHALL be \(10^{-6}\).
5. The default \(\lambda_{ES}\) SHALL be 0.
6. In Wrapper Mode, \(\Phi_{ES}\) MAY be used as an outer-loop ranking or filtering score.
7. In Integrated Mode, \(\Phi_{ES}\) SHALL require analytical or finite-difference-validated gradients.
8. If no active low-to-high transition exists, the hard switch-up count SHALL be zero.
9. The system SHALL report both:
   - smooth sequencing penalty;
   - hard switch-up event count.

---

### Requirement 7: Delivery Sequence Definition

**User story:** As a researcher, I want a deterministic delivery sequence definition so that BDT estimates are reproducible.

#### Acceptance Criteria

1. The system SHALL extract active layers from the optimized weights using the configured layer threshold.
2. The system SHALL extract active spots using the configured spot threshold.
3. Control points SHALL be ordered by increasing gantry angle unless a custom order is provided.
4. If multiple active energy layers exist at the same control point, the default order SHALL be descending energy.
5. The system SHALL allow user-defined within-control-point energy ordering as an optional override.
6. The system SHALL record the final delivery sequence used for BDT estimation.
7. The sequence record SHALL include, for each delivered layer:
   - control-point index;
   - gantry angle;
   - energy;
   - active spot indices;
   - active spot weights;
   - number of active spots;
   - layer total weight;
   - preceding and following energy values.

---

### Requirement 8: Estimated Machine-Specific BDT Estimator

**User story:** As a researcher, I want an estimated machine-specific delivery-time model that accounts for energy switching, spot switching, irradiation, gantry motion, and buffer times.

#### Acceptance Criteria

1. The BDT estimator SHALL compute

\[
T_{\mathrm{machine}}
=
T_{\mathrm{energy}}
+
T_{\mathrm{spot}}
+
T_{\mathrm{irradiation}}
+
T_{\mathrm{gantry}}
+
T_{\mathrm{buffer}}.
\]

2. The default implementation SHALL use a serial-time approximation.
3. The system SHALL explicitly label this as an estimated BDT proxy unless validated machine-specific timing data are supplied.
4. Energy switching time SHALL be direction-dependent:
   - upward switch: \(t_{\mathrm{up,base}} + k_{\mathrm{up}}|\Delta E|\);
   - downward switch: \(t_{\mathrm{down,base}} + k_{\mathrm{down}}|\Delta E|\);
   - same energy: 0.
5. Default timing values SHALL be configurable and stored in `machine`.
6. Spot switching time SHALL be computed as

\[
T_{\mathrm{spot}}
=
\sum_{\ell}
\max(n_{\ell}-1,0)t_{\mathrm{spot}},
\]

where \(n_{\ell}\) is the number of active spots in delivered layer \(\ell\).
7. Irradiation time SHALL use MU-equivalent spot weights.
8. If matRad weights are not in MU-equivalent units and no conversion data are available, the system SHALL report irradiation time as unavailable or approximate, depending on configuration.
9. Gantry time SHALL use a trapezoidal motion profile:

If \(\theta \leq v_{\max}^2/a_{\max}\),

\[
t=2\sqrt{\theta/a_{\max}}.
\]

Otherwise,

\[
t=\theta/v_{\max}+v_{\max}/a_{\max}.
\]

10. If two consecutive control points are within the gantry tolerance window, gantry transition time SHALL be set to zero.
11. Buffer time SHALL be

\[
T_{\mathrm{buffer}}
=
N_{\mathrm{layers}}t_{\mathrm{layer-buffer}}
+
N_{\mathrm{CP}}t_{\mathrm{CP-buffer}}.
\]

12. The system SHALL report all BDT components separately and not only total BDT.
13. The system SHALL support replacement of the default estimator by a user-provided machine-log-derived timing model.

---

### Requirement 9: Smooth BDT Objective Approximation

**User story:** As a researcher, I want a differentiable BDT approximation for Integrated Mode so that delivery time can be optimized directly.

#### Acceptance Criteria

1. Smooth BDT objective terms SHALL only be enabled in Integrated Mode.
2. Active/inactive decisions SHALL use sigmoid or softplus-based approximations.
3. The smooth active spot approximation MAY be

\[
\hat{z}_i
=
\sigma(k(x_i-\tau_S)).
\]

4. The smooth active layer approximation MAY be

\[
\hat{z}_{b,e}
=
\sigma(k(\|\mathbf{x}_{b,e}\|_2-\tau_{EL})).
\]

5. The steepness parameter \(k\) SHALL be configurable.
6. The smooth BDT objective SHALL be validated against hard-threshold BDT evaluation after optimization.
7. The system SHALL report the discrepancy between smooth BDT and hard-threshold BDT.
8. Integrated optimization SHALL fail validation if smooth and hard BDT rankings are inconsistent across tested hyperparameter settings.

---

### Requirement 10: Spot Sparsity and MMU-Aware Re-optimization

**User story:** As a researcher, I want to penalize low-weight spots and enforce deliverability through thresholding and re-optimization.

#### Acceptance Criteria

1. The spot sparsity penalty SHALL be

\[
R_{\mathrm{spot}}(\mathbf{x})
=
\sum_i
\frac{x_i}{x_i+\epsilon_S}.
\]

2. The gradient SHALL be

\[
\frac{\partial R_{\mathrm{spot}}}{\partial x_i}
=
\frac{\epsilon_S}{(x_i+\epsilon_S)^2}.
\]

3. The default \(\epsilon_S\) SHALL be \(10^{-3}\).
4. The default \(\lambda_S\) SHALL be 0.
5. The MMU threshold \(x_{\min}\) SHALL be configurable.
6. The system SHALL identify all spots satisfying

\[
0 < x_i < x_{\min}.
\]

7. Sub-MMU spots SHALL be fixed to zero and removed from the candidate set before re-optimization.
8. The system SHALL re-optimize on the reduced `dij`/`stf` candidate set.
9. If CTV D95 degrades by more than the configured tolerance relative to the pre-threshold plan, the system SHALL warn and revert to the pre-threshold plan unless override is enabled.
10. If all spots are removed, the system SHALL return an error and preserve the previous feasible solution.
11. The system SHALL report:
    - number of spots before filtering;
    - number of removed sub-MMU spots;
    - number of active spots after re-optimization;
    - final minimum, median, and percentile distribution of active spot weights.

---

### Requirement 11: Dose-Influence Matrix and STF Reduction

**User story:** As a researcher, I want to remove inactive layers and spots from `dij` and `stf` safely so that reduced optimization problems remain consistent.

#### Acceptance Criteria

1. The system SHALL support reduction by:
   - energy-layer mask;
   - spot mask;
   - combined mask.
2. The reduced `dij` SHALL preserve correct correspondence between columns and remaining bixels.
3. The reduced `stf` SHALL preserve correct beam, ray, energy, and spot structure.
4. The system SHALL validate that the reduced weight vector length equals reduced `dij.totalNumOfBixels`.
5. The system SHALL validate that forward dose calculation with the reduced structure can be executed.
6. The system SHALL preserve a mapping from reduced bixel indices back to original bixel indices.
7. The system SHALL fail safely if a mask removes all bixels from any required control point, unless empty control points are explicitly allowed.

---

### Requirement 12: Robust Optimization and Robust Evaluation

**User story:** As a researcher, I want robust optimization and evaluation under setup and range uncertainty scenarios.

#### Acceptance Criteria

1. The system SHALL support robust optimization through matRad-supported robust modes, including at least COWC and OWC where available.
2. The system SHALL support robust evaluation even when robust optimization is disabled.
3. The system SHALL configure setup uncertainty and range uncertainty through the configuration structure.
4. Default setup uncertainty SHALL be \(\pm 3\) mm in each axis.
5. Default range uncertainty SHALL be \(\pm 3.5\%\), if supported by the selected matRad scenario model.
6. The system SHALL report the total number of evaluated uncertainty scenarios.
7. Robust metrics SHALL include:
   - worst-case CTV D95;
   - worst-case CTV D98;
   - worst-case OAR Dmax;
   - scenario pass rate;
   - dose standard-deviation cube where available.
8. A robust plan SHALL be flagged as failing if worst-case CTV D95 falls more than 2 percentage points below the configured target threshold.
9. If robust optimization is enabled, robust re-optimization SHALL occur after final layer/spot reduction unless configured otherwise.
10. The system SHALL separately report nominal and robust metrics.

---

### Requirement 13: Planning Methods for Comparison

**User story:** As a researcher, I want to compare RMS-ELO-SPAT against appropriate baselines.

#### Acceptance Criteria

The system SHALL support the following methods:

1. **IMPT baseline:** 2–5 fixed proton beams.
2. **Dense Arc:** many-angle discrete proton arc without layer/spot sparsification.
3. **ELO-SPAT-like:** Dense Arc with \(L_{2,1/2}\) layer sparsity and optional energy sequencing.
4. **Post-filter spot removal:** threshold-based spot removal followed by re-optimization.
5. **RMS-ELO-SPAT Wrapper:** layer sparsification + sequencing score + estimated BDT-aware layer/spot selection + MMU re-optimization + robust evaluation or robust re-optimization.
6. **RMS-ELO-SPAT Integrated:** custom differentiable objective including dose, sparsity, sequencing, spot sparsity, and smooth BDT terms.

The system SHALL make clear which implementation mode was used for each method.

---

### Requirement 14: Evaluation Metrics

**User story:** As a researcher, I want comprehensive dose, robustness, delivery, and solver metrics for every method.

#### Acceptance Criteria

1. Dose metrics SHALL include:
   - CTV D98;
   - CTV D95;
   - CTV V95;
   - CTV Dmean;
   - homogeneity index \((D5-D95)/D_{\mathrm{prescribed}}\);
   - OAR Dmax;
   - OAR Dmean;
   - OAR-specific Dx/Vx where defined.
2. DVH curves SHALL be exportable for each structure.
3. Delivery metrics SHALL include:
   - number of control points;
   - number of active energy layers;
   - number of active spots;
   - number of low-to-high energy switch-up events;
   - total energy switch magnitude;
   - \(T_{\mathrm{energy}}\);
   - \(T_{\mathrm{spot}}\);
   - \(T_{\mathrm{irradiation}}\);
   - \(T_{\mathrm{gantry}}\);
   - \(T_{\mathrm{buffer}}\);
   - \(T_{\mathrm{machine}}\).
4. Robustness metrics SHALL include nominal and worst-case versions where applicable.
5. Solver metrics SHALL include:
   - final objective value;
   - dose objective component;
   - regularization components;
   - iteration count;
   - runtime;
   - convergence status.
6. The system SHALL export:
   - MATLAB struct;
   - `.mat` file;
   - CSV table;
   - LaTeX table.
7. The system SHALL generate tradeoff data pairing dose quality against BDT.

---

### Requirement 15: Iterative Wrapper Optimization Loop

**User story:** As a researcher, I want an outer-loop wrapper that alternates selection, reduction, and re-optimization so that a sparse deliverable plan can be obtained without modifying matRad core code.

#### Acceptance Criteria

1. Wrapper Mode SHALL use an outer-loop structure:
   - dense arc optimization;
   - layer scoring;
   - layer reduction;
   - re-optimization;
   - spot scoring;
   - spot/MMU reduction;
   - final re-optimization;
   - robust evaluation or robust re-optimization.
2. The wrapper SHALL use warm-start weights whenever dimensions are compatible.
3. When dimensions change due to reduction, the wrapper SHALL map previous weights to the reduced candidate set.
4. Each iteration SHALL log:
   - objective or composite score;
   - active layers;
   - active spots;
   - estimated BDT;
   - dose metrics;
   - removed layer/spot count.
5. The maximum number of wrapper iterations per stage SHALL be configurable.
6. Default maximum iterations SHALL be 20 for layer stage and 20 for spot stage.
7. The default relative convergence tolerance SHALL be \(10^{-4}\).
8. If convergence is not reached, the wrapper SHALL return the best solution according to a configured composite selection score.
9. The composite selection score SHALL be reported explicitly and SHALL NOT be confused with the matRad internal optimization objective.

---

### Requirement 16: Integrated Optimization Mode

**User story:** As a researcher, I want an optional custom optimization problem that adds RMS-ELO-SPAT penalties directly into the optimization objective.

#### Acceptance Criteria

1. Integrated Mode SHALL define a custom optimization problem class compatible with matRad optimizer callbacks.
2. The objective callback SHALL compute dose objective plus enabled RMS-ELO-SPAT penalties.
3. The gradient callback SHALL compute the corresponding gradient for all enabled differentiable terms.
4. Non-differentiable hard thresholding SHALL NOT be used inside IPOPT objective callbacks.
5. Integrated Mode SHALL include gradient tests for:
   - \(L_{2,1/2}\) group sparsity;
   - control-point barrier;
   - spot sparsity;
   - smooth BDT approximation;
   - energy-sequencing penalty if enabled.
6. A gradient test SHALL pass if relative finite-difference error is below \(10^{-4}\) for at least 95% of randomly sampled tested variables.
7. If gradient validation fails, Integrated Mode SHALL refuse to run and return diagnostic output.
8. Integrated Mode SHALL be optional and SHALL NOT block Wrapper Mode development.

---

### Requirement 17: Configuration and Reproducibility

**User story:** As a researcher, I want all parameters stored in a single configuration object so that experiments are reproducible.

#### Acceptance Criteria

1. All tunable settings SHALL be stored in a single `cfg` structure.
2. The `cfg` structure SHALL include:
   - patient file path;
   - structure names;
   - gantry angle settings;
   - optimizer settings;
   - all \(\lambda\) values;
   - all epsilon values;
   - thresholds;
   - MMU settings;
   - robust scenario settings;
   - machine timing model;
   - output path;
   - random seed;
   - implementation mode.
3. A default configuration function SHALL provide valid defaults for all fields.
4. A validation function SHALL check all required fields before execution.
5. Invalid configuration SHALL fail before optimization begins.
6. Each run SHALL save:
   - configuration;
   - git commit hash if available;
   - matRad version if available;
   - timestamp;
   - input patient identifier;
   - output metrics;
   - result objects.
7. Parameter sweeps SHALL run one independent optimization per parameter combination.
8. If one sweep setting fails, the system SHALL log the failure and continue remaining settings.

---

### Requirement 18: Safety and Clinical Disclaimer

**User story:** As a researcher, I want all outputs clearly labeled as research-only so that results are not mistaken for clinical treatment plans.

#### Acceptance Criteria

1. All reports SHALL include the statement: "Research use only. Not for clinical treatment planning."
2. Machine timing parameters SHALL be labeled as estimated unless validated.
3. The system SHALL not export DICOM RT plans for clinical delivery.
4. If DICOM export is added in future work, it SHALL be explicitly restricted to research export and require a separate validation specification.
5. Reports SHALL state that clinical deliverability requires independent machine-specific QA and vendor-validated delivery modeling.

---

## 8. Default Configuration

The default configuration SHALL be:

```matlab
cfg.mode = 'wrapper';

cfg.patientFile = 'PROSTATE.mat';
cfg.radiationMode = 'protons';
cfg.machine = 'Generic';
cfg.bioModel = 'constRBE';

cfg.imptAngles = [90 270];
cfg.arcSpacing = 20;
cfg.couchAngles = 0;

cfg.lambdaEL = 0;
cfg.lambdaB  = 0;
cfg.lambdaES = 0;
cfg.lambdaS  = 0;
cfg.lambdaT  = 0;

cfg.epsEL = 1e-8;
cfg.epsB  = 1e-8;
cfg.epsZ  = 1e-6;
cfg.epsS  = 1e-3;

cfg.layerThreshold = 1e-5;
cfg.spotThreshold  = 1e-6;
cfg.mmuThreshold   = 1e-5;

cfg.robust.enabled = false;
cfg.robust.mode = 'COWC';
cfg.robust.setupMM = [3 3 3];
cfg.robust.rangePercent = 3.5;

cfg.wrapper.maxLayerIterations = 20;
cfg.wrapper.maxSpotIterations = 20;
cfg.wrapper.relTol = 1e-4;

cfg.evaluation.ctvD95TolerancePercent = 2.0;
cfg.evaluation.robustD95TolerancePercent = 2.0;

cfg.outputDir = 'userdata/elospat_time/results';
cfg.saveIntermediate = true;
cfg.randomSeed = 1;
```

The default configuration is conservative and intended for debugging. Publication-quality experiments SHALL use parameter sweeps.

---

## 9. Default Machine Model

The default machine model SHALL be configurable and SHALL NOT be interpreted as clinically validated.

```matlab
machine.energy.upBase = 1.5;
machine.energy.downBase = 0.5;
machine.energy.perMeVUp = 0.02;
machine.energy.perMeVDown = 0.005;

machine.spot.switchTime = 0.002;

machine.irradiation.muRate = 1000;
machine.irradiation.weightUnit = 'unknown_or_MU_equivalent';

machine.gantry.maxSpeed = 6.0;
machine.gantry.maxAccel = 1.0;
machine.gantry.toleranceWindow = 0.5;

machine.buffer.perLayer = 0.05;
machine.buffer.perControlPoint = 0.02;

machine.validation.status = 'proxy_unvalidated';
machine.validation.source = 'default_placeholder';
```

If machine-specific data are later available, the system SHALL update:

```matlab
machine.validation.status = 'validated';
machine.validation.source = '<vendor/model/log/reference>';
```

---

## 10. Validation Plan

### 10.1 Unit Tests

The system SHALL include tests for:

1. bixel mapping integrity;
2. group index completeness;
3. \(L_{2,1/2}\) penalty value;
4. \(L_{2,1/2}\) gradient finite-difference agreement;
5. spot sparsity value and gradient;
6. barrier gradient;
7. BDT estimator decomposition;
8. layer ordering reproducibility;
9. `dij` reduction consistency;
10. robust scenario metric aggregation.

### 10.2 Integration Tests

The system SHALL run end-to-end on:

1. `PROSTATE.mat` with 2-beam IMPT;
2. `PROSTATE.mat` with \(20^\circ\) discrete arc;
3. robust phantom example if available;
4. at least one reduced candidate set after layer/spot filtering.

### 10.3 Research Validation Experiments

The initial research experiment SHALL compare:

1. IMPT baseline;
2. Dense Arc;
3. ELO-SPAT-like;
4. Post-filter spot removal;
5. RMS-ELO-SPAT Wrapper.

The advanced experiment MAY include:

6. RMS-ELO-SPAT Integrated.

---

## 11. Success Criteria

A result SHALL be considered successful only if all of the following hold:

1. Estimated total BDT decreases relative to Dense Arc.
2. Estimated total BDT decreases relative to ELO-SPAT-like in at least one clinically acceptable tradeoff setting.
3. Nominal CTV D95 degradation is within configured tolerance.
4. OAR metrics do not exceed configured tolerance.
5. Worst-case CTV D95 remains within configured robustness tolerance when robust evaluation is enabled.
6. Active spot count decreases or MMU distribution improves relative to Dense Arc.
7. No required control point is unintentionally emptied unless explicitly allowed.
8. The system reports all BDT components, not only total time.
9. The system clearly labels whether timing data are proxy or validated.
10. All output plans are labeled research-only.

---

## 12. Non-Goals

The following are not required for the initial implementation:

1. Clinical TPS integration.
2. Clinical DICOM RT export.
3. Vendor-specific machine delivery control.
4. Real-time gantry control.
5. Biological variable-RBE optimization.
6. Motion interplay modeling beyond available matRad scenario support.
7. Full mixed-integer optimization with exact binary \(z_{b,e}\).
8. Proof of global optimality.
9. Clinical deliverability certification.

---

## 13. Known Risks and Mitigations

| Risk | Impact | Mitigation |
|---|---|---|
| `dij` too large for fine arc spacing | Memory failure | Start with \(20^\circ\) and \(10^\circ\); reduce dose grid |
| BDT model not clinically validated | Overclaiming | Label as estimated proxy unless validated |
| Gradient errors in Integrated Mode | IPOPT failure | Require finite-difference gradient tests |
| Hard thresholding degrades dose | Loss of target coverage | Re-optimize and revert if CTV D95 tolerance exceeded |
| Spot weights not in MU | Incorrect irradiation time | Require MU conversion or mark irradiation time approximate |
| Empty control points after sparsification | Poor arc continuity | Use barrier or minimum active-layer rules |
| Robust optimization too slow | Long runtime | Run nominal first; robust final; use coarse angle spacing |
| Sequencing penalty without gradient | Invalid integrated optimization | Use only in wrapper score unless gradient is implemented |
| Serial BDT overestimates dynamic delivery | Conservative timing bias | Report assumption and allow future synchronized model |

---

## 14. Reporting Requirements

Every completed run SHALL produce a report containing:

1. configuration summary;
2. implementation mode;
3. dataset and geometry summary;
4. method comparison table;
5. DVH metrics;
6. robust metrics if available;
7. BDT decomposition;
8. active layer/spot statistics;
9. solver status;
10. warnings and failures;
11. research-only disclaimer.

The LaTeX table export SHALL include at minimum:

| Method | CTV D95 | CTV D98 | OAR Dmax | Active EL | Active Spots | Switch-up | \(T_{\mathrm{energy}}\) | \(T_{\mathrm{spot}}\) | \(T_{\mathrm{gantry}}\) | \(T_{\mathrm{total}}\) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|

---

## 15. Final Statement

This specification defines RMS-ELO-SPAT as a research framework that operationalizes the literature gap left by ELO-SPAT. The framework preserves the mathematically correct \(L_{2,1/2}\) energy-layer group sparsity of ELO-SPAT, but extends the evaluation and optimization target toward estimated machine-specific beam delivery time, spot sparsity, MMU-aware re-optimization, robust scenario performance, and deliverability constraints.

The first implementation SHALL be a wrapper-based prototype in `userdata/elospat_time/`. A later integrated version MAY add differentiable RMS-ELO-SPAT penalties directly into a custom matRad optimization problem class after gradient validation.
