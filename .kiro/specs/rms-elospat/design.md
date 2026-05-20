# Technical Design: RMS-ELO-SPAT

## Overview

This document describes the technical design for the RMS-ELO-SPAT (Robust Machine-specific Sparse ELO-SPAT) research framework. The system is implemented as an external MATLAB module layer on top of matRad, located in `userdata/elospat_time/`.

The design follows a modular architecture with clear separation between planning, mapping, objective computation, delivery modeling, reduction, robustness, and evaluation concerns.

---

## Architecture

### 1.1 High-Level Component Diagram

```
┌─────────────────────────────────────────────────────────────────────┐
│                        run_rms_elospat_main.m                       │
│                     (Orchestrator / Entry Point)                     │
└───────────┬─────────────────────────────────────────────────────────┘
            │
            ▼
┌───────────────────┐   ┌──────────────────┐   ┌──────────────────────┐
│   config/         │   │   planning/      │   │   mapping/           │
│                   │   │                  │   │                      │
│ default_config    │──▶│ build_impt_plan  │──▶│ map_bixel_to_cp_     │
│ default_machine   │   │ build_arc_plan   │   │   energy_spot        │
│ validate_config   │   │ run_dense_arc    │   │ build_group_index    │
└───────────────────┘   │ run_elospat_like │   │ validate_mapping     │
                        │ run_spot_sparse  │   │ get_group_indices    │
                        │ run_final_reopt  │   └──────────────────────┘
                        └──────────────────┘
            │
            ▼
┌───────────────────┐   ┌──────────────────┐   ┌──────────────────────┐
│   objective/      │   │   delivery/      │   │   reduction/         │
│                   │   │                  │   │                      │
│ l21half_penalty   │   │ extract_sequence │   │ build_layer_mask     │
│ l21half_gradient  │   │ order_layers     │   │ build_spot_mask      │
│ cp_barrier        │   │ estimate_bdt     │   │ cut_dij_by_mask      │
│ cp_barrier_grad   │   │ estimate_energy  │   │ cut_stf_by_mask      │
│ sequencing_pen    │   │ estimate_spot    │   │ reoptimize_reduced   │
│ spot_sparsity     │   │ estimate_gantry  │   └──────────────────────┘
│ spot_sparsity_gr  │   │ validate_gantry  │
│ composite_score   │   └──────────────────┘
└───────────────────┘
            │
            ▼
┌───────────────────┐   ┌──────────────────────────────────────────────┐
│   robustness/     │   │   evaluation/                                │
│                   │   │                                              │
│ configure_robust  │   │ compute_dvh_metrics    compare_methods_table │
│ configure_obj     │   │ compute_delivery_met   export_latex_table    │
│ evaluate_robust   │   │ compute_solver_met     export_tradeoff_data  │
│ worst_case_met    │   └──────────────────────────────────────────────┘
└───────────────────┘
```


### 1.2 Data Flow Diagram

```
┌──────────┐    ┌──────────┐    ┌──────────┐    ┌──────────┐
│ Patient  │───▶│ matRad   │───▶│ matRad   │───▶│ matRad   │
│ Data     │    │ Plan     │    │ STF Gen  │    │ Dij Calc │
│(ct,cst)  │    │ Setup    │    │          │    │          │
└──────────┘    └──────────┘    └──────────┘    └────┬─────┘
                                                     │
                                                     ▼
┌──────────┐    ┌──────────┐    ┌──────────┐    ┌──────────┐
│ Delivery │◀───│ Layer    │◀───│ Bixel    │◀───│ Baseline │
│ Sequence │    │ Scoring  │    │ Mapping  │    │ FluOpt   │
│ Extract  │    │ & Filter │    │ & Groups │    │          │
└────┬─────┘    └────┬─────┘    └──────────┘    └──────────┘
     │               │
     ▼               ▼
┌──────────┐    ┌──────────┐    ┌──────────┐    ┌──────────┐
│ BDT      │    │ Dij/STF  │───▶│ Re-Opt   │───▶│ Spot     │
│ Estimate │    │ Reduction│    │ (Layer)  │    │ Scoring  │
│          │    │          │    │          │    │ & Filter │
└──────────┘    └──────────┘    └──────────┘    └────┬─────┘
                                                     │
                                                     ▼
┌──────────┐    ┌──────────┐    ┌──────────┐    ┌──────────┐
│ Final    │◀───│ Robust   │◀───│ Final    │◀───│ Dij/STF  │
│ Metrics  │    │ Eval     │    │ Re-Opt   │    │ Reduction│
│ & Export │    │          │    │(Robust)  │    │ (Spots)  │
└──────────┘    └──────────┘    └──────────┘    └──────────┘
```

### 1.3 Dependency on matRad

The system depends on the following matRad functions (called but NOT modified):

| matRad Function | Purpose | Note |
|---|---|---|
| `matRad_rc` | Initialize matRad runtime context | Required |
| `matRad_generateStf` | Generate beam steering information | Required |
| `matRad_calcDoseInfluence` | Compute dose influence matrix | Required |
| `matRad_fluenceOptimization` | Run IPOPT-based fluence optimization | Required |
| `matRad_calcDoseForward` / `matRad_calcCubes` / equivalent | Compute dose from weights for evaluation | Version-dependent; exact function name to be confirmed after inspecting installed matRad version |
| `matRad_WorstCaseScenarios` | Create robust scenario set | Required for robust mode |

**Note on dose recalculation:** matRad's dose forward-calculation API varies across versions. The implementation SHALL detect the available function at runtime (checking for `matRad_calcDoseForward`, `matRad_calcCubes`, or `matRad_calcDoseDirect` in order of preference) and use the first available. If none is found, the system SHALL fall back to manual matrix-vector multiplication `dij.physicalDose{1} * w` for nominal dose computation.

---

## Components and Interfaces

### 2.1 Configuration Module (`config/`)

#### 2.1.1 `default_config.m`

```matlab
function cfg = default_config()
% Returns the default configuration structure with all fields initialized.
%
% Output:
%   cfg - struct with fields for all tunable parameters
```

**Key design decisions:**
- Single struct contains ALL parameters (no nested config files)
- Every field has a documented default and valid range
- Mode field selects 'wrapper' or 'integrated'
- Random seed ensures reproducibility


**Data model:**

```matlab
cfg.mode              % 'wrapper' | 'integrated'
cfg.patientFile       % string path to .mat patient file
cfg.radiationMode     % 'protons'
cfg.machine           % 'Generic' (matRad machine name)
cfg.bioModel          % 'constRBE'

% Geometry
cfg.imptAngles        % [90 270] default IMPT angles
cfg.arcSpacing        % 20 (degrees)
cfg.couchAngles       % 0

% Penalty weights
cfg.lambdaEL          % 0 (L2,1/2 group sparsity)
cfg.lambdaB           % 0 (log barrier)
cfg.lambdaES          % 0 (energy sequencing)
cfg.lambdaS           % 0 (spot sparsity)
cfg.lambdaT           % 0 (BDT proxy)

% Smoothing parameters
cfg.epsEL             % 1e-8
cfg.epsB              % 1e-8
cfg.epsZ              % 1e-6
cfg.epsS              % 1e-3

% Thresholds
cfg.layerThreshold    % 1e-5
cfg.spotThreshold     % 1e-6
cfg.mmuThreshold      % 1e-5

% Robustness
cfg.robust.enabled    % false
cfg.robust.mode       % 'COWC'
cfg.robust.setupMM    % [3 3 3]
cfg.robust.rangePercent % 3.5

% Wrapper loop
cfg.wrapper.maxLayerIterations  % 20
cfg.wrapper.maxSpotIterations   % 20
cfg.wrapper.relTol              % 1e-4

% Evaluation tolerances
cfg.evaluation.ctvD95TolerancePercent    % 2.0
cfg.evaluation.robustD95TolerancePercent % 2.0

% Output
cfg.outputDir         % 'userdata/elospat_time/results'
cfg.saveIntermediate  % true
cfg.randomSeed        % 1
```

#### 2.1.2 `default_machine_model.m`

```matlab
function machine = default_machine_model()
% Returns the default machine timing model (proxy/unvalidated).
%
% Output:
%   machine - struct with energy, spot, irradiation, gantry, buffer fields
```

**Data model:**

```matlab
machine.energy.upBase       % 1.5 s
machine.energy.downBase     % 0.5 s
machine.energy.perMeVUp     % 0.02 s/MeV
machine.energy.perMeVDown   % 0.005 s/MeV

machine.spot.switchTime     % 0.002 s

machine.irradiation.muRate  % 1000 MU/s
machine.irradiation.weightUnit % 'unknown_or_MU_equivalent'

machine.gantry.maxSpeed     % 6.0 deg/s
machine.gantry.maxAccel     % 1.0 deg/s^2
machine.gantry.toleranceWindow % 0.5 deg

machine.buffer.perLayer         % 0.05 s
machine.buffer.perControlPoint  % 0.02 s

machine.validation.status   % 'proxy_unvalidated'
machine.validation.source   % 'default_placeholder'
```


#### 2.1.3 `validate_config.m`

```matlab
function validate_config(cfg, machine)
% Validates all configuration fields before optimization begins.
% Throws error with field name if invalid.
%
% Checks:
%   - All required fields exist
%   - Numeric fields are within valid ranges
%   - String fields match allowed values
%   - Patient file exists on disk
%   - Arc spacing divides 360 (or custom angles provided)
```

---

### 2.2 Planning Module (`planning/`)

#### 2.2.1 `build_impt_plan.m`

```matlab
function [pln, stf, dij, result] = build_impt_plan(ct, cst, cfg)
% Builds and optimizes a standard IMPT plan as baseline.
%
% Inputs:
%   ct  - matRad CT structure
%   cst - matRad structure/objective table
%   cfg - configuration struct
%
% Outputs:
%   pln    - matRad plan structure
%   stf    - steering information
%   dij    - dose influence matrix
%   result - optimization result (resultGUI)
%
% Algorithm:
%   1. Set pln fields (radiationMode, machine, bioModel, gantryAngles)
%   2. Call matRad_generateStf
%   3. Call matRad_calcDoseInfluence
%   4. Call matRad_fluenceOptimization
%   5. Check convergence status
%   6. Return all structures
```

#### 2.2.2 `build_discrete_arc_plan.m`

```matlab
function [pln, stf, dij] = build_discrete_arc_plan(ct, cst, cfg)
% Builds a discrete proton arc geometry and computes dij.
%
% Inputs:
%   ct  - matRad CT structure
%   cst - matRad structure/objective table
%   cfg - configuration struct (uses cfg.arcSpacing)
%
% Outputs:
%   pln - matRad plan with many gantry angles
%   stf - steering information for all control points
%   dij - dose influence matrix for full arc
%
% Algorithm:
%   1. Generate gantry angles: 0:cfg.arcSpacing:(360-cfg.arcSpacing)
%   2. Set couch angles to zeros
%   3. Configure pln for protons
%   4. Call matRad_generateStf
%   5. Call matRad_calcDoseInfluence
%   6. Validate dij.totalNumOfBixels > 0
```

#### 2.2.3 `run_dense_arc_baseline.m`

```matlab
function result = run_dense_arc_baseline(dij, cst, pln)
% Runs matRad fluence optimization on the full dense arc.
%
% This is the "Dense Arc" baseline method.
```


#### 2.2.4 `run_elospat_like_stage.m`

```matlab
function [result, layerMask] = run_elospat_like_stage(dij, cst, pln, ...
    wInit, groups, cpGroups, cfg, machine)
% ELO-SPAT-inspired wrapper layer sparsification stage.
%
% IMPORTANT: This is NOT an exact reproduction of the original ELO-SPAT
% solver. The original ELO-SPAT uses a unified differentiable objective
% with L2,1/2 + log barrier + sequencing penalty solved jointly by IPOPT.
% This wrapper approximation uses iterative layer scoring, filtering, and
% re-optimization as an outer loop around matRad's standard optimizer.
%
% For exact integrated ELO-SPAT-like optimization within IPOPT, see
% Future Integrated Mode (Section 6) with custom
% matRad_OptimizationProblemELOSPAT class.
%
% Algorithm (Wrapper Mode — ELO-SPAT-like layer filtering):
%   1. Initialize w = wInit
%   2. For iter = 1:cfg.wrapper.maxLayerIterations:
%      a. Compute composite score:
%         score = lambdaEL*R_EL + lambdaB*B + lambdaES*Phi_ES
%      b. Compute layer importance scores (see Section 4.1)
%      c. Compute estimated BDT saving per candidate layer
%      d. Rank layers by removal score (low importance + high BDT saving)
%      e. Identify candidate layers for removal
%      f. Build layer mask (keep active layers)
%      g. Expand layer mask to bixel mask
%      h. Cut dij/stf by bixel mask
%      i. Re-optimize on reduced problem (warm-start)
%      j. Check convergence (relative change in score < relTol)
%      k. Log iteration metrics
%   3. Return final result and layer mask
```

#### 2.2.5 `run_spot_sparse_stage.m`

```matlab
function [result, spotMask] = run_spot_sparse_stage(dij, cst, pln, ...
    wInit, map, groups, cfg, machine)
% Iterative spot sparsity and MMU-aware filtering stage.
%
% Algorithm (Wrapper Mode):
%   1. Initialize w = wInit
%   2. For iter = 1:cfg.wrapper.maxSpotIterations:
%      a. Compute spot sparsity score
%      b. Compute BDT estimate
%      c. Identify sub-MMU spots
%      d. Build spot mask
%      e. Cut dij/stf by spot mask
%      f. Re-optimize on reduced problem
%      g. Check CTV D95 degradation
%      h. If degradation > tolerance: revert and stop
%      i. Check convergence
%      j. Log iteration metrics
%   3. Return final result and spot mask
```

#### 2.2.6 `run_final_reoptimization.m`

```matlab
function result = run_final_reoptimization(dij, cst, pln, wInit, cfg)
% Final re-optimization, optionally with robust objectives.
%
% If cfg.robust.enabled:
%   1. Configure robust scenarios on pln
%   2. Set robustness mode on cst objectives
%   3. Call matRad_fluenceOptimization with robust settings
% Else:
%   1. Call matRad_fluenceOptimization with nominal settings
```

---

### 2.3 Mapping Module (`mapping/`)

#### 2.3.1 `map_bixel_to_cp_energy_spot.m`

```matlab
function map = map_bixel_to_cp_energy_spot(dij, stf, pln)
% Maps every bixel to its (control_point, energy_layer, spot) identity.
%
% Inputs:
%   dij - matRad dose influence matrix (uses .beamNum, .rayNum, .bixelNum)
%   stf - matRad steering info (uses stf(b).ray(r).energy)
%   pln - matRad plan (uses pln.propStf.gantryAngles)
%
% Output:
%   map - struct array of length dij.totalNumOfBixels with fields:
%         .globalIx    - 1:N sequential index
%         .beamNum     - beam/control-point number
%         .rayNum      - ray number within beam
%         .bixelNum    - bixel number within ray
%         .gantryAngle - gantry angle in degrees
%         .energy      - physical energy (MeV)
%         .cpIndex     - control point index (= beamNum)
%         .energyIndex - energy layer index within control point
%
% Algorithm:
%   For each bixel i = 1:N:
%     b = dij.beamNum(i)
%     r = dij.rayNum(i)
%     s = dij.bixelNum(i)
%     energy = stf(b).ray(r).energy(s)
%     map(i) = {i, b, r, s, gantryAngles(b), energy, b, TBD}
%
%   For each beam b:
%     Find unique energies among bixels in beam b
%     Assign energyIndex = position in sorted unique energy list
```


#### 2.3.2 `build_group_index.m`

```matlab
function [groups, cpGroups] = build_group_index(map)
% Builds group index structures from the bixel mapping.
%
% Outputs:
%   groups   - struct array, one per (cpIndex, energyIndex) pair
%              .cpIndex, .energyIndex, .ix (bixel indices), .nSpots
%   cpGroups - struct array, one per control point
%              .cpIndex, .ix (all bixel indices at this CP)
%
% Algorithm:
%   1. Extract unique (cpIndex, energyIndex) pairs
%   2. For each pair, collect all bixel indices
%   3. Extract unique cpIndex values
%   4. For each CP, collect all bixel indices
```

#### 2.3.3 `validate_mapping.m`

```matlab
function validate_mapping(map, dij)
% Validates mapping integrity.
%
% Checks:
%   1. numel(map) == dij.totalNumOfBixels
%   2. [map.globalIx] == 1:dij.totalNumOfBixels (no gaps/duplicates)
%   3. Every bixel appears in exactly one group
%
% Throws error with affected indices if validation fails.
```

#### 2.3.4 `get_group_indices.m`

```matlab
function ix = get_group_indices(groups, cpIndex, energyIndex)
% Returns bixel indices for a given (cpIndex, energyIndex) group.
% Returns empty [] if group does not exist.
```

#### 2.3.5 `get_control_point_indices.m`

```matlab
function ix = get_control_point_indices(cpGroups, cpIndex)
% Returns all bixel indices at a given control point.
% Returns empty [] if control point does not exist.
```

---

### 2.4 Objective Module (`objective/`)

#### 2.4.1 `l21half_group_penalty.m`

```matlab
function [val, grad] = l21half_group_penalty(w, groups, epsEL)
% Computes L2,1/2 group sparsity penalty and gradient.
%
% Formula:
%   val = sum_g (sum_{i in g} w_i^2 + epsEL)^(1/4)
%   grad_i = 0.5 * w_i * (sum_{j in g} w_j^2 + epsEL)^(-3/4)
%
% Inputs:
%   w      - spot weight vector (N x 1)
%   groups - struct array with .ix field
%   epsEL  - smoothing parameter (default 1e-8)
%
% Outputs:
%   val  - scalar penalty value
%   grad - gradient vector (N x 1)
```

**Implementation detail:** Loop over groups, vectorize within each group.

#### 2.4.2 `control_point_barrier.m`

```matlab
function [val, grad] = control_point_barrier(w, cpGroups, epsB)
% Computes log-barrier penalty for control-point activity preservation.
%
% Formula:
%   val = -sum_b log(sqrt(sum_{i in C_b} w_i^2 + epsB))
%   grad_i = -w_i / (sum_{j in C_b} w_j^2 + epsB)  for i in C_b
%
% Inputs:
%   w        - spot weight vector (N x 1)
%   cpGroups - struct array with .ix field (one per control point)
%   epsB     - barrier parameter (default 1e-8)
%
% Outputs:
%   val  - scalar barrier value
%   grad - gradient vector (N x 1)
```


#### 2.4.3 `energy_sequencing_penalty.m`

```matlab
function [val, switchUpCount] = energy_sequencing_penalty(w, groups, epsZ)
% Computes energy sequencing penalty (switch-up cost).
%
% Formula:
%   a_{b,e} = ||x_{b,e}||_2
%   z_tilde_{b,e} = a_{b,e} / (a_{b,e} + epsZ)
%   Phi_ES = sum_{b=1}^{B-1} sum_{e,e'} z_tilde_{b,e} * z_tilde_{b+1,e'}
%            * max(E_{b+1,e'} - E_{b,e}, 0)
%
% Inputs:
%   w      - spot weight vector
%   groups - struct array with .ix, .cpIndex, .energyIndex, .energy fields
%   epsZ   - soft activity smoothing (default 1e-6)
%
% Outputs:
%   val          - smooth sequencing penalty value
%   switchUpCount - hard count of switch-up events (for reporting)
%
% Note: Gradient not computed here. Use finite difference for validation.
```

#### 2.4.4 `spot_sparsity_penalty.m`

```matlab
function [val, grad] = spot_sparsity_penalty(w, epsS)
% Computes smooth L0-like spot sparsity penalty and gradient.
%
% Formula:
%   val = sum_i w_i / (w_i + epsS)
%   grad_i = epsS / (w_i + epsS)^2
%
% Inputs:
%   w    - spot weight vector (N x 1), assumed >= 0
%   epsS - smoothing parameter (default 1e-3)
%
% Outputs:
%   val  - scalar penalty value
%   grad - gradient vector (N x 1)
```

#### 2.4.5 `composite_score.m`

```matlab
function [score, components] = composite_score(w, groups, cpGroups, ...
    cfg, machine, map, stf)
% Computes the full composite objective score for wrapper mode.
%
% score = lambdaEL * R_EL + lambdaB * B + lambdaES * Phi_ES
%         + lambdaS * R_spot + lambdaT * T_machine
%
% Outputs:
%   score      - total composite score (scalar)
%   components - struct with individual penalty values:
%                .l21half, .barrier, .sequencing, .spotSparsity, .bdt
```

---

### 2.5 Delivery Module (`delivery/`)

#### 2.5.1 `extract_delivery_sequence.m`

```matlab
function seq = extract_delivery_sequence(w, map, groups, cfg)
% Extracts the delivery sequence from optimized weights.
%
% Algorithm:
%   1. Identify active layers (||x_{b,e}||_2 > cfg.layerThreshold)
%   2. Identify active spots (x_i > cfg.spotThreshold)
%   3. Order control points by gantry angle (ascending)
%   4. Within each CP, order layers by energy (descending, default)
%   5. Build sequence struct array
%
% Output:
%   seq - struct array with fields per delivered layer:
%         .cpIndex, .gantryAngle, .energy, .energyIndex
%         .spotIndices, .spotWeights, .numActiveSpots
%         .layerWeightSum, .prevEnergy, .nextEnergy
```

#### 2.5.2 `order_layers_within_control_point.m`

```matlab
function orderedEnergies = order_layers_within_control_point(energies, mode)
% Orders energy layers within a control point.
%
% Inputs:
%   energies - vector of active energy values at this CP
%   mode     - 'descending' (default) | 'ascending' | 'custom'
%
% Output:
%   orderedEnergies - sorted energy values
```


#### 2.5.3 `estimate_bdt.m`

```matlab
function metrics = estimate_bdt(seq, machine)
% Estimates total beam delivery time from delivery sequence using the
% configured machine timing model.
%
% TERMINOLOGY: This function computes an "estimated machine-specific BDT
% proxy" based on configurable timing parameters. It should only be called
% "validated BDT" if machine.validation.status == 'validated' (i.e., timing
% parameters come from vendor specifications or measured delivery logs).
%
% Inputs:
%   seq     - delivery sequence struct array
%   machine - machine timing model struct
%
% Output:
%   metrics - struct with fields:
%             .T_energy, .T_spot, .T_irradiation, .T_gantry, .T_buffer
%             .T_total (sum of all components)
%             .nActiveLayers, .nActiveSpots, .nSwitchUp
%             .totalSwitchMagnitude
%
% Algorithm:
%   T_energy: sum direction-dependent switching times between consecutive
%             layers (up: upBase + perMeVUp*|dE|, down: downBase + perMeVDown*|dE|)
%   T_spot:   sum max(nSpots-1, 0) * spotSwitchTime per layer
%   T_irrad:  sum(all spot weights) / muRate
%   T_gantry: sum trapezoidal motion time between consecutive CPs
%   T_buffer: nLayers * perLayer + nCP * perControlPoint
```

#### 2.5.4 `estimate_gantry_time.m`

```matlab
function t = estimate_gantry_time(theta, machine)
% Computes gantry rotation time using trapezoidal motion profile.
%
% Formula:
%   If theta <= vmax^2/amax:  t = 2*sqrt(theta/amax)
%   Else:                     t = theta/vmax + vmax/amax
%
% Inputs:
%   theta   - angular difference in degrees
%   machine - machine model (uses .gantry.maxSpeed, .gantry.maxAccel)
%
% Output:
%   t - estimated rotation time in seconds
```

#### 2.5.5 `validate_gantry_constraints.m`

```matlab
function [feasible, violations] = validate_gantry_constraints(seq, machine)
% Checks if the delivery sequence respects gantry motion constraints.
%
% Checks:
%   - Gantry speed does not exceed maxSpeed at any transition
%   - Angular differences are within tolerance window where needed
%
% Outputs:
%   feasible   - boolean
%   violations - struct array of constraint violations (if any)
```

---

### 2.6 Reduction Module (`reduction/`)

#### 2.6.1 `build_layer_mask.m`

```matlab
function layerMask = build_layer_mask(w, groups, threshold)
% Identifies active energy layers based on group L2 norm.
%
% A layer group g is active if: ||w(groups(g).ix)||_2 > threshold
%
% Output:
%   layerMask - logical vector (nGroups x 1), true = keep
```

#### 2.6.2 `build_spot_mask.m`

```matlab
function spotMask = build_spot_mask(w, threshold)
% Identifies active spots based on individual weight threshold.
%
% A spot i is active if: w(i) > threshold
%
% Output:
%   spotMask - logical vector (N x 1), true = keep
```

#### 2.6.2b `expand_layer_mask_to_bixel_mask.m`

```matlab
function bixelMask = expand_layer_mask_to_bixel_mask(layerMask, groups, totalNumOfBixels)
% Converts a group-level layer mask to a bixel-level mask.
%
% Inputs:
%   layerMask       - logical vector (nGroups x 1), true = keep group
%   groups          - struct array with .ix field (bixel indices per group)
%   totalNumOfBixels - total number of bixels (for output size)
%
% Output:
%   bixelMask - logical vector (totalNumOfBixels x 1), true = keep bixel
%
% Algorithm:
%   bixelMask = false(totalNumOfBixels, 1);
%   for g = 1:numel(groups)
%       if layerMask(g)
%           bixelMask(groups(g).ix) = true;
%       end
%   end
%
% IMPORTANT: layerMask has length nGroups (one per energy-layer group),
% while bixelMask has length totalNumOfBixels. These are different
% dimensions and must NOT be confused. cut_dij_by_mask requires a
% bixel-level mask.
```

#### 2.6.3 `cut_dij_by_mask.m`

```matlab
function [dijReduced, mapReduced, wReduced] = cut_dij_by_mask(dij, mask, map, w)
% Removes inactive bixels from dij based on a bixel-level boolean mask.
%
% Inputs:
%   dij  - original dose influence matrix
%   mask - logical vector (N x 1) at BIXEL level, true = keep
%          NOTE: this must be a bixel-level mask, NOT a group-level mask.
%          Use expand_layer_mask_to_bixel_mask() to convert if needed.
%   map  - original bixel mapping
%   w    - original weight vector
%
% Outputs:
%   dijReduced - reduced dij with only active columns
%   mapReduced - reduced mapping (preserves original globalIx for tracing)
%   wReduced   - reduced weight vector (warm-start init)
%
% Algorithm:
%   1. For each dose-influence-like field that stores per-bixel data as
%      columns, keep only columns where mask is true. Fields handled:
%        - dij.physicalDose{s}  (cell array, one per scenario)
%        - dij.mAlphaDose{s}    (if exists; cell array per scenario)
%        - dij.mSqrtBetaDose{s} (if exists; cell array per scenario)
%        - dij.mLETDose{s}      (if exists; cell array per scenario)
%      For each field: iterate over all scenarios s = 1:numel(field).
%      If the field does not exist, skip safely.
%
%   2. For per-bixel vector fields, keep only entries where mask is true:
%        - dij.beamNum(mask)
%        - dij.rayNum(mask)
%        - dij.bixelNum(mask)
%        - dij.minMU(mask)          (if exists)
%        - dij.maxMU(mask)          (if exists)
%        - dij.numParticlesPerMU(mask) (if exists)
%
%   3. Update dij.totalNumOfBixels = sum(mask)
%
%   4. Rebuild mapReduced from map(mask), preserving original globalIx
%      for reverse traceability.
%
%   5. wReduced = w(mask)
%
% Robustness:
%   - If a field is not present in dij, it is skipped without error.
%   - Cell array fields are iterated over all cells (scenarios).
%   - Sparse matrices are sliced using column indexing: M(:, mask).
%   - The function asserts numel(mask) == dij.totalNumOfBixels before
%     proceeding; throws error if mismatch.
```


#### 2.6.4 `cut_stf_by_mask.m`

```matlab
function stfReduced = cut_stf_by_mask(stf, mask, map)
% Removes inactive rays/bixels from stf based on mask.
%
% Must maintain consistency with cut_dij_by_mask output.
% Removes empty rays and adjusts ray/bixel counts per beam.
```

#### 2.6.5 `reoptimize_reduced_problem.m`

```matlab
function result = reoptimize_reduced_problem(dijReduced, cst, pln, wInit)
% Re-optimizes on a reduced dij with warm-start.
%
% Calls matRad_fluenceOptimization with the reduced problem.
% Uses wInit as starting point if dimensions match.
```

---

### 2.7 Robustness Module (`robustness/`)

#### 2.7.1 `configure_robust_scenarios.m`

```matlab
function pln = configure_robust_scenarios(pln, ct, cfg)
% Configures matRad worst-case scenarios on the plan.
%
% Sets:
%   pln.multScen = matRad_WorstCaseScenarios(ct)
%   pln.multScen.shiftSD = cfg.robust.setupMM / cfg.robust.wcFactor
%   pln.multScen.wcFactor = 3 (or configured)
```

#### 2.7.2 `configure_robust_objectives.m`

```matlab
function cst = configure_robust_objectives(cst, cfg)
% Sets robustness mode on all objectives in cst.
%
% For each structure with objectives:
%   cst{i,6}{j}.robustness = cfg.robust.mode  (e.g., 'COWC')
```

#### 2.7.3 `evaluate_robust_scenarios.m`

```matlab
function robustEval = evaluate_robust_scenarios(ct, cst, stf, dij, pln, w, cfg)
% Evaluates dose across all uncertainty scenarios.
%
% Inputs:
%   ct  - matRad CT structure
%   cst - matRad structure/objective table
%   stf - matRad steering info (may be needed for scenario dose calc;
%          optional if dij already contains scenario-specific matrices)
%   dij - dose influence matrix (may contain multi-scenario cells)
%   pln - matRad plan structure (with multScen configured)
%   w   - optimized spot weight vector
%   cfg - configuration struct
%
% Note: If dij.physicalDose is a cell array with numel > 1, each cell
% corresponds to a pre-computed uncertainty scenario and stf is not
% required for dose recalculation. If only nominal dij is available,
% stf is needed to recompute scenario-shifted dose influence matrices.
%
% Algorithm:
%   1. Generate scenario doses using matRad scenario framework
%   2. Compute DVH metrics per scenario
%   3. Identify worst-case scenario
%   4. Compute scenario pass rate
%   5. Compute dose standard deviation cube
%
% Output:
%   robustEval - struct with:
%     .worstCase.ctvD95, .worstCase.ctvD98, .worstCase.oarDmax
%     .scenarioPassRate
%     .stdDoseCube
%     .nScenarios
%     .allScenarioMetrics (table)
%     .isRobustOptimized (bool: was robust optimization used, or eval only?)
```

#### 2.7.4 `compute_worst_case_metrics.m`

```matlab
function wcMetrics = compute_worst_case_metrics(scenarioMetrics)
% Extracts worst-case values from per-scenario metric table.
%
% Worst-case CTV D95 = min across scenarios
% Worst-case OAR Dmax = max across scenarios
```

---

### 2.8 Evaluation Module (`evaluation/`)

#### 2.8.1 `compute_dvh_metrics.m`

```matlab
function doseMetrics = compute_dvh_metrics(doseCube, ct, cst, cfg)
% Computes dose-volume histogram metrics for all structures.
%
% Output fields per structure:
%   .D98, .D95, .D5, .V95, .Dmean, .Dmax
%   .homogeneityIndex = (D5 - D95) / Dprescribed
%   .dvhCurve (dose bins at 0.1 Gy resolution)
```


#### 2.8.2 `compute_delivery_metrics.m`

```matlab
function delivMetrics = compute_delivery_metrics(w, map, groups, machine, cfg)
% Computes all delivery-related metrics.
%
% Output:
%   .nControlPoints, .nActiveLayers, .nActiveSpots
%   .nSwitchUpEvents, .totalSwitchMagnitude
%   .T_energy, .T_spot, .T_irradiation, .T_gantry, .T_buffer, .T_total
%   .minSpotWeight, .medianSpotWeight, .pctSubMMU
```

#### 2.8.3 `compute_solver_metrics.m`

```matlab
function solverMetrics = compute_solver_metrics(result, startTime)
% Extracts solver performance metrics from matRad result.
%
% Output:
%   .finalObjective, .doseObjective, .regularizationValue
%   .nIterations, .runtime, .convergenceStatus
```

#### 2.8.4 `compare_methods_table.m`

```matlab
function T = compare_methods_table(allResults)
% Builds a comparison table across all methods.
%
% Input:
%   allResults - struct with fields: .impt, .denseArc, .elospatLike,
%                .postFilter, .rmsElospat
%
% Output:
%   T - MATLAB table with one row per method, columns for all metrics
```

#### 2.8.5 `export_latex_table.m`

```matlab
function export_latex_table(T, filepath)
% Exports comparison table as LaTeX tabular environment.
%
% Columns: Method | CTV D95 | CTV D98 | OAR Dmax | Active EL |
%          Active Spots | Switch-up | T_energy | T_spot | T_gantry | T_total
```

#### 2.8.6 `export_tradeoff_data.m`

```matlab
function export_tradeoff_data(allResults, filepath)
% Exports dose quality vs BDT tradeoff data for plotting.
%
% Output: .mat file with paired (doseMetric, bdtMetric) per method
%         and per hyperparameter setting for Pareto front visualization.
```

---

## 3. Orchestrator Design (`run_rms_elospat_main.m`)

### 3.1 Main Script Pseudocode

```matlab
function results = run_rms_elospat_main(cfgOverride)
% Main entry point for RMS-ELO-SPAT research pipeline.

%% 0. Setup
matRad_rc;
cfg = default_config();
if nargin > 0, cfg = merge_config(cfg, cfgOverride); end
machine = default_machine_model();
validate_config(cfg, machine);
rng(cfg.randomSeed);

%% 1. Load patient data
load(cfg.patientFile, 'ct', 'cst');

%% 2. IMPT Baseline (Method 1)
[plnIMPT, stfIMPT, dijIMPT, resultIMPT] = build_impt_plan(ct, cst, cfg);

%% 3. Discrete Arc Setup
[plnArc, stfArc, dijArc] = build_discrete_arc_plan(ct, cst, cfg);

%% 4. Dense Arc Baseline (Method 2)
resultDenseArc = run_dense_arc_baseline(dijArc, cst, plnArc);

%% 5. Build Mapping and Groups
map = map_bixel_to_cp_energy_spot(dijArc, stfArc, plnArc);
validate_mapping(map, dijArc);
[groups, cpGroups] = build_group_index(map);

%% 6. ELO-SPAT-like Stage (Method 3)
[resultELO, layerMask] = run_elospat_like_stage(...
    dijArc, cst, plnArc, resultDenseArc.w, groups, cpGroups, cfg, machine);

%% 7. Layer Reduction
bixelMaskLayer = expand_layer_mask_to_bixel_mask(layerMask, groups, dijArc.totalNumOfBixels);
[dijLayer, mapLayer, wLayerInit] = cut_dij_by_mask(dijArc, bixelMaskLayer, map, resultELO.w);
stfLayer = cut_stf_by_mask(stfArc, bixelMaskLayer, map);
resultLayer = reoptimize_reduced_problem(dijLayer, cst, plnArc, wLayerInit);

%% 8. Spot Sparsity + BDT Stage (Method 5 - RMS-ELO-SPAT)
[groupsLayer, cpGroupsLayer] = build_group_index(mapLayer);
[resultSparse, spotMask] = run_spot_sparse_stage(...
    dijLayer, cst, plnArc, resultLayer.w, mapLayer, groupsLayer, cfg, machine);

%% 9. Spot Reduction + Final Re-optimization
[dijFinal, mapFinal, wFinalInit] = cut_dij_by_mask(dijLayer, spotMask, mapLayer, resultSparse.w);
stfFinal = cut_stf_by_mask(stfLayer, spotMask, mapLayer);
resultFinal = run_final_reoptimization(dijFinal, cst, plnArc, wFinalInit, cfg);

%% 10. Post-filter Spot Removal Baseline (Method 4)
resultPostFilter = run_post_filter_baseline(dijArc, cst, plnArc, resultDenseArc.w, cfg);

%% 11. Robust Evaluation
if cfg.robust.enabled
    robustEval = evaluate_robust_scenarios(ct, cst, stfFinal, dijFinal, plnArc, resultFinal.w, cfg);
else
    robustEval = [];
end

%% 12. Compute All Metrics
results.impt       = evaluate_method(resultIMPT, dijIMPT, map, machine, cfg, 'IMPT');
results.denseArc   = evaluate_method(resultDenseArc, dijArc, map, machine, cfg, 'Dense Arc');
results.elospatLike = evaluate_method(resultELO, dijArc, map, machine, cfg, 'ELO-SPAT-like');
results.postFilter = evaluate_method(resultPostFilter, dijArc, map, machine, cfg, 'Post-filter');
results.rmsElospat = evaluate_method(resultFinal, dijFinal, mapFinal, machine, cfg, 'RMS-ELO-SPAT');
results.robustEval = robustEval;
results.cfg = cfg;
results.machine = machine;

%% 13. Export
T = compare_methods_table(results);
export_latex_table(T, fullfile(cfg.outputDir, 'comparison_table.tex'));
export_tradeoff_data(results, fullfile(cfg.outputDir, 'tradeoff_data.mat'));
save(fullfile(cfg.outputDir, 'rms_elospat_results.mat'), 'results');

fprintf('RMS-ELO-SPAT pipeline complete. Results saved to: %s\n', cfg.outputDir);
end
```


---

## 4. Key Algorithms

### 4.1 Layer Sparsification (Wrapper Mode)

```
Input: dij, w_init, groups, cpGroups, cfg, machine
Output: w_final, layerMask

w = w_init
best_score = Inf
best_w = w
D = dij.physicalDose{1}  // nominal scenario dose influence matrix

FOR iter = 1 to cfg.wrapper.maxLayerIterations:
    // Compute group norms
    FOR each group g:
        groupNorm(g) = ||w(groups(g).ix)||_2

    // Compute layer dose importance (heuristic wrapper score, NOT the
    // IPOPT dose objective). This is a proxy for how much dose quality
    // would degrade if this layer were removed.
    FOR each group g:
        ix_g = groups(g).ix
        doseContribution(g) = norm(D(:, ix_g) * w(ix_g), 2)

    // Compute estimated BDT saving per layer removal.
    // bdtSaving(g) = T_current - T_without_group_g
    // This estimates how much delivery time would decrease if group g
    // were removed (accounts for energy switching, spot count, etc.)
    T_current = estimate_bdt(extract_delivery_sequence(w, ...), machine).T_total
    FOR each active group g:
        w_without_g = w; w_without_g(groups(g).ix) = 0
        T_without_g = estimate_bdt(extract_delivery_sequence(w_without_g, ...), machine).T_total
        bdtSaving(g) = T_current - T_without_g

    // Compute removal score: layers with LOW dose importance and HIGH
    // BDT saving are best candidates for removal.
    // removalScore(g) = doseContribution(g) / (bdtSaving(g) + eps)
    // Low removalScore => good candidate for removal.
    FOR each active group g:
        removalScore(g) = doseContribution(g) / (bdtSaving(g) + eps)

    // Identify bottom-K layers for removal by removalScore
    K = ceil(nActiveGroups * cfg.wrapper.removalFraction)  // e.g., 10%
    candidates = bottom K layers by removalScore

    // Build mask: keep all except candidates
    layerMask = true(nGroups, 1)
    layerMask(candidates) = false

    // Expand to bixel-level mask and reduce
    bixelMask = expand_layer_mask_to_bixel_mask(layerMask, groups, N)
    [dijR, mapR, wR] = cut_dij_by_mask(dij, bixelMask, map, w)
    result = matRad_fluenceOptimization(dijR, cst, pln, wR)
    w_reduced = result.w

    // Compute composite score on reduced plan
    score = composite_score(w_reduced, ...)

    // Track best
    IF score < best_score:
        best_score = score
        best_w = w_reduced
        best_mask = layerMask

    // Check convergence
    IF |score - prev_score| / |prev_score| < cfg.wrapper.relTol AND iter >= 3:
        BREAK

    // Update for next iteration
    w = embed_reduced_to_full(w_reduced, bixelMask)
    prev_score = score

RETURN best_w, best_mask
```

**Note on dose importance:** `doseContribution(g) = norm(D(:, ix_g) * w(ix_g), 2)` is a heuristic wrapper-level score measuring the L2 norm of the dose deposited by group g. It is NOT the matRad/IPOPT internal dose objective. Alternative formulations (e.g., normalized by total weight: `doseContribution(g) / (sum(abs(w(ix_g))) + eps)`) may be explored as hyperparameter choices.

**Note on BDT saving:** Computing `bdtSaving` for every group at every iteration is O(nGroups) BDT evaluations. For large problems, an approximate version may cache the current sequence and estimate marginal savings analytically (removing one layer reduces T_energy by its switching contribution and T_spot by its spot-switch contribution).

### 4.2 Spot Filtering (Wrapper Mode)

```
Input: dij, w_init, map, groups, cfg, machine
Output: w_final, spotMask

w = w_init
prev_ctvD95 = compute_ctvD95(dij, w, cst)
D = dij.physicalDose{1}  // nominal scenario

FOR iter = 1 to cfg.wrapper.maxSpotIterations:
    // Compute per-spot dose importance (heuristic, not IPOPT objective)
    FOR each spot i where w(i) > 0:
        spotDoseImportance(i) = norm(D(:, i) * w(i), 2)

    // Compute per-spot BDT saving estimate
    // Removing a spot reduces T_spot (one fewer transition) and
    // T_irradiation (less weight to deliver). If removing the last spot
    // in a layer, also saves T_energy for that layer's switching.
    FOR each spot i where w(i) > 0:
        spotBdtSaving(i) = machine.spot.switchTime  // at minimum
                         + w(i) / machine.irradiation.muRate
        // If spot is the only active spot in its layer:
        IF is_last_spot_in_layer(i, w, groups):
            spotBdtSaving(i) += estimated_layer_switch_saving(i, ...)

    // Compute spot removal score: low importance + high BDT saving = remove
    FOR each spot i where w(i) > 0:
        spotRemovalScore(i) = spotDoseImportance(i) / (spotBdtSaving(i) + eps)

    // Identify candidates: sub-MMU spots OR lowest removal scores
    subMMU = (w > 0) & (w < cfg.mmuThreshold)
    lowScore = spotRemovalScore < percentile(spotRemovalScore(w>0), cfg.wrapper.spotRemovalPercentile)
    candidates = subMMU | lowScore

    // Build spot mask
    spotMask = ~candidates

    // Reduce and re-optimize
    [dijR, mapR, wR] = cut_dij_by_mask(dij, spotMask, map, w)
    result = matRad_fluenceOptimization(dijR, cst, pln, wR)

    // Check dose degradation
    new_ctvD95 = compute_ctvD95(dijR, result.w, cst)
    IF (prev_ctvD95 - new_ctvD95) / prev_ctvD95 > cfg.evaluation.ctvD95TolerancePercent/100:
        WARNING: dose degradation exceeded tolerance
        REVERT to previous w
        BREAK

    w = result.w
    prev_ctvD95 = new_ctvD95

    // Check convergence (estimated BDT change)
    bdt = estimate_bdt(extract_delivery_sequence(w, mapR, groups, cfg), machine)
    IF |bdt.T_total - prev_bdt| / prev_bdt < cfg.wrapper.relTol:
        BREAK

RETURN w, spotMask
```

**Note on spot removal score:** The spot removal score `spotDoseImportance / (spotBdtSaving + eps)` implements a Pareto-like tradeoff: spots that contribute little dose but save significant delivery time are removed first. Sub-MMU spots are always candidates regardless of their score, since they are not machine-deliverable.

### 4.3 BDT Estimation Algorithm

```
Input: seq (delivery sequence), machine (timing model)
Output: metrics (BDT decomposition)

T_energy = 0
T_spot = 0
T_irradiation = 0
T_gantry = 0
T_buffer = 0
nSwitchUp = 0

FOR k = 1 to length(seq):
    // Energy switching (from previous layer)
    IF k > 1:
        dE = seq(k).energy - seq(k-1).energy
        IF dE > 0:  // switch up
            T_energy += machine.energy.upBase + machine.energy.perMeVUp * dE
            nSwitchUp += 1
        ELSEIF dE < 0:  // switch down
            T_energy += machine.energy.downBase + machine.energy.perMeVDown * |dE|

        // Gantry motion (only between different CPs)
        IF seq(k).cpIndex != seq(k-1).cpIndex:
            dTheta = |seq(k).gantryAngle - seq(k-1).gantryAngle|
            IF dTheta > machine.gantry.toleranceWindow:
                T_gantry += estimate_gantry_time(dTheta, machine)

    // Spot switching within this layer
    nSpots = seq(k).numActiveSpots
    T_spot += max(nSpots - 1, 0) * machine.spot.switchTime

    // Irradiation time
    T_irradiation += seq(k).layerWeightSum / machine.irradiation.muRate

    // Buffer
    T_buffer += machine.buffer.perLayer

// Add per-CP buffer
nCP = numel(unique([seq.cpIndex]))
T_buffer += nCP * machine.buffer.perControlPoint

metrics.T_energy = T_energy
metrics.T_spot = T_spot
metrics.T_irradiation = T_irradiation
metrics.T_gantry = T_gantry
metrics.T_buffer = T_buffer
metrics.T_total = T_energy + T_spot + T_irradiation + T_gantry + T_buffer
metrics.nSwitchUp = nSwitchUp
```


---

## Data Models

### 5.1 Core Data Structures

#### Bixel Map Entry

```matlab
map(i).globalIx     % int: 1-based sequential index
map(i).beamNum      % int: beam number (= control point)
map(i).rayNum       % int: ray number within beam
map(i).bixelNum     % int: bixel number within ray
map(i).gantryAngle  % double: gantry angle in degrees
map(i).energy       % double: physical energy in MeV
map(i).cpIndex      % int: control point index (= beamNum)
map(i).energyIndex  % int: energy layer index within CP
```

#### Energy-Layer Group

```matlab
groups(g).cpIndex      % int: control point index
groups(g).energyIndex  % int: energy layer index within CP
groups(g).energy       % double: physical energy (MeV)
groups(g).ix           % int vector: global bixel indices in this group
groups(g).nSpots       % int: number of spots in group
```

#### Control-Point Group

```matlab
cpGroups(b).cpIndex    % int: control point index
cpGroups(b).ix         % int vector: all bixel indices at this CP
cpGroups(b).gantryAngle % double: gantry angle
```

#### Delivery Sequence Entry

```matlab
seq(k).cpIndex        % int: control point index
seq(k).gantryAngle    % double: gantry angle (degrees)
seq(k).energy         % double: energy (MeV)
seq(k).energyIndex    % int: energy layer index
seq(k).spotIndices    % int vector: active spot global indices
seq(k).spotWeights    % double vector: active spot weights
seq(k).numActiveSpots % int: count of active spots
seq(k).layerWeightSum % double: sum of spot weights in layer
seq(k).prevEnergy     % double: energy of previous layer (NaN if first)
seq(k).nextEnergy     % double: energy of next layer (NaN if last)
```

#### BDT Metrics

```matlab
metrics.T_energy          % double: total energy switching time (s)
metrics.T_spot            % double: total spot switching time (s)
metrics.T_irradiation     % double: total irradiation time (s)
metrics.T_gantry          % double: total gantry rotation time (s)
metrics.T_buffer          % double: total buffer time (s)
metrics.T_total           % double: sum of all components (s)
metrics.nActiveLayers     % int: number of active energy layers
metrics.nActiveSpots      % int: number of active spots
metrics.nSwitchUp         % int: number of switch-up events
metrics.totalSwitchMagnitude % double: sum of |dE| across all switches
```

#### Method Result

```matlab
result.w              % double vector: optimized spot weights
result.doseCube       % 3D array: dose distribution
result.doseMetrics    % struct: DVH metrics
result.delivMetrics   % struct: BDT decomposition
result.solverMetrics  % struct: optimizer performance
result.robustMetrics  % struct: worst-case metrics (if available)
result.methodName     % string: method identifier
result.cfg            % struct: configuration used
```

### 5.2 File I/O

| File | Format | Content |
|---|---|---|
| `results/rms_elospat_results.mat` | MAT | Full results struct |
| `results/comparison_table.tex` | LaTeX | Method comparison table |
| `results/tradeoff_data.mat` | MAT | Dose vs BDT pairs |
| `results/sweep_<params>.mat` | MAT | Per-combination sweep results |
| `results/intermediate_iter<N>.mat` | MAT | Per-iteration snapshots |

---

## 6. Integrated Mode Design (Mode B - Future)

### 6.1 Custom Optimization Problem Class

```matlab
classdef matRad_OptimizationProblemELOSPAT < matRad_OptimizationProblem
    properties
        groups      % energy-layer group index
        cpGroups    % control-point group index
        map         % bixel mapping
        machine     % machine timing model
        cfg         % RMS-ELO-SPAT configuration
    end

    methods
        function obj = matRad_OptimizationProblemELOSPAT(...)
            % Constructor: store groups, map, machine, cfg
        end

        function f = objectiveFunction(obj, w)
            % f = f_dose + lambdaEL*R_EL + lambdaB*B + lambdaES*Phi_ES
            %     + lambdaS*R_spot + lambdaT*T_smooth
        end

        function g = objectiveGradient(obj, w)
            % g = g_dose + lambdaEL*g_EL + lambdaB*g_B + lambdaES*g_ES
            %     + lambdaS*g_spot + lambdaT*g_T_smooth
        end
    end
end
```

### 6.2 Gradient Validation Protocol

Before enabling Integrated Mode, each penalty gradient must pass:

1. Generate random w vector (N elements, uniform [0, 1])
2. For each component (R_EL, B, Phi_ES, R_spot, T_smooth):
   a. Compute analytical gradient g_analytic
   b. For 100 randomly sampled indices i:
      - Compute finite-difference: g_fd(i) = (f(w+h*e_i) - f(w-h*e_i)) / (2h)
      - h = max(1e-7, 1e-7 * |w_i|)
   c. Relative error = |g_analytic - g_fd| / max(|g_analytic|, |g_fd|, 1e-10)
   d. PASS if 95% of samples have relative error < 1e-4

---

## Error Handling

### Configuration Errors
- Missing or invalid config fields: `validate_config` throws descriptive error before optimization starts
- Invalid arc spacing (doesn't divide 360): rejected with supported values listed
- Patient file not found: error with file path

### Optimization Errors
- matRad convergence failure: result stored with `converged = false`, iteration count, last objective
- Memory exceeded (large dij): error recommending coarser spacing or smaller dose grid
- All spots removed by MMU filter: error, revert to pre-threshold solution
- CTV D95 degradation > tolerance: warning, revert to previous solution

### Data Integrity Errors
- Bixel unmapped or multiply-mapped: `validate_mapping` throws error with affected indices
- Reduced dij dimension mismatch: error before re-optimization
- Empty control point after reduction: warning (or error if barrier not enabled)

### Graceful Degradation
- Parameter sweep: single combination failure logged, remaining combinations continue
- Robust evaluation unavailable: skipped with warning, nominal metrics still reported
- Irradiation time unit unknown: reported as approximate

---

## Correctness Properties

### Property 1: Mapping Bijectivity
Every bixel maps to exactly one (CP, energy, spot) triplet. Verified by `validate_mapping`. `numel(map) == dij.totalNumOfBixels` and all globalIx values are unique consecutive integers.

**Validates: Requirements 3.1, 3.4, 3.5**

### Property 2: Gradient Consistency
All analytical gradients must agree with finite-difference approximation to relative error < 1e-4 for 95% of sampled variables. Tested before enabling Integrated Mode.

**Validates: Requirements 4.2, 5.2, 10.2**

### Property 3: Reduction Consistency
After `cut_dij_by_mask`, forward dose calculation `dij_reduced * w_reduced` must equal the corresponding subset of `dij_full * w_full` to machine precision.

**Validates: Requirements 11.2, 11.4, 11.5**

### Property 4: BDT Decomposition Additivity
`T_total == T_energy + T_spot + T_irradiation + T_gantry + T_buffer` (exact floating-point equality).

**Validates: Requirements 8.1, 8.12**

### Property 5: Monotonic Layer Ordering
When descending order is selected, `seq(k).energy >= seq(k+1).energy` within same control point for all consecutive layers.

**Validates: Requirements 7.4, 7.5**

### Property 6: Dose Preservation
After re-optimization on reduced problem, dose must be recomputable from reduced dij and new weights without NaN or Inf values.

**Validates: Requirements 11.2, 11.5, 15.2**

### Property 7: Reproducibility
Same cfg + same patient + same random seed produces identical or numerically equivalent results (within configured solver tolerances) on the same software/hardware environment. IPOPT version, MATLAB version, BLAS implementation, and multithreading settings may affect exact floating-point reproducibility across different platforms.

**Validates: Requirements 17.1, 17.6**

---

## Testing Strategy

### 7.1 Unit Tests

| Test | Module | Validates |
|---|---|---|
| `test_mapping_integrity` | mapping | Req 3: bijective mapping |
| `test_l21half_gradient` | objective | Req 4: gradient correctness |
| `test_barrier_gradient` | objective | Req 5: barrier gradient |
| `test_spot_sparsity_gradient` | objective | Req 10: spot gradient |
| `test_bdt_estimator` | delivery | Req 8: BDT decomposition |
| `test_dij_reduction` | reduction | Req 11: consistent reduction |
| `test_dij_reduction_multiscenario` | reduction | Req 11: cell-array/multi-scenario field handling |
| `test_gantry_time` | delivery | Req 8: trapezoidal profile |
| `test_sequencing_penalty` | objective | Req 6: zero when monotonic |
| `test_expand_layer_mask` | reduction | Verify layer→bixel mask expansion correctness |
| `test_mask_dimension_safety` | reduction | Verify layerMask (nGroups) cannot be passed where bixelMask (N) is expected |
| `test_bdt_aware_removal_score` | planning | Verify removal score ranks low-importance/high-BDT-saving layers first |

### 7.2 Integration Tests

| Test | Validates |
|---|---|
| PROSTATE 2-beam IMPT | Req 1: baseline runs |
| PROSTATE 20° arc | Req 2: arc geometry |
| Full pipeline 20° | Req 15: wrapper loop |
| Robust evaluation | Req 12: scenario metrics |

### 7.3 Gradient Test Template

```matlab
function pass = test_gradient(func, w, h, nSamples, relTolPass)
    [~, g_analytic] = func(w);
    N = numel(w);
    indices = randperm(N, min(nSamples, N));
    errors = zeros(numel(indices), 1);

    for k = 1:numel(indices)
        i = indices(k);
        wp = w; wp(i) = wp(i) + h;
        wm = w; wm(i) = wm(i) - h;
        g_fd = (func(wp) - func(wm)) / (2*h);
        errors(k) = abs(g_analytic(i) - g_fd) / max(abs(g_analytic(i)), abs(g_fd), 1e-10);
    end

    passRate = mean(errors < relTolPass);
    pass = passRate >= 0.95;
    fprintf('Gradient test: %.1f%% pass (threshold 95%%)\n', passRate*100);
end
```


---

## 8. Requirements Traceability

| Requirement | Design Component(s) |
|---|---|
| Req 1: Baseline IMPT | `planning/build_impt_plan.m` |
| Req 2: Discrete Arc | `planning/build_discrete_arc_plan.m` |
| Req 3: Bixel Mapping | `mapping/map_bixel_to_cp_energy_spot.m`, `build_group_index.m`, `validate_mapping.m` |
| Req 4: L2,1/2 Sparsity | `objective/l21half_group_penalty.m` |
| Req 5: CP Barrier | `objective/control_point_barrier.m` |
| Req 6: Sequencing Penalty | `objective/energy_sequencing_penalty.m` |
| Req 7: Delivery Sequence | `delivery/extract_delivery_sequence.m`, `order_layers_within_control_point.m` |
| Req 8: BDT Estimator | `delivery/estimate_bdt.m`, `estimate_energy_switch_time.m`, `estimate_gantry_time.m` |
| Req 9: Smooth BDT | Integrated Mode class (future) |
| Req 10: Spot Sparsity/MMU | `objective/spot_sparsity_penalty.m`, `planning/run_spot_sparse_stage.m` |
| Req 11: Dij Reduction | `reduction/cut_dij_by_mask.m`, `cut_stf_by_mask.m` |
| Req 12: Robustness | `robustness/configure_robust_scenarios.m`, `evaluate_robust_scenarios.m` |
| Req 13: Method Comparison | `evaluation/compare_methods_table.m` |
| Req 14: Metrics | `evaluation/compute_dvh_metrics.m`, `compute_delivery_metrics.m`, `compute_solver_metrics.m` |
| Req 15: Wrapper Loop | `planning/run_elospat_like_stage.m`, `run_spot_sparse_stage.m`, orchestrator |
| Req 16: Integrated Mode | `@matRad_OptimizationProblemELOSPAT` (future) |
| Req 17: Configuration | `config/default_config.m`, `validate_config.m` |
| Req 18: Safety Disclaimer | All report outputs, README.md |

---

## 9. Implementation Phases

### Phase 0: Infrastructure (MVP Foundation)
- `config/default_config.m`
- `config/default_machine_model.m`
- `config/validate_config.m`
- `planning/build_impt_plan.m`
- Run matRad baseline on PROSTATE.mat
- Verify: resultGUI.w, .RBExDose exist

### Phase 1: Arc Geometry + Mapping
- `planning/build_discrete_arc_plan.m`
- `mapping/map_bixel_to_cp_energy_spot.m`
- `mapping/build_group_index.m`
- `mapping/validate_mapping.m`
- `planning/run_dense_arc_baseline.m`
- Verify: mapping integrity, group counts

### Phase 2: Objective Functions
- `objective/l21half_group_penalty.m`
- `objective/control_point_barrier.m`
- `objective/energy_sequencing_penalty.m`
- `objective/spot_sparsity_penalty.m`
- `objective/composite_score.m`
- Verify: gradient tests pass

### Phase 3: Delivery Model
- `delivery/extract_delivery_sequence.m`
- `delivery/order_layers_within_control_point.m`
- `delivery/estimate_bdt.m` (+ sub-functions)
- `delivery/validate_gantry_constraints.m`
- Verify: BDT decomposition on known sequence

### Phase 4: Reduction + Re-optimization
- `reduction/build_layer_mask.m`
- `reduction/build_spot_mask.m`
- `reduction/cut_dij_by_mask.m`
- `reduction/cut_stf_by_mask.m`
- `reduction/reoptimize_reduced_problem.m`
- Verify: reduced dij produces valid dose

### Phase 5: Wrapper Loop Integration
- `planning/run_elospat_like_stage.m`
- `planning/run_spot_sparse_stage.m`
- `planning/run_final_reoptimization.m`
- `run_rms_elospat_main.m`
- Verify: full pipeline runs end-to-end on 20° arc

### Phase 6: Robustness
- `robustness/configure_robust_scenarios.m`
- `robustness/configure_robust_objectives.m`
- `robustness/evaluate_robust_scenarios.m`
- `robustness/compute_worst_case_metrics.m`
- Verify: worst-case metrics computed

### Phase 7: Evaluation + Export
- `evaluation/compute_dvh_metrics.m`
- `evaluation/compute_delivery_metrics.m`
- `evaluation/compute_solver_metrics.m`
- `evaluation/compare_methods_table.m`
- `evaluation/export_latex_table.m`
- `evaluation/export_tradeoff_data.m`
- Verify: LaTeX table generated, tradeoff data exported

### Phase 8 (Future): Integrated Mode
- `@matRad_OptimizationProblemELOSPAT` class
- Gradient validation suite
- Smooth BDT approximation

---

## 10. Design Decisions and Rationale

| Decision | Rationale |
|---|---|
| External module in `userdata/` | Avoid breaking matRad core; easy to version separately |
| Wrapper Mode first | Lower risk, faster iteration, sufficient for proof-of-concept paper |
| Serial BDT model | Conservative estimate; synchronized model is future work |
| Descending energy order default | Reduces switch-up events (high→low is faster) |
| 2-stage loop (layer then spot) | Decouples concerns; layer reduction first reduces problem size |
| Threshold + reopt for MMU | Avoids forcing all spots nonzero via lower bounds |
| Composite score (not IPOPT obj) | Wrapper mode uses matRad's optimizer internally; outer score guides filtering |
| Gradient tests before Integrated Mode | IPOPT requires accurate gradients; non-smooth terms cause failure |
| Single cfg struct | Reproducibility; easy parameter sweeps; single save point |
| matRad functions called, not modified | Maintainability; can update matRad independently |

---

## 11. Limitations and Future Work

1. **Serial BDT model** — Real machines may overlap gantry motion with irradiation. Future: synchronized timeline model.
2. **No interplay effects** — Motion interplay during delivery not modeled. Future: 4D dose calculation.
3. **Proxy machine parameters** — Default timing values are placeholders. Future: validate against vendor logs.
4. **Wrapper mode is heuristic** — Not a single unified optimization. Future: Integrated Mode with validated gradients.
5. **Memory constraints** — Fine arc spacing (2°) may exceed available RAM. Future: sparse dij storage or out-of-core computation.
6. **Single patient initially** — PROSTATE.mat only. Future: multi-site evaluation (H&N, lung).


---

## Revision Notes

Changes applied to this document to improve implementation correctness, mathematical consistency, and alignment with the literature gap:

1. **Dependency table fixed:** Replaced hard-coded `matRad_calcDoseDirect` with version-aware fallback logic (`matRad_calcDoseForward` / `matRad_calcCubes` / manual `dij * w`). Exact function selected at runtime after inspecting installed matRad version.

2. **`cut_dij_by_mask` expanded:** Now explicitly handles all dose-influence-like fields (`physicalDose`, `mAlphaDose`, `mSqrtBetaDose`, `mLETDose`) as cell arrays iterated over scenarios, plus per-bixel vector fields (`beamNum`, `rayNum`, `bixelNum`, `minMU`, `maxMU`, `numParticlesPerMU`). Missing fields are skipped safely.

3. **Layer mask → bixel mask dimension mismatch resolved:** Added `expand_layer_mask_to_bixel_mask(layerMask, groups, totalNumOfBixels)` function. Orchestrator pseudocode updated to call this before `cut_dij_by_mask`. Documentation explicitly warns that `layerMask` (nGroups) and `bixelMask` (totalNumOfBixels) are different dimensions.

4. **`run_elospat_like_stage` clarified:** Explicitly stated this is an "ELO-SPAT-inspired wrapper approximation", NOT an exact reproduction of the original ELO-SPAT unified solver. Exact integrated optimization deferred to Future Integrated Mode.

5. **Layer importance / dose contribution proxy defined precisely:** Replaced vague `dij_columns contribution` with explicit formula: `doseContribution(g) = norm(D(:, ix_g) * w(ix_g), 2)` where `D = dij.physicalDose{1}`. Clearly labeled as heuristic wrapper score, not IPOPT dose objective.

6. **BDT-saving integrated into layer/spot selection:** Layer removal now uses `removalScore(g) = doseContribution(g) / (bdtSaving(g) + eps)` where `bdtSaving(g) = T_current - T_without_group_g`. Spot filtering uses analogous `spotRemovalScore`. This directly implements the RMS-ELO-SPAT novelty of optimizing for estimated delivery time reduction, not just layer count.

7. **Robust evaluation interface fixed:** `evaluate_robust_scenarios` now accepts `stf` as parameter. Documentation clarifies when `stf` is needed vs. when multi-scenario `dij` cells suffice. Distinguishes robust optimization from robust evaluation.

8. **Reproducibility claim corrected:** Changed from "identical results" to "identical or numerically equivalent within configured tolerances on the same software/hardware environment." Added note about IPOPT/MATLAB/BLAS/platform variability.

9. **BDT terminology made consistent:** `estimate_bdt.m` docstring explicitly states output is an "estimated machine-specific BDT proxy" unless `machine.validation.status == 'validated'`. Term "true BDT" reserved only for validated machine-specific data.

10. **Testing strategy expanded:** Added tests for `expand_layer_mask_to_bixel_mask`, multi-scenario `cut_dij_by_mask`, BDT-aware removal score, and mask dimension safety (preventing layerMask/bixelMask confusion).
