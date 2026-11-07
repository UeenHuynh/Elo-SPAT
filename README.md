# RMS-ELO-SPAT

MATLAB research code for discrete proton arc planning with energy-layer and spot sparsity, dose-quality checks, and an estimated beam delivery time (BDT) model. The implementation builds on [matRad](https://github.com/e0404/matRad) without changing matRad core files.

> **Research use only.** This repository does not provide a clinical treatment-planning system. BDT is a configurable model estimate, not a measured treatment-machine delivery time. Plans and metrics here must not be used for patient care.

## Start reading here

| If you want to understand… | Read |
| --- | --- |
| The intended behavior and design | [Requirements](.kiro/specs/rms-elospat/requirements.md) and [design](.kiro/specs/rms-elospat/design.md) |
| The top-level wrapper workflow | [Main orchestrator](userdata/elospat_time/run_rms_elospat_main.m) and [default configuration](userdata/elospat_time/config/default_config.m) |
| How a discrete arc and dose problem are built | [`planning/`](userdata/elospat_time/planning/) and [`mapping/`](userdata/elospat_time/mapping/) |
| The integrated optimization path | [`integrated/`](userdata/elospat_time/integrated/) |
| Stage 2 support-selection experiments | [`experiments/`](userdata/elospat_time/experiments/) and their [`slurm/`](userdata/elospat_time/slurm/) launch scripts |
| Dose, robustness, and delivery metrics | [`evaluation/`](userdata/elospat_time/evaluation/), [`robustness/`](userdata/elospat_time/robustness/), and [`delivery/`](userdata/elospat_time/delivery/) |
| Small checks and smoke tests | [`tests/`](userdata/elospat_time/tests/) |

The wrapper's numbered stages in `run_rms_elospat_main.m` describe execution order. The **integrated Stage 1 / Stage 2** terminology used in experiment scripts refers to layer selection followed by spot/support optimization; it is a separate research path. Check the configuration and experiment entry point before comparing outputs from these paths.

## Source layout

```text
userdata/elospat_time/
  config/       Defaults, machine timing proxy, validation
  planning/     IMPT baseline, discrete arc construction, wrapper stages
  mapping/      Spot, energy-layer, and control-point indexing
  objective/    Sparsity and delivery-related objective terms
  integrated/   Integrated optimizer and numerical solvers
  experiments/  Isolated Stage 2 method studies and repair workflows
  reduction/    Restricted dose-influence and support operations
  robustness/   Uncertainty scenarios and robust evaluation helpers
  delivery/     Delivery sequence and estimated BDT
  evaluation/   DVH, solver, and comparison metrics
  tests/        MATLAB smoke and unit tests
  slurm/        Batch launch scripts
```

## Running locally

You need MATLAB, a compatible matRad installation on the MATLAB path, and patient data available to matRad. The default configuration expects `PROSTATE.mat` and a generic proton machine. Large full-arc runs can require substantial RAM and time; start with a small test or the baseline workflow.

```matlab
% Run from the repository root after initializing matRad.
cd('userdata/elospat_time')
check_matrad_path
run_rms_elospat_main
```

`run_rms_elospat_main(cfgOverride)` accepts a configuration override struct. See [`default_config.m`](userdata/elospat_time/config/default_config.m) for parameters. Experiment-specific jobs have their own `run_*.m` entry points and Slurm scripts; they may require prior checkpoints and should not be treated as standalone reproductions from a fresh clone.

## Data and results

This Git repository contains source code and launch scripts. Patient data, dose-influence matrices, checkpoints, generated `.mat` files, CSV tables, logs, figures, and PDFs are excluded. Local runs write generated artifacts under `userdata/elospat_time/results/` by default. A fresh clone therefore **does not contain the historical experiment results** or enough data to reproduce their numerical values immediately.

In the original research workspace, the current-result summaries are `CURRENT_RESULTS.md` and `KET_QUA_DEN_HIEN_TAI_2026-11-05.md`; the detailed Stage 2 ledger is `userdata/elospat_time/docs/STAGE2_COMPLETE_EXPERIMENT_LEDGER.md`. These reports and the underlying results are not part of this code-only upload. Solver completion, clinical-gate acceptance, numerical stationarity, and robust acceptance are separate checks; do not infer one from another.

## Terms

- **Spot:** One pencil-beam irradiation position and its weight.
- **Energy layer:** Spots delivered at one beam energy within a control point.
- **Control point (CP):** A sampled gantry angle in the discrete arc.
- **Support:** The spots and layers retained with nonzero weight.
- **BDT:** Estimated beam delivery time from the configurable timing model.
- **`dij`:** matRad dose-influence data used to calculate dose from spot weights.
