function cfg = default_config()
% default_config  Returns the default RMS-ELO-SPAT configuration struct.
%
% All tunable parameters are defined here with documented defaults and
% valid ranges. Pass the returned struct directly to the main orchestrator
% or to validate_config() before use.
%
% OUTPUT
%   cfg - struct with all configuration fields initialised to defaults
%
% USAGE
%   cfg = default_config();
%   cfg.arcSpacing = 10;          % override individual fields as needed
%   validate_config(cfg, machine);
%
% VALID MODES
%   cfg.mode = 'wrapper'     - Wrapper Prototype Mode (Mode A, required)
%   cfg.mode = 'integrated'  - MATLAB-native integrated RMS-ELO-SPAT
%
% =========================================================================
% RESEARCH-ONLY DISCLAIMER:
%   This configuration is part of RMS-ELO-SPAT, a research-only module.
%   Not for clinical use. Machine timing values are estimated proxies only.
% =========================================================================

% -------------------------------------------------------------------------
% Implementation mode
% -------------------------------------------------------------------------
% 'wrapper'    : Mode A — external iterative filtering around matRad.
%                matRad core is NOT modified. Required minimum viable mode.
% 'integrated' : Mode B — native MATLAB proximal optimization followed by
%                support-constrained matRad/IPOPT reoptimization.
cfg.mode = 'wrapper';

% -------------------------------------------------------------------------
% Patient and dataset
% -------------------------------------------------------------------------
% Full path to the patient .mat file, or a filename resolvable via
% MATLAB's path (e.g. 'PROSTATE.mat' if matRad data/ is on the path).
% validate_config() will verify the file exists before optimisation starts.
cfg.patientFile = 'PROSTATE.mat';

% -------------------------------------------------------------------------
% matRad plan settings
% -------------------------------------------------------------------------
cfg.radiationMode = 'protons';   % radiation modality passed to matRad pln
cfg.machine       = 'Generic';   % matRad machine name (must match matRad DB)
cfg.bioModel      = 'constRBE';  % biological model ('none' | 'constRBE' | ...)

% -------------------------------------------------------------------------
% Dose calculation and discretization
% -------------------------------------------------------------------------
cfg.numOfFractions = 1;          % Number of fractions (research default)
cfg.bixelWidth = 5;              % Bixel width in mm (matRad default)

% Dose grid resolution (defaults to CT resolution if not set)
% Can be coarser than CT for faster computation
cfg.doseGrid = struct();
cfg.doseGrid.x = [];  % empty = use ct.resolution.x
cfg.doseGrid.y = [];  % empty = use ct.resolution.y
cfg.doseGrid.z = [];  % empty = use ct.resolution.z

% Isocenter (optional, matRad auto-computes if empty)
cfg.isoCenter = [];  % [x, y, z] in mm, or [] for auto

% -------------------------------------------------------------------------
% Geometry — IMPT baseline
% -------------------------------------------------------------------------
% Gantry angles (degrees) for the IMPT multi-beam baseline plan.
% Must be a 1-D vector with 2–5 elements. Default: [90, 270].
cfg.imptAngles = [90, 270];

% -------------------------------------------------------------------------
% Geometry — discrete proton arc
% -------------------------------------------------------------------------
% Angular spacing between consecutive control points (degrees).
% Must divide 360 exactly (e.g. 20, 10, 5, 4, 3, 2).
% Couch angles are fixed at 0 for all control points.
%
% Valid spacing values: any divisor of 360.
% Recommended starting point: 20 (18 CPs, low memory).
% Valid range: 1 <= arcSpacing <= 180.
cfg.arcSpacing  = 20;
cfg.couchAngles = 0;   % scalar — applied uniformly to all control points

% -------------------------------------------------------------------------
% Penalty weights  (lambda_*)
% -------------------------------------------------------------------------
% All weights default to 0 so the nominal matRad dose objective drives
% optimisation by default. Set non-zero values to activate each term.
%
%   lambdaEL  : L_{2,1/2} energy-layer group sparsity  [0, 1e-2]
%   lambdaB   : log-barrier CP activity preservation   [0, 1e-2]
%   lambdaES  : energy-sequencing penalty              [0, 1e-2]
%   lambdaS   : spot sparsity                          [0, 1e-2]
%   lambdaT   : estimated BDT proxy term               [0, 1e-2]
cfg.lambdaEL = 0;
cfg.lambdaB  = 0;
cfg.lambdaES = 0;
cfg.lambdaS  = 0;
cfg.lambdaT  = 0;

% -------------------------------------------------------------------------
% Smoothing / regularisation parameters  (epsilon_*)
% -------------------------------------------------------------------------
% Numerical smoothing constants that prevent singularities in gradients.
%
%   epsEL : group sparsity smoothing    valid range [1e-10, 1e-6]
%   epsB  : log-barrier smoothing       valid range [1e-10, 1e-6]
%   epsZ  : activity soft-step (seqPen) valid range [1e-8, 1e-4]
%   epsS  : spot sparsity smoothing     valid range [1e-5, 1e-1]
cfg.epsEL = 1e-8;
cfg.epsB  = 1e-8;
cfg.epsZ  = 1e-6;
cfg.epsS  = 1e-3;

% -------------------------------------------------------------------------
% Activity / deliverability thresholds
% -------------------------------------------------------------------------
% layerThreshold : group L2-norm below which an energy layer is considered
%                  inactive for delivery evaluation and scoring.
%                  Valid range: [1e-6, 1e-4].
% spotThreshold  : individual spot weight below which a spot is considered
%                  inactive for delivery sequence extraction.
%                  Valid range: [1e-8, 1e-4].
% mmuThreshold   : minimum monitor-unit threshold (x_min).
%                  Spots satisfying 0 < w < mmuThreshold are sub-MMU and
%                  will be zeroed and removed before re-optimisation.
%                  Valid range: [0, Inf) — set to 0 to disable MMU filtering.
cfg.layerThreshold = 1e-5;
cfg.spotThreshold  = 1e-6;
cfg.mmuThreshold   = 1e-5;

% -------------------------------------------------------------------------
% Robustness settings
% -------------------------------------------------------------------------
% enabled      : false = nominal only; true = enable robust scenarios.
% mode         : matRad robust optimisation paradigm.
%                Supported: 'COWC' (Composite OWC), 'OWC' (Objective-Wise WC).
% setupMM      : 1x3 vector of setup uncertainty [x,y,z] in mm.
%                Default: ±3 mm in each axis.
% rangePercent : range uncertainty in percent.
%                Default: 3.5% (if supported by the selected scenario model).
cfg.robust.enabled      = false;
cfg.robust.mode         = 'COWC';
cfg.robust.setupMM      = [3, 3, 3];
cfg.robust.rangePercent = 3.5;

% -------------------------------------------------------------------------
% Wrapper loop settings
% -------------------------------------------------------------------------
% maxLayerIterations : maximum outer-loop iterations for layer filtering
%                      stage (run_elospat_like_stage).
% maxSpotIterations  : maximum outer-loop iterations for spot filtering
%                      stage (run_spot_sparse_stage).
% relTol             : relative change in composite score below which the
%                      outer loop is considered converged.
%                      Valid range: [1e-6, 1e-2].
cfg.wrapper.maxLayerIterations = 20;
cfg.wrapper.maxSpotIterations  = 20;
cfg.wrapper.relTol             = 1e-4;

% -------------------------------------------------------------------------
% Integrated MATLAB optimizer
% -------------------------------------------------------------------------
% componentStage controls incremental verification of one optimizer:
%   layer : F_matRad + lambdaEL*L_{2,1/2}
%   spot  : layer foundation + smooth-L0 spot activation
%   full  : all configured ES and delivery-time surrogate terms
cfg.integrated.componentStage = 'layer';
cfg.integrated.maxIterations = 200;
cfg.integrated.initialStep = 1e-3;
cfg.integrated.minimumStep = 1e-12;
cfg.integrated.maximumStep = 1;
cfg.integrated.backtrackFactor = 0.5;
cfg.integrated.backtrackGrowth = 0.8;
cfg.integrated.maxBacktracking = 30;
cfg.integrated.relativeStepTolerance = 1e-5;
cfg.integrated.relativeObjectiveTolerance = 1e-6;
% Convergence is declared from the normalized proximal-gradient mapping,
% not from the two relative-change diagnostics above. Dividing the proximal
% displacement by the accepted step prevents a small initial step from
% masquerading as stationarity.
cfg.integrated.stationarityTolerance = 1e-7;
% Disabled until a separate, scale-aware calibration is completed. The
% residual is always recorded as a diagnostic. Practical stopping is assessed
% by evaluate_numerical_stabilization at periodic BDT checkpoints.
cfg.integrated.stationarityStoppingEnabled = false;
cfg.integrated.minimumIterations = 10;
cfg.integrated.stagnationWindow = 5;
cfg.integrated.stagnationRelativeStepTolerance = 1e-12;
cfg.integrated.stagnationRelativeObjectiveTolerance = 1e-12;
cfg.integrated.sufficientDecreaseTolerance = 1e-12;
cfg.integrated.constraintPenalty = 1e3;
cfg.integrated.enforceFinalMinimumWeight = true;
cfg.integrated.verbose = true;
cfg.integrated.printEvery = 10;
% Number of accepted proximal updates for which smooth-L0 MM coefficients
% remain frozen. 1 reproduces the historical per-iteration MM algorithm.
% Values >1 implement controlled block/outer-MM schedules without changing
% the mathematical smooth-L0 or group-half objective.
cfg.integrated.spotMMUpdatePeriod = 1;

% -------------------------------------------------------------------------
% Evaluation tolerances
% -------------------------------------------------------------------------
% ctvD95TolerancePercent    : maximum allowed CTV D95 degradation (%) relative
%                             to the reference plan before a warning is raised
%                             and the pre-threshold plan is restored.
% robustD95TolerancePercent : same threshold applied to worst-case CTV D95
%                             in robust evaluation.
cfg.evaluation.ctvD95TolerancePercent    = 2.0;
cfg.evaluation.robustD95TolerancePercent = 2.0;

% -------------------------------------------------------------------------
% Output settings
% -------------------------------------------------------------------------
% outputDir        : directory where results (.mat, LaTeX tables) are saved.
%                    Relative paths are resolved from the working directory.
% saveIntermediate : if true, save intermediate per-iteration results to
%                    outputDir during the wrapper loop stages.
% captureConsoleOutput : if true, store raw MATLAB command-window output in
%                        console.log for long-running stages.
cfg.outputDir        = fullfile('userdata', 'elospat_time', 'results');
cfg.saveIntermediate = true;
cfg.captureConsoleOutput = true;

% -------------------------------------------------------------------------
% Runtime / stability settings
% -------------------------------------------------------------------------
% headless              : if true, suppress interactive figures and request
%                         matRad GUI/progress plotting to stay disabled.
% forceSoftwareOpenGL   : if true, request software OpenGL for the current
%                         MATLAB session. This can reduce graphics-driver
%                         instability on long optimization runs.
cfg.runtime.headless            = false;
cfg.runtime.forceSoftwareOpenGL = false;

% -------------------------------------------------------------------------
% Reproducibility
% -------------------------------------------------------------------------
% randomSeed : integer seed passed to rng() at the start of the main
%              orchestrator. Set to 0 to use the current random state.
cfg.randomSeed = 1;

end
