function machine = default_machine_model()
% default_machine_model  Returns the default machine timing model (proxy/unvalidated).
%
% This function provides a complete set of estimated timing parameters for
% a generic pencil-beam scanning (PBS) proton therapy machine. All values
% are order-of-magnitude estimates drawn from published literature and are
% NOT validated against any real machine.  They are labelled as proxy
% estimates via machine.validation.status = 'proxy_unvalidated'.
%
% OUTPUT
%   machine - struct with fields:
%               .energy      - energy-switching timing model
%               .spot        - intra-layer spot-switching timing
%               .irradiation - monitor-unit delivery rate
%               .gantry      - rotational motion profile
%               .buffer      - per-layer and per-CP overhead times
%               .validation  - provenance and validation status
%
% USAGE
%   machine = default_machine_model();
%   metrics = estimate_bdt(seq, machine);
%
% FIELD SUMMARY
%   machine.energy.upBase         1.5   s    (base cost, energy increase)
%   machine.energy.downBase       0.5   s    (base cost, energy decrease)
%   machine.energy.perMeVUp       0.02  s/MeV
%   machine.energy.perMeVDown     0.005 s/MeV
%   machine.spot.switchTime       0.002 s
%   machine.irradiation.muRate    1000  MU/s
%   machine.gantry.maxSpeed       6.0   deg/s
%   machine.gantry.maxAccel       1.0   deg/s^2
%   machine.gantry.toleranceWindow 0.5  deg
%   machine.buffer.perLayer       0.05  s
%   machine.buffer.perControlPoint 0.02 s
%
% =========================================================================
% RESEARCH-ONLY DISCLAIMER:
%   This machine model is part of RMS-ELO-SPAT, a research-only module.
%   Not for clinical use. All timing values are ESTIMATED PROXIES based
%   on published literature ranges. They have NOT been validated against
%   vendor specifications or measured delivery log data. BDT estimates
%   produced using this model should be treated as indicative only.
% =========================================================================

% -------------------------------------------------------------------------
% Energy switching — direction-dependent timing model
% -------------------------------------------------------------------------
% Changing the beam energy in a PBS system requires re-tuning the
% accelerator and beam transport elements.  This is the dominant source of
% delivery overhead in proton arc therapy.
%
% The switching time is modelled as a linear function of the energy
% difference |dE| in MeV, with separate base costs and per-MeV slopes for
% upward (energy increase) and downward (energy decrease) transitions:
%
%   T_switch = base + perMeV * |dE|
%
% upBase   (s)     : fixed overhead for any upward energy switch.
%                    Upward switches are slower because the accelerator
%                    must ramp up RF power and re-tune magnets from scratch.
%                    Typical reported range: 1–3 s.
% downBase (s)     : fixed overhead for any downward energy switch.
%                    Downward switches are faster on degrader-based systems
%                    (the degrader position is simply reduced).
%                    Typical reported range: 0.2–1 s.
% perMeVUp (s/MeV) : additional time per MeV of upward energy change.
% perMeVDown(s/MeV): additional time per MeV of downward energy change.
%
% NOTE: These values are proxy estimates for a synchrocyclotron-class
% machine with a range-degrader-based energy selection system. They may
% differ substantially for isochronous-cyclotron or synchrotron systems.

machine.energy.upBase    = 1.5;    % s     — base cost for energy increase
machine.energy.downBase  = 0.5;    % s     — base cost for energy decrease
machine.energy.perMeVUp  = 0.02;   % s/MeV — slope for upward switches
machine.energy.perMeVDown = 0.005; % s/MeV — slope for downward switches

% -------------------------------------------------------------------------
% Spot switching — intra-layer lateral scanning time
% -------------------------------------------------------------------------
% Within a single energy layer the beam is deflected laterally by fast
% scanning magnets to cover all pencil-beam spots. The time between
% consecutive spots (spot switching time) is much shorter than the
% energy-switching overhead.
%
% switchTime (s) : time between consecutive spot positions within one
%                  energy layer (magnetic deflection + settling).
%                  This is the per-spot transition cost, so the total
%                  intra-layer scanning overhead for a layer with N spots
%                  is max(N-1, 0) * switchTime.
%                  Typical values: 1–5 ms.  Default: 2 ms.

machine.spot.switchTime = 0.002;   % s (2 ms) — per spot lateral transition

% -------------------------------------------------------------------------
% Irradiation — monitor-unit delivery rate
% -------------------------------------------------------------------------
% Once the beam is positioned at a spot, the dose is deposited at a rate
% determined by the beam intensity expressed in monitor units per second.
%
% muRate     (MU/s)  : average dose-rate in monitor units per second.
%                      The irradiation time for a spot with weight w is
%                      T_irrad = w / muRate.
%                      Typical clinical PBS rates: 200–2000 MU/s.
%                      Default: 1000 MU/s (mid-range proxy estimate).
%
% weightUnit (string): describes what unit the optimised spot weights are
%                      expressed in.  matRad weights are dimensionless
%                      fluence equivalents; the exact correspondence to
%                      physical MU depends on the machine calibration.
%                      Until a calibrated conversion factor is available,
%                      weights are treated as 'unknown_or_MU_equivalent'.

machine.irradiation.muRate     = 1000;                    % MU/s
machine.irradiation.weightUnit = 'unknown_or_MU_equivalent';

% -------------------------------------------------------------------------
% Gantry motion — trapezoidal rotation profile
% -------------------------------------------------------------------------
% In proton arc delivery the gantry rotates continuously (or step-and-shoot)
% between consecutive control points. The time to traverse an angular
% separation theta (deg) is estimated using a trapezoidal velocity profile:
%
%   If theta <= vmax^2 / amax  (triangle profile, no constant-speed phase):
%     t = 2 * sqrt(theta / amax)
%   Else (trapezoid profile with cruise phase):
%     t = theta / vmax + vmax / amax
%
% maxSpeed (deg/s)           : peak gantry angular velocity.
%                              IEC 60601 and vendor data suggest 4–6 deg/s.
%                              Default: 6 deg/s (upper proxy estimate).
% maxAccel (deg/s^2)         : peak angular acceleration.
%                              Typical values: 0.5–2 deg/s^2.
%                              Default: 1 deg/s^2.
% toleranceWindow (deg)      : angular tolerance window within which the
%                              gantry is considered to have reached its
%                              target position (used in constraint
%                              validation and feasibility checks).
%                              Default: 0.5 deg.

machine.gantry.maxSpeed        = 6.0;  % deg/s    — peak angular velocity
machine.gantry.maxAccel        = 1.0;  % deg/s^2  — peak angular acceleration
machine.gantry.toleranceWindow = 0.5;  % deg      — positional tolerance

% -------------------------------------------------------------------------
% Buffer / overhead times
% -------------------------------------------------------------------------
% Miscellaneous per-event overhead times that account for verification,
% interlock checks, beam current stabilisation, and other ancillary delays
% not captured by the energy/spot/gantry sub-models.
%
% perLayer         (s) : fixed overhead applied once per delivered energy
%                        layer (e.g. dosimetric verification, layer setup).
% perControlPoint  (s) : fixed overhead applied once per control point
%                        (e.g. gantry-angle confirmation, CP transition
%                        bookkeeping).
%
% These are additive buffer terms; set to 0 to suppress.

machine.buffer.perLayer        = 0.05;  % s — overhead per energy layer
machine.buffer.perControlPoint = 0.02;  % s — overhead per control point

% -------------------------------------------------------------------------
% Validation metadata
% -------------------------------------------------------------------------
% Provenance and validation status of the machine timing model.
%
% status (string) : lifecycle label.
%   'proxy_unvalidated'  — parameters are order-of-magnitude estimates;
%                          BDT output is indicative only.
%   'literature_sourced' — parameters sourced from peer-reviewed data but
%                          not measured on a specific machine.
%   'vendor_spec'        — parameters provided by machine vendor.
%   'validated'          — parameters measured from actual delivery logs.
%
% source (string) : human-readable description of where the values come
%                   from.  Update this field when replacing proxy values
%                   with measured data.

machine.validation.status = 'proxy_unvalidated';
machine.validation.source = 'default_placeholder';

end
