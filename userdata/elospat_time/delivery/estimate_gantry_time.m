function t = estimate_gantry_time(theta, machine)
% estimate_gantry_time  Compute gantry rotation time via trapezoidal profile.
%
% INPUTS
%   theta   - angular difference in degrees
%   machine - machine model with gantry.maxSpeed and gantry.maxAccel
%
% OUTPUT
%   t - estimated gantry rotation time in seconds

narginchk(2, 2);

if ~isscalar(theta) || ~isnumeric(theta) || ~isfinite(theta)
    error('estimate_gantry_time:invalidInput', ...
        'theta must be a finite numeric scalar.');
end

if ~isstruct(machine) || ~isfield(machine, 'gantry') || ~isstruct(machine.gantry) || ...
        ~isfield(machine.gantry, 'maxSpeed') || ~isfield(machine.gantry, 'maxAccel')
    error('estimate_gantry_time:missingField', ...
        'machine.gantry.maxSpeed and machine.gantry.maxAccel are required.');
end

vmax = double(machine.gantry.maxSpeed);
amax = double(machine.gantry.maxAccel);

if ~isscalar(vmax) || ~isnumeric(vmax) || ~isfinite(vmax) || vmax <= 0
    error('estimate_gantry_time:invalidInput', ...
        'machine.gantry.maxSpeed must be a positive finite scalar.');
end

if ~isscalar(amax) || ~isnumeric(amax) || ~isfinite(amax) || amax <= 0
    error('estimate_gantry_time:invalidInput', ...
        'machine.gantry.maxAccel must be a positive finite scalar.');
end

theta = abs(double(theta));

if theta <= (vmax ^ 2) / amax
    t = 2 * sqrt(theta / amax);
else
    t = theta / vmax + vmax / amax;
end

end
