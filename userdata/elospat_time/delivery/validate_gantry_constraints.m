function [feasible, violations] = validate_gantry_constraints(seq, machine)
% validate_gantry_constraints  Check basic gantry consistency constraints.
%
% INPUTS
%   seq     - delivery sequence struct array
%   machine - machine timing model with gantry fields
%
% OUTPUTS
%   feasible   - true if no violations are found
%   violations - struct array describing each detected issue

narginchk(2, 2);

if ~isstruct(seq)
    error('validate_gantry_constraints:invalidInput', ...
        'seq must be a struct array.');
end

if isempty(seq)
    feasible = true;
    violations = repmat(struct( ...
        'transitionIndex', 0, ...
        'type', '', ...
        'message', '', ...
        'theta', 0, ...
        'transitionTime', 0), 0, 1);
    return;
end

requiredFields = {'cpIndex', 'gantryAngle'};
for i = 1:numel(requiredFields)
    if ~isfield(seq, requiredFields{i})
        error('validate_gantry_constraints:missingField', ...
            'seq.%s is required.', requiredFields{i});
    end
end

if ~isstruct(machine) || ~isfield(machine, 'gantry') || ~isstruct(machine.gantry) || ...
        ~isfield(machine.gantry, 'maxSpeed') || ~isfield(machine.gantry, 'toleranceWindow')
    error('validate_gantry_constraints:missingField', ...
        'machine.gantry.maxSpeed and machine.gantry.toleranceWindow are required.');
end

violations = repmat(struct( ...
    'transitionIndex', 0, ...
    'type', '', ...
    'message', '', ...
    'theta', 0, ...
    'transitionTime', 0), 0, 1);

for i = 2:numel(seq)
    theta = abs(double(seq(i).gantryAngle) - double(seq(i - 1).gantryAngle));
    transitionTime = estimate_gantry_time(theta, machine);

    if double(seq(i).cpIndex) == double(seq(i - 1).cpIndex) && ...
            theta > double(machine.gantry.toleranceWindow)
        violations(end + 1, 1) = make_violation( ... %#ok<AGROW>
            i - 1, 'same_cp_angle_mismatch', ...
            sprintf('Layers within CP %d differ by %.3f deg, exceeding tolerance %.3f deg.', ...
            double(seq(i).cpIndex), theta, double(machine.gantry.toleranceWindow)), ...
            theta, transitionTime);
    end

    if theta > 0
        avgSpeed = theta / max(transitionTime, eps);
        if avgSpeed > double(machine.gantry.maxSpeed) + double(machine.gantry.toleranceWindow)
            violations(end + 1, 1) = make_violation( ... %#ok<AGROW>
                i - 1, 'speed_exceeded', ...
                sprintf('Average gantry speed %.3f deg/s exceeds limit %.3f deg/s.', ...
                avgSpeed, double(machine.gantry.maxSpeed)), ...
                theta, transitionTime);
        end
    end
end

feasible = isempty(violations);

end

function violation = make_violation(transitionIndex, type, message, theta, transitionTime)
violation = struct( ...
    'transitionIndex', transitionIndex, ...
    'type', type, ...
    'message', message, ...
    'theta', theta, ...
    'transitionTime', transitionTime);
end
