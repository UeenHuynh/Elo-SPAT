function [val, switchUpCount] = energy_sequencing_penalty(w, groups, epsZ)
% energy_sequencing_penalty  Compute the smooth energy switch-up penalty.
%
% INPUTS
%   w      - spot-weight vector (N x 1 or 1 x N)
%   groups - struct array with fields:
%            .ix, .cpIndex, .energy
%   epsZ   - soft activity smoothing constant (must be positive)
%
% OUTPUTS
%   val            - smooth sequencing penalty value
%   switchUpCount  - hard count of active low-to-high transitions

narginchk(3, 3);

if ~isnumeric(w) || isempty(w) || ~isvector(w)
    error('energy_sequencing_penalty:invalidInput', ...
        'w must be a non-empty numeric vector.');
end

if ~isstruct(groups) || isempty(groups)
    error('energy_sequencing_penalty:invalidInput', ...
        'groups must be a non-empty struct array.');
end

requiredFields = {'ix', 'cpIndex', 'energy'};
for i = 1:numel(requiredFields)
    if ~isfield(groups, requiredFields{i})
        error('energy_sequencing_penalty:missingField', ...
            'groups.%s is required.', requiredFields{i});
    end
end

if ~isscalar(epsZ) || ~isnumeric(epsZ) || ~isfinite(epsZ) || epsZ <= 0
    error('energy_sequencing_penalty:invalidInput', ...
        'epsZ must be a positive finite scalar.');
end

w = reshape(double(w), [], 1);
n = numel(w);
nGroups = numel(groups);

activity = zeros(nGroups, 1);
cpIndex = zeros(nGroups, 1);
energy = zeros(nGroups, 1);

for g = 1:nGroups
    ix = reshape(double(groups(g).ix), [], 1);

    if isempty(ix)
        error('energy_sequencing_penalty:emptyGroup', ...
            'groups(%d).ix must not be empty.', g);
    end

    if any(~isfinite(ix)) || any(ix < 1) || any(ix > n) || any(floor(ix) ~= ix)
        error('energy_sequencing_penalty:indexOutOfRange', ...
            'groups(%d).ix contains invalid indices for w of length %d.', g, n);
    end

    if ~isscalar(groups(g).cpIndex) || ~isfinite(groups(g).cpIndex)
        error('energy_sequencing_penalty:invalidInput', ...
            'groups(%d).cpIndex must be a finite scalar.', g);
    end

    if ~isscalar(groups(g).energy) || ~isfinite(groups(g).energy)
        error('energy_sequencing_penalty:invalidInput', ...
            'groups(%d).energy must be a finite scalar.', g);
    end

    activity(g) = norm(w(ix), 2);
    cpIndex(g) = double(groups(g).cpIndex);
    energy(g) = double(groups(g).energy);
end

softActivity = activity ./ (activity + epsZ);
uniqueCp = unique(cpIndex, 'sorted');

val = 0;
switchUpCount = 0;

for k = 1:(numel(uniqueCp) - 1)
    currentMask = cpIndex == uniqueCp(k);
    nextMask = cpIndex == uniqueCp(k + 1);

    currentSoft = softActivity(currentMask);
    nextSoft = softActivity(nextMask);
    currentActivity = activity(currentMask);
    nextActivity = activity(nextMask);
    currentEnergy = energy(currentMask);
    nextEnergy = energy(nextMask);

    for i = 1:numel(currentEnergy)
        deltaEnergy = nextEnergy - currentEnergy(i);
        switchMagnitude = max(deltaEnergy, 0);

        val = val + currentSoft(i) * sum(nextSoft .* switchMagnitude);
        switchUpCount = switchUpCount + sum( ...
            currentActivity(i) > 0 & nextActivity > 0 & deltaEnergy > 0);
    end
end

switchUpCount = double(switchUpCount);

end
