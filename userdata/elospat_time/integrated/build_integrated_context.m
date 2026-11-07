function context = build_integrated_context(groups, cpGroups, map, machine, cfg)
% build_integrated_context  Precompute sparse group and delivery metadata.

narginchk(5, 5);
if ~isstruct(groups) || isempty(groups) || ~isstruct(cpGroups) || isempty(cpGroups)
    error('build_integrated_context:invalidGroups', ...
        'groups and cpGroups must be non-empty struct arrays.');
end
n = numel(map);
groupOfSpot = zeros(n, 1);
for g = 1:numel(groups)
    ix = reshape(double(groups(g).ix), [], 1);
    if any(ix < 1 | ix > n | ix ~= floor(ix)) || any(groupOfSpot(ix) ~= 0)
        error('build_integrated_context:invalidPartition', ...
            'Energy-layer groups must form a non-overlapping spot partition.');
    end
    groupOfSpot(ix) = g;
end
if any(groupOfSpot == 0)
    error('build_integrated_context:invalidPartition', ...
        'Every spot must belong to exactly one energy-layer group.');
end

cpValues = reshape(double([groups.cpIndex]), [], 1);
cpOrder = unique(cpValues, 'stable');
edgeFrom = zeros(0, 1);
edgeTo = zeros(0, 1);
edgeSeqCost = zeros(0, 1);
edgeTimeCost = zeros(0, 1);
for k = 1:(numel(cpOrder) - 1)
    from = find(cpValues == cpOrder(k));
    to = find(cpValues == cpOrder(k + 1));
    for ii = 1:numel(from)
        for jj = 1:numel(to)
            dE = double(groups(to(jj)).energy) - double(groups(from(ii)).energy);
            edgeFrom(end + 1, 1) = from(ii); %#ok<AGROW>
            edgeTo(end + 1, 1) = to(jj); %#ok<AGROW>
            edgeSeqCost(end + 1, 1) = max(dE, 0); %#ok<AGROW>
            if dE > 0
                t = machine.energy.upBase + machine.energy.perMeVUp * dE;
            elseif dE < 0
                t = machine.energy.downBase + machine.energy.perMeVDown * abs(dE);
            else
                t = 0;
            end
            edgeTimeCost(end + 1, 1) = t; %#ok<AGROW>
        end
    end
end

gantryEdgeCost = zeros(max(numel(cpGroups) - 1, 0), 1);
for k = 1:numel(gantryEdgeCost)
    theta = abs(double(cpGroups(k + 1).gantryAngle) - ...
        double(cpGroups(k).gantryAngle));
    theta = min(theta, 360 - theta);
    threshold = machine.gantry.maxSpeed ^ 2 / machine.gantry.maxAccel;
    if theta <= threshold
        gantryEdgeCost(k) = 2 * sqrt(theta / machine.gantry.maxAccel);
    else
        gantryEdgeCost(k) = theta / machine.gantry.maxSpeed + ...
            machine.gantry.maxSpeed / machine.gantry.maxAccel;
    end
end

context = struct();
context.groups = groups;
context.cpGroups = cpGroups;
context.groupOfSpot = groupOfSpot;
context.nWeights = n;
context.nGroups = numel(groups);
context.nControlPoints = numel(cpGroups);
context.edgeFrom = edgeFrom;
context.edgeTo = edgeTo;
context.edgeSeqCost = edgeSeqCost;
context.edgeTimeCost = edgeTimeCost;
context.gantryEdgeCost = gantryEdgeCost;
context.machine = machine;
context.epsActivity = cfg.epsZ;
context.epsSpot = cfg.epsS;
end
