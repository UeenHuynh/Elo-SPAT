function [phiES, gradES, timeValue, gradTime, details] = ...
    evaluate_delivery_surrogates(w, context)
% evaluate_delivery_surrogates  Smooth optimization surrogates for ES/BDT.
%
% The returned timeValue is deliberately not the final discrete BDT. It is
% a differentiable activation surrogate correlated with energy switches,
% active layers/spots, irradiation, gantry traversal and fixed buffers.

w = reshape(double(w), [], 1);
if numel(w) ~= context.nWeights || any(w < 0) || any(~isfinite(w))
    error('evaluate_delivery_surrogates:invalidWeights', ...
        'w must be a finite non-negative vector of the expected length.');
end
epsA = context.epsActivity;
[zGroup, dzGroup] = group_activity(w, context.groups, epsA);
[zCp, dzCp] = group_activity(w, context.cpGroups, epsA);

phiES = 0;
timeEnergy = 0;
dPhiDz = zeros(context.nGroups, 1);
dTimeDz = zeros(context.nGroups, 1);
for e = 1:numel(context.edgeFrom)
    i = context.edgeFrom(e);
    j = context.edgeTo(e);
    pair = zGroup(i) * zGroup(j);
    cSeq = context.edgeSeqCost(e);
    cTime = context.edgeTimeCost(e);
    phiES = phiES + cSeq * pair;
    timeEnergy = timeEnergy + cTime * pair;
    dPhiDz(i) = dPhiDz(i) + cSeq * zGroup(j);
    dPhiDz(j) = dPhiDz(j) + cSeq * zGroup(i);
    dTimeDz(i) = dTimeDz(i) + cTime * zGroup(j);
    dTimeDz(j) = dTimeDz(j) + cTime * zGroup(i);
end
gradES = scatter_activity_gradient(dPhiDz, dzGroup, context.groups, context.nWeights);
gradTime = scatter_activity_gradient(dTimeDz, dzGroup, context.groups, context.nWeights);

epsS = context.epsSpot;
spotActivity = w ./ (w + epsS);
dSpot = epsS ./ (w + epsS) .^ 2;
timeSpot = context.machine.spot.switchTime * sum(spotActivity);
gradTime = gradTime + context.machine.spot.switchTime * dSpot;

timeLayerBuffer = context.machine.buffer.perLayer * sum(zGroup);
gradTime = gradTime + scatter_activity_gradient( ...
    context.machine.buffer.perLayer * ones(context.nGroups, 1), ...
    dzGroup, context.groups, context.nWeights);

timeCpBuffer = context.machine.buffer.perControlPoint * sum(zCp);
gradTime = gradTime + scatter_activity_gradient( ...
    context.machine.buffer.perControlPoint * ones(context.nControlPoints, 1), ...
    dzCp, context.cpGroups, context.nWeights);

timeGantry = 0;
dGantryDz = zeros(context.nControlPoints, 1);
for k = 1:numel(context.gantryEdgeCost)
    c = context.gantryEdgeCost(k);
    timeGantry = timeGantry + c * zCp(k) * zCp(k + 1);
    dGantryDz(k) = dGantryDz(k) + c * zCp(k + 1);
    dGantryDz(k + 1) = dGantryDz(k + 1) + c * zCp(k);
end
gradTime = gradTime + scatter_activity_gradient( ...
    dGantryDz, dzCp, context.cpGroups, context.nWeights);

timeIrradiation = sum(w) / context.machine.irradiation.muRate;
gradTime = gradTime + ones(size(w)) / context.machine.irradiation.muRate;
timeValue = timeEnergy + timeSpot + timeLayerBuffer + timeCpBuffer + ...
    timeGantry + timeIrradiation;
details = struct('energy', timeEnergy, 'spot', timeSpot, ...
    'layerBuffer', timeLayerBuffer, 'controlPointBuffer', timeCpBuffer, ...
    'gantry', timeGantry, 'irradiation', timeIrradiation, ...
    'groupActivity', zGroup, 'spotActivity', spotActivity);
end

function [z, dz] = group_activity(w, groups, epsA)
z = zeros(numel(groups), 1);
dz = cell(numel(groups), 1);
for g = 1:numel(groups)
    ix = groups(g).ix;
    a = sqrt(sum(w(ix) .^ 2) + epsA ^ 2);
    z(g) = (a - epsA) / (a + epsA);
    dz{g} = (2 * epsA / (a + epsA) ^ 2) * (w(ix) / a);
end
end

function grad = scatter_activity_gradient(dFdz, dz, groups, n)
grad = zeros(n, 1);
for g = 1:numel(groups)
    grad(groups(g).ix) = grad(groups(g).ix) + dFdz(g) * dz{g};
end
end
