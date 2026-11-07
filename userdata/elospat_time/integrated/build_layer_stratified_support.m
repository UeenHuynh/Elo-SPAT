function [supportMask, diagnostics] = build_layer_stratified_support( ...
    referenceW, map, maxSpotsPerControlPoint)
% build_layer_stratified_support  Nested Dense-reference support by CP/layer.
%
% The strongest reference-weight spot from every energy layer is retained
% first. Remaining capacity at each control point is filled in descending
% reference-weight order. For caps larger than the number of layers, supports
% are nested as the cap increases.

narginchk(3, 3);
referenceW = reshape(double(referenceW), [], 1);
if ~isstruct(map) || numel(map) ~= numel(referenceW)
    error('build_layer_stratified_support:dimensionMismatch', ...
        'map and referenceW must have the same number of entries.');
end
if ~isscalar(maxSpotsPerControlPoint) || ...
        ~isfinite(maxSpotsPerControlPoint) || ...
        maxSpotsPerControlPoint < 1 || ...
        floor(maxSpotsPerControlPoint) ~= maxSpotsPerControlPoint
    error('build_layer_stratified_support:invalidCap', ...
        'maxSpotsPerControlPoint must be a positive integer.');
end

cpIndex = reshape(double([map.cpIndex]), [], 1);
energyIndex = reshape(double([map.energyIndex]), [], 1);
cpOrder = unique(cpIndex, 'stable');
supportMask = false(numel(referenceW), 1);
diagnostics = repmat(struct('cpIndex', 0, 'availableSpots', 0, ...
    'availableLayers', 0, 'selectedLayers', 0, 'selectedSpots', 0, ...
    'selectedReferenceWeightFraction', NaN), numel(cpOrder), 1);

for k = 1:numel(cpOrder)
    candidates = find(cpIndex == cpOrder(k));
    layers = unique(energyIndex(candidates), 'stable');
    layerScore = zeros(numel(layers), 1);
    for j = 1:numel(layers)
        layerIx = candidates(energyIndex(candidates) == layers(j));
        layerScore(j) = sum(max(referenceW(layerIx), 0));
    end
    [~, layerOrder] = sort(layerScore, 'descend');
    nLayerKeep = min([numel(layers), maxSpotsPerControlPoint]);
    selectedLayers = layers(layerOrder(1:nLayerKeep));

    keep = zeros(0, 1);
    pool = zeros(0, 1);
    for j = 1:numel(selectedLayers)
        layerIx = candidates(energyIndex(candidates) == selectedLayers(j));
        [~, localOrder] = sort(referenceW(layerIx), 'descend');
        keep(end + 1, 1) = layerIx(localOrder(1)); %#ok<AGROW>
        pool = [pool; layerIx(localOrder(2:end))]; %#ok<AGROW>
    end
    remaining = min(maxSpotsPerControlPoint - numel(keep), numel(pool));
    if remaining > 0
        [~, poolOrder] = sort(referenceW(pool), 'descend');
        keep = [keep; pool(poolOrder(1:remaining))]; %#ok<AGROW>
    end
    supportMask(keep) = true;

    totalWeight = sum(max(referenceW(candidates), 0));
    diagnostics(k).cpIndex = cpOrder(k);
    diagnostics(k).availableSpots = numel(candidates);
    diagnostics(k).availableLayers = numel(layers);
    diagnostics(k).selectedLayers = numel(unique(energyIndex(keep)));
    diagnostics(k).selectedSpots = numel(keep);
    diagnostics(k).selectedReferenceWeightFraction = ...
        sum(max(referenceW(keep), 0)) / max(totalWeight, eps);
end
end
