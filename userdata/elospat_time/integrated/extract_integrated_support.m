function support = extract_integrated_support(w, referenceW, context, cfg)
% extract_integrated_support  Convert optimizer weights to fixed final support.

w = reshape(double(w), [], 1);
referenceW = reshape(double(referenceW), [], 1);
if numel(w) ~= context.nWeights || numel(referenceW) ~= context.nWeights
    error('extract_integrated_support:dimensionMismatch', ...
        'w and referenceW must match the mapping dimension.');
end
groupNorm = zeros(context.nGroups, 1);
layerMask = false(context.nGroups, 1);
for g = 1:context.nGroups
    groupNorm(g) = norm(w(context.groups(g).ix), 2);
    layerMask(g) = groupNorm(g) > cfg.layerThreshold;
end

% Enforce control-point coverage at support level. If a CP vanished, retain
% the strongest layer according to the optimized weights, with Dense Arc as
% a deterministic tie/fallback reference.
cpValues = reshape(double([context.groups.cpIndex]), [], 1);
for b = 1:context.nControlPoints
    cp = context.cpGroups(b).cpIndex;
    candidates = find(cpValues == cp);
    if ~any(layerMask(candidates))
        scores = groupNorm(candidates);
        if all(scores == 0)
            for j = 1:numel(candidates)
                g = candidates(j);
                scores(j) = norm(referenceW(context.groups(g).ix), 2);
            end
        end
        [~, j] = max(scores);
        layerMask(candidates(j)) = true;
    end
end

supportMask = false(context.nWeights, 1);
threshold = max(cfg.spotThreshold, cfg.mmuThreshold);
for g = find(layerMask)'
    ix = context.groups(g).ix;
    local = w(ix) >= threshold;
    if ~any(local)
        [~, j] = max(w(ix));
        if w(ix(j)) <= 0
            [~, j] = max(referenceW(ix));
        end
        local(j) = true;
    end
    supportMask(ix(local)) = true;
end

wSupport = zeros(size(w));
wSupport(supportMask) = w(supportMask);
if cfg.mmuThreshold > 0
    wSupport(supportMask) = max(wSupport(supportMask), cfg.mmuThreshold);
end
support = struct('mask', supportMask, 'layerMask', layerMask, ...
    'w', wSupport, 'groupNorm', groupNorm, ...
    'nSpots', nnz(supportMask), 'nLayers', nnz(layerMask));
end
