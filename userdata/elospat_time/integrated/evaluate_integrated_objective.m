function [value, parts] = evaluate_integrated_objective(w, doseModel, context, cfg)
% evaluate_integrated_objective  Evaluate the actual nonconvex objective.

[smoothValue, ~, parts] = evaluate_integrated_smooth( ...
    w, doseModel, context, cfg, false);
layerPenalty = 0;
for g = 1:numel(context.groups)
    layerPenalty = layerPenalty + sqrt(norm(w(context.groups(g).ix), 2));
end
spotPenalty = sum(w ./ (w + cfg.epsS));
value = smoothValue + cfg.lambdaEL * layerPenalty + cfg.lambdaS * spotPenalty;
parts.layerL21half = layerPenalty;
parts.spotSmoothL0 = spotPenalty;
parts.total = value;
end
