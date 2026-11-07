function report = validate_matrad_objective_gradient(doseModel, w, sampleIndices, step)
% validate_matrad_objective_gradient  Directional finite-difference audit.
%
% Use this on a reduced dij before any full arc optimization. It checks the
% complete native matRad objective plus configured constraint merit.

narginchk(2, 4);
w = reshape(double(w), [], 1);
if nargin < 3 || isempty(sampleIndices)
    sampleIndices = unique(round(linspace(1, numel(w), min(20, numel(w)))))';
end
if nargin < 4 || isempty(step)
    step = 1e-6 * max(1, norm(w, inf));
end
[value, gradient, diagnostics] = evaluate_matrad_dose_model(doseModel, w, true);
sampleIndices = reshape(double(sampleIndices), [], 1);
finiteDifference = zeros(size(sampleIndices));
analytic = gradient(sampleIndices);
for k = 1:numel(sampleIndices)
    i = sampleIndices(k);
    if w(i) > step
        wp = w; wm = w;
        wp(i) = wp(i) + step;
        wm(i) = wm(i) - step;
        fp = evaluate_matrad_dose_model(doseModel, wp, false);
        fm = evaluate_matrad_dose_model(doseModel, wm, false);
        finiteDifference(k) = (fp - fm) / (2 * step);
    else
        wp = w;
        wp(i) = wp(i) + step;
        fp = evaluate_matrad_dose_model(doseModel, wp, false);
        finiteDifference(k) = (fp - value) / step;
    end
end
absoluteError = abs(analytic - finiteDifference);
relativeError = absoluteError ./ max(1, max(abs(analytic), abs(finiteDifference)));
report = struct('value', value, 'indices', sampleIndices, ...
    'analytic', analytic, 'finiteDifference', finiteDifference, ...
    'absoluteError', absoluteError, 'relativeError', relativeError, ...
    'maximumRelativeError', max(relativeError), ...
    'doseDiagnostics', diagnostics, 'step', step);
end
