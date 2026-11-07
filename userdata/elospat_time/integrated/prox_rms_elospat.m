function w = prox_rms_elospat(z, step, groups, lambdaEL, spotWeights)
% prox_rms_elospat  Nonnegative spot-weighted/group-half proximal map.
%
% Solves each disjoint group subproblem for
%   0.5||x-z||^2 + step*sum_i spotWeights_i*x_i
%                    + step*lambdaEL*sqrt(||x||_2),  x >= 0.

z = reshape(double(z), [], 1);
spotWeights = reshape(double(spotWeights), [], 1);
if numel(z) ~= numel(spotWeights) || step <= 0 || lambdaEL < 0 || ...
        any(spotWeights < 0)
    error('prox_rms_elospat:invalidInput', 'Invalid proximal-map inputs.');
end
v = max(0, z - step * spotWeights);
w = zeros(size(v));
tau = step * lambdaEL;
for g = 1:numel(groups)
    ix = groups(g).ix;
    vg = v(ix);
    radius = norm(vg, 2);
    if radius == 0
        continue;
    end
    if tau == 0
        w(ix) = vg;
        continue;
    end
    threshold = 1.5 * tau ^ (2 / 3);
    if radius <= threshold
        continue;
    end
    arg = -3 * sqrt(3) * tau / (4 * radius ^ (3 / 2));
    arg = min(1, max(-1, arg));
    u = 2 * sqrt(radius / 3) * cos(acos(arg) / 3);
    r = u ^ 2;
    w(ix) = (r / radius) * vg;
end
end
