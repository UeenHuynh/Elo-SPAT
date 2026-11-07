function [rhoNew, uRescaled, action] = balance_admm_rho( ...
    rho, u, primalNorm, dualNorm, settings)
%BALANCE_ADMM_RHO Bounded residual balancing with scaled-dual correction.
%
% The unscaled multiplier y=rho*u is preserved whenever rho changes.

narginchk(5, 5);
rhoNew = rho;
action = 'unchanged';
if primalNorm > settings.mu * dualNorm
    rhoNew = min(settings.maximum, rho * settings.increaseFactor);
    action = 'increased';
elseif dualNorm > settings.mu * primalNorm
    rhoNew = max(settings.minimum, rho / settings.decreaseFactor);
    action = 'decreased';
end
if rhoNew == rho
    action = 'unchanged';
    uRescaled = u;
else
    uRescaled = (rho / rhoNew) * u;
end
end
