function test_numerical_stabilization_rule()
% Unit checks for separation of stabilization, stationarity and plan gates.

root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'integrated'));
rule = default_stabilization_rule();
n = 120;
h = struct('iteration',(1:n)', ...
    'objective',1000 - 0.001*(1:n)', ...
    'relativeObjective',1e-6*ones(n,1), ...
    'stepSize',ones(n,1), ...
    'normalizedStationarity',1e-3*ones(n,1), ...
    'layerSupportChanged',false(n,1), ...
    'spotSupportChanged',false(n,1));
bdt = struct('iteration',[20;70;120], 'value',[5000;5001;5002]);
r = evaluate_numerical_stabilization(h, bdt, rule);
assert(r.numerically_stabilized);
assert(~r.stationary_converged);

h.spotSupportChanged(end) = true;
r = evaluate_numerical_stabilization(h, bdt, rule);
assert(~r.numerically_stabilized);
assert(strcmp(r.reason, 'spot_support_not_stabilized'));
h.spotSupportChanged(end) = false;

h.layerSupportChanged(end) = true;
r = evaluate_numerical_stabilization(h, bdt, rule);
assert(~r.numerically_stabilized);
assert(strcmp(r.reason, 'layer_support_not_stabilized'));

p = evaluate_plan_acceptance(true, 0, true, true, rule);
assert(p.plan_accepted);
p = evaluate_plan_acceptance(false, 0, true, true, rule);
assert(~p.plan_accepted);
fprintf('test_numerical_stabilization_rule: PASS\n');
end
