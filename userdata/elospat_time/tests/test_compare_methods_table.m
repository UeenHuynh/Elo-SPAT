function test_compare_methods_table()
% test_compare_methods_table  Focused validation for Task 33 helper.

fprintf('\n========================================\n');
fprintf('TEST: compare_methods_table (Task 33)\n');
fprintf('========================================\n');

methods(1).name = 'Dense Arc';
methods(1).doseMetrics = struct('summary', struct('targetD95Min', 73.1, 'oarDmaxMax', 75.0));
methods(1).deliveryMetrics = struct('nActiveLayers', 50, 'T_total', 120.5);
methods(1).solverMetrics = struct('finalObjective', 3000, 'iterations', 500, 'converged', false);

methods(2).name = 'RMS-ELO-SPAT';
methods(2).doseMetrics = struct('summary', struct('targetD95Min', 72.4, 'oarDmaxMax', 74.1));
methods(2).deliveryMetrics = struct('nActiveLayers', 18, 'T_total', 61.2);
methods(2).solverMetrics = struct('finalObjective', 2800, 'iterations', 185, 'converged', true);

tbl = compare_methods_table(methods);

assert(istable(tbl));
assert(height(tbl) == 2);
assert(any(strcmp(tbl.Properties.VariableNames, 'dose_targetD95Min')));
assert(any(strcmp(tbl.Properties.VariableNames, 'delivery_T_total')));
assert(any(strcmp(tbl.Properties.VariableNames, 'solver_iterations')));
assert(strcmp(tbl.Method(1), "Dense Arc"));
assert(abs(tbl.delivery_T_total(2) - 61.2) < 1e-12);
assert(tbl.solver_converged(1) == 0);
assert(tbl.solver_converged(2) == 1);

fprintf('[PASS] Comparison table assembled across methods.\n');
end
