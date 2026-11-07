function test_compute_solver_metrics()
% test_compute_solver_metrics  Focused validation for Task 32 helper.

fprintf('\n========================================\n');
fprintf('TEST: compute_solver_metrics (Task 32)\n');
fprintf('========================================\n');

resultGUI = struct();
resultGUI.info = struct('objective', 123.4, 'iter', 17, 'status', 1, 'cpu', 8.5);
resultGUI.usedOptimizer = 'IPOPT';
resultGUI.converged = true;
resultGUI.statusMessage = 'Solved To Acceptable Level';

meta = struct();
meta.runtime = 9.1;
meta.exitStatus = 1;
meta.converged = true;
meta.optimizerName = 'IPOPT';
meta.optimizerResultInfo = struct('iterations', 19);

metrics = compute_solver_metrics(resultGUI, meta);

assert(abs(metrics.finalObjective - 123.4) < 1e-12);
assert(metrics.iterations == 19);
assert(abs(metrics.runtimeSeconds - 9.1) < 1e-12);
assert(metrics.converged == true);
assert(metrics.exitStatus == 1);
assert(strcmp(metrics.optimizerName, 'IPOPT'));
assert(strcmp(metrics.statusMessage, 'Solved To Acceptable Level'));

fallbackMetrics = compute_solver_metrics(struct('info', struct('objective', 50, 'iter', 7, 'status', -1, 'cpu', 3.2)));
assert(fallbackMetrics.converged == false);
assert(fallbackMetrics.exitStatus == -1);
assert(strcmp(fallbackMetrics.statusMessage, 'Maximum Number of Iterations Exceeded'));
assert(fallbackMetrics.acceptableButNotConverged == true);

fprintf('[PASS] Solver metrics extracted from matRad-style metadata.\n');
end
