function tests = test_admm_rho_balance
tests = functiontests(localfunctions);
end

function testIncreasePreservesUnscaledDual(testCase)
s = settings();
rho = 1; u = [2; -3];
[rhoNew, uNew, action] = balance_admm_rho(rho, u, 30, 1, s);
verifyEqual(testCase, rhoNew, 2);
verifyEqual(testCase, action, 'increased');
verifyEqual(testCase, rhoNew*uNew, rho*u, 'AbsTol', 1e-14);
end

function testDecreaseAndLowerBound(testCase)
s = settings();
[rhoNew, uNew, action] = balance_admm_rho(0.125, [2; 1], 0, 10, s);
verifyEqual(testCase, rhoNew, 0.125);
verifyEqual(testCase, uNew, [2; 1]);
verifyEqual(testCase, action, 'unchanged');
end

function testNoChangeInsideBalanceBand(testCase)
s = settings();
[rhoNew, uNew, action] = balance_admm_rho(1, [1; 2], 3, 2, s);
verifyEqual(testCase, rhoNew, 1);
verifyEqual(testCase, uNew, [1; 2]);
verifyEqual(testCase, action, 'unchanged');
end

function testSwitchingBoundaryBothBranches(testCase)
rho = 2;
% First entry has sqrt(2 lambda/rho)>=m; second uses sub-MMU branch.
lambda = [4; 0.01];
m = [1; 1];
b = l0_mmu_switching_boundary(lambda, rho, m);
verifyEqual(testCase, b(1), 2, 'AbsTol', 1e-14);
verifyEqual(testCase, b(2), 0.505, 'AbsTol', 1e-14);
for i = 1:2
    [~, info] = prox_l0_mmu(b(i), lambda(i), rho, m(i));
    verifyEqual(testCase, info.offCost, info.onCost, 'AbsTol', 1e-12);
end
end

function s = settings()
s = struct('mu', 2, 'increaseFactor', 2, 'decreaseFactor', 2, ...
    'minimum', 0.125, 'maximum', 8);
end
