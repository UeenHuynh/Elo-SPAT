function tests = test_prox_l0_mmu
tests = functiontests(localfunctions);
end

function testNegativeAndSubMMUInputs(testCase)
a = [-2; 0.2; 0.8; 3];
m = [1; 1; 1; 1];
[z, info] = prox_l0_mmu(a, 0.1, 1, m);
verifyEqual(testCase, z, [0; 0; 1; 3], 'AbsTol', 1e-14);
verifyTrue(testCase, info.mmuFeasible);
end

function testExactTieChoosesZeroAboveMMU(testCase)
lambda0 = 2;
rho = 1;
a = sqrt(2 * lambda0 / rho);
[z, info] = prox_l0_mmu(a, lambda0, rho, 1);
verifyEqual(testCase, z, 0, 'AbsTol', 1e-14);
verifyTrue(testCase, info.tie);
end

function testExactTieChoosesZeroBelowMMU(testCase)
% For a<m, equality solves lambda + rho/2(m-a)^2 = rho/2 a^2,
% hence a = m/2 + lambda/(rho*m).
lambda0 = 0.1;
rho = 2;
m = 1;
a = m/2 + lambda0/(rho*m);
verifyLessThan(testCase, a, m);
[z, info] = prox_l0_mmu(a, lambda0, rho, m);
verifyEqual(testCase, z, 0, 'AbsTol', 1e-14);
verifyTrue(testCase, info.tie);
end

function testHeterogeneousMinimumWeights(testCase)
a = [0.4; 1.4; 2.4; -0.1];
m = [0.5; 1.5; 3.0; 0.01];
lambda0 = [0.01; 0.2; 0.1; 0];
rho = 2;
[z, info] = prox_l0_mmu(a, lambda0, rho, m);

% Verify against exhaustive evaluation of the only two candidate states.
expected = zeros(size(a));
for i = 1:numel(a)
    on = max(a(i), m(i));
    c0 = 0.5 * rho * a(i)^2;
    c1 = lambda0(i) + 0.5 * rho * (on-a(i))^2;
    if c1 < c0
        expected(i) = on;
    end
end
verifyEqual(testCase, z, expected, 'AbsTol', 1e-14);
verifyTrue(testCase, info.mmuFeasible);
end

function testZeroLambdaIsExactProjection(testCase)
a = [-1; 0.2; 0.6; 2];
m = [0.1; 1; 1; 3];
z = prox_l0_mmu(a, 0, 1, m);
verifyEqual(testCase, z, [0; 0; 1; 3], 'AbsTol', 1e-14);
end

function testNoMandatoryLayerRepair(testCase)
% Both spots may be eliminated. This is required by revised Variant H.
z = prox_l0_mmu([0.01; 0.02], 1, 1, [0.1; 0.1]);
verifyEqual(testCase, z, zeros(2, 1));
end
