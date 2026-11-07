function tests = test_pdasc_l0_ls
tests = functiontests(localfunctions);
end

function testOrthogonalContinuation(testCase)
U = speye(4); p = [4; 3; 2; 1];
s = struct('lambdaMinRatio', 1e-4, 'maxContinuationLevels', 25, ...
    'maxInnerIterations', 10, 'verbose', false);
r = run_pdasc_l0_ls(U, p, s);
verifyEqual(testCase, r.lambdaMax, 8, 'AbsTol', 1e-12);
verifyEqual(testCase, r.history(1).projectedActiveSpots, 0);
verifyTrue(testCase, all(diff([r.history.projectedActiveSpots]) >= 0));
verifyEqual(testCase, r.w, p, 'AbsTol', 1e-10);
verifyTrue(testCase, all(strcmp({r.history.innerStop}, 'active_set_stable')));
end

function testColumnNormalizationIsEquivalent(testCase)
U = sparse([2 0; 0 0.5; 1 1]); p = [2; 1; 2];
base = struct('lambdaMinRatio', 1e-8, 'maxContinuationLevels', 31, ...
    'maxInnerIterations', 20, 'verbose', false);
a = run_pdasc_l0_ls(U, p, base);
base.normalizeColumns = false;
b = run_pdasc_l0_ls(U, p, base);
verifyEqual(testCase, U*a.w, U*b.w, 'AbsTol', 1e-8);
verifyEqual(testCase, a.support, b.support);
end

function testStrictTieRemainsInactive(testCase)
U = speye(2); p = [2; 1];
s = struct('lambdaMinRatio', 1, 'maxContinuationLevels', 2, ...
    'maxInnerIterations', 5, 'verbose', false);
r = run_pdasc_l0_ls(U, p, s);
verifyEqual(testCase, nnz(r.support), 0);
end

function testZeroColumnNeverActivates(testCase)
U = sparse([1 0; 1 0]); p = [1;1];
s = struct('lambdaMinRatio', 1e-6, 'maxContinuationLevels', 15, ...
    'verbose', false);
r = run_pdasc_l0_ls(U, p, s);
verifyFalse(testCase, r.support(2));
end
