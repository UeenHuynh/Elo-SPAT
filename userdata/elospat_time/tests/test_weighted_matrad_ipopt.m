% test_weighted_matrad_ipopt  Exact linear-term and KKT smoke test.
rootDir=fullfile(fileparts(mfilename('fullpath')),'..');
addpath(fullfile(rootDir,'integrated'));
b=[1.0;0.25;2.0;0.1]; q=[0.2;0.5;0.1;0.2];
model=struct('nWeights',4,'hasConstraints',false, ...
    'objectiveFunction',@(w)0.5*sum((w-b).^2), ...
    'objectiveGradient',@(w)w-b);
expected=max(0,b-q); w0=ones(4,1)*0.4;
settings=struct('ipopt',struct('print_level',0,'max_iter',100));
r=solve_weighted_matrad_ipopt(model,w0,q,settings);
fprintf('status=%d %s iter=%g accepted=%d\n',r.status,r.statusMessage, ...
    r.iterations,r.innerAccepted);
fprintf('expected:'); fprintf(' %.12g',expected); fprintf('\n');
fprintf('actual:  '); fprintf(' %.12g',r.w); fprintf('\n');
fprintf('KKT primal=%.6e dual=%.6e compl=%.6e window=%.6e\n', ...
    r.postKKT.primalInfeasibility,r.postKKT.dualInfeasibility, ...
    r.postKKT.complementarity,r.stability.relativeRange);
assert(norm(r.w-expected,Inf)<=1e-5,'Weighted IPOPT solution mismatch.');
assert(r.postKKT.primalInfeasibility<=1e-6);
assert(r.postKKT.dualInfeasibility<=1e-3);
assert(r.postKKT.complementarity<=1e-3);
assert(r.callbackCounts.objective>0 && r.callbackCounts.gradient>0);
fprintf('[OK] Weighted matRad/IPOPT linear-term and KKT smoke test passed.\n');
