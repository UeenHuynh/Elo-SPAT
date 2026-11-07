function test_secant_witness_core()
%TEST_SECANT_WITNESS_CORE Deterministic correctness tests for local mechanics.
rng(7); n=12; w=0.2+rand(n,1); anchor=4;
hFeasible=0.1*randn(n,1); hFeasible=max(hFeasible,-0.5*w);
hFeasible(anchor)=-w(anchor); hFeasible(5)=hFeasible(5)+2*w(anchor);
R=randn(3,n); C=R-(R*hFeasible)*(hFeasible'/(hFeasible'*hFeasible));
[h,qp]=solve_secant_exchange_qp(C,w,anchor,1e-10);
assert(qp.success,'Synthetic exchange QP failed.');
assert(norm(C*h)<=1e-8*max(1,norm(C,'fro')*norm(h)));
assert(min(w+h)>=-1e-9);
assert(abs(w(anchor)+h(anchor))<=1e-8);

Q=orth(randn(30,4)); a0=randn(30,1); v=randn(30,1);
delta=0.37;
[Q2,s]=append_secant_witness(Q,a0,v,delta,1e-10);
assert(s.secantError<=1e-11,'Secant identity failed.');
assert(s.independent && size(Q2,2)==size(Q,2)+1);
assert(abs(s.abar'*v-delta)<=1e-10*max(1,abs(delta)));

% A direction in the old kernel is excluded by a non-orthogonal witness.
A=randn(30,n); Cold=Q'*A; Z=null(Cold); assert(~isempty(Z)); oldh=Z(:,1);
trialV=A*oldh; [Q3,s2]=append_secant_witness(Q,a0,trialV,1.2,1e-10);
assert(s2.independent);
assert(norm(Q'*A*oldh)<=1e-9*max(1,norm(A)*norm(oldh)));
assert(norm(Q3'*A*oldh)>1e-8,'New witness did not exclude old direction.');
disp('SECANT_WITNESS_CORE_TESTS_PASS');
end
