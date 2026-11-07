function test_collective_exchange_core()
rng(31); n=24; deaths=[2;7;13;19]; w=0.5+rand(n,1);
h0=0.2*randn(n,1); h0=max(h0,-0.4*w); h0(deaths)=-w(deaths);
h0(20:24)=h0(20:24)+sum(w(deaths))/5;
R=randn(5,n); C=R-(R*h0)*(h0'/(h0'*h0));
[h,info]=solve_collective_exchange_qp(C,w,deaths,1e-10);
assert(info.success); assert(info.kernelResidual<=1e-8);
assert(info.deathResidual<=1e-8); assert(min(w+h)>=-1e-8);
assert(all(abs(w(deaths)+h(deaths))<=1e-8));
disp('COLLECTIVE_EXCHANGE_CORE_TESTS_PASS');
end
