function tests=test_working_set_quadratic,tests=functiontests(localfunctions);end
function testReducedCostTermination(t)
% Pure pricing-level surrogate for min 0.5||w-b||^2, w>=0.
b=[3;0;2;-1]; w=zeros(4,1); active=false(4,1);
groups=struct('ix',{[1;2],[3;4]});
for k=1:4
    g=w-b; p=price_working_set_candidates(g,active,groups, ...
        struct('absoluteTolerance',1e-12,'relativeTolerance',0));
    ix=find(p.spotEligible);
    if isempty(ix),break;end
    [~,q]=max(p.spotViolation(ix));active(ix(q))=true;w(active)=max(b(active),0);
end
verifyEqual(t,w,[3;0;2;0],'AbsTol',1e-12);
verifyEqual(t,max([-g(~active);0]),0,'AbsTol',1e-12);
end
