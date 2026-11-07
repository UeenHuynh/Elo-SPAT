function tests=test_working_set_pricing, tests=functiontests(localfunctions); end
function testSpotReducedCosts(t)
g=[-3;-1;2;-4]; active=[false;true;false;false];
groups=struct('ix',{[1;2],[3;4]});
p=price_working_set_candidates(g,active,groups,struct('absoluteTolerance',0,'relativeTolerance',0));
verifyEqual(t,p.spotViolation,[3;0;0;4]);
verifyEqual(t,p.spotEligible,[true;false;false;true]);
end
function testLayerSteepestDirection(t)
g=[-3;4;-4]; active=false(3,1); groups=struct('ix',{[1;2;3]});
p=price_working_set_candidates(g,active,groups,struct('absoluteTolerance',0,'relativeTolerance',0));
verifyEqual(t,p.layerViolation,5,'AbsTol',1e-12);
end
function testActiveLayerNotPriced(t)
g=[-3;-4]; active=[true;false]; groups=struct('ix',{[1;2]});
p=price_working_set_candidates(g,active,groups);
verifyFalse(t,p.layerEligible); verifyEqual(t,p.layerViolation,0);
end
function testScaleAwareTolerance(t)
g=[-1e-8;1e3]; groups=struct('ix',{1,2});
p=price_working_set_candidates(g,false(2,1),groups, ...
    struct('absoluteTolerance',1e-7,'relativeTolerance',1e-6));
verifyFalse(t,p.spotEligible(1));
end
