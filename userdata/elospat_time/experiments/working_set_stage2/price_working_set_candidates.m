function pricing = price_working_set_candidates(gradient, activeMask, groups, settings)
%PRICE_WORKING_SET_CANDIDATES First-order spot and layer reduced costs.
narginchk(3,4);
if nargin<4, settings=struct(); end
settings=with_defaults(settings);
g=reshape(double(gradient),[],1); activeMask=reshape(logical(activeMask),[],1);
if numel(g)~=numel(activeMask) || any(~isfinite(g))
    error('price_working_set_candidates:invalidInput','Gradient/mask mismatch.');
end
scale=max(1,norm(g,Inf)); tolerance=settings.absoluteTolerance+ ...
    settings.relativeTolerance*scale;
inactive=~activeMask;
spotViolation=zeros(size(g));
spotViolation(inactive)=max(-g(inactive),0);
spotEligible=inactive & spotViolation>tolerance;

nGroups=numel(groups); layerViolation=zeros(nGroups,1);
layerEligible=false(nGroups,1); layerInactive=false(nGroups,1);
for q=1:nGroups
    ix=reshape(double(groups(q).ix),[],1);
    if any(ix<1|ix>numel(g)), error('price_working_set_candidates:groups','Bad group index.'); end
    layerInactive(q)=~any(activeMask(ix));
    if layerInactive(q)
        layerViolation(q)=norm(max(-g(ix),0),2);
        layerEligible(q)=layerViolation(q)>tolerance;
    end
end
pricing=struct('gradient',g,'scale',scale,'tolerance',tolerance, ...
    'spotViolation',spotViolation,'spotEligible',spotEligible, ...
    'layerViolation',layerViolation,'layerEligible',layerEligible, ...
    'layerInactive',layerInactive,'maximumInactiveSpotViolation', ...
    max([0;spotViolation(inactive)]),'maximumInactiveLayerViolation', ...
    max([0;layerViolation(layerInactive)]));
end
function s=with_defaults(s)
if ~isfield(s,'absoluteTolerance'),s.absoluteTolerance=1e-7;end
if ~isfield(s,'relativeTolerance'),s.relativeTolerance=1e-6;end
end
