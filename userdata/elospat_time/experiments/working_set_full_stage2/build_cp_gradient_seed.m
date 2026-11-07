function [mask,diagnostics]=build_cp_gradient_seed(gradient,wReference,map,perCp)
%BUILD_CP_GRADIENT_SEED CP-balanced, clinical-gradient working-set seed.
% All candidates outside this seed remain available to subsequent pricing.
narginchk(4,4);
g=reshape(double(gradient),[],1);w=reshape(double(wReference),[],1);
assert(numel(g)==numel(w)&&numel(map)==numel(w));
cp=reshape(double([map.cpIndex]),[],1); cpOrder=unique(cp,'stable');
mask=false(size(w));
diagnostics=repmat(struct('cpIndex',0,'availableSpots',0, ...
    'negativeGradientSpots',0,'selectedSpots',0, ...
    'selectedNegativeGradientSpots',0,'coverageFillSpots',0, ...
    'minimumSelectedViolation',NaN),numel(cpOrder),1);
for k=1:numel(cpOrder)
    ix=find(cp==cpOrder(k)); violation=max(-g(ix),0);
    negative=violation>0; take=min(perCp,numel(ix));
    % First-order improvement is primary. Stage-1 weight only breaks ties and
    % fills a CP quota when fewer than perCp columns have negative g_j at zero.
    ranking=[-violation,-max(w(ix),0),double(ix)];
    [~,order]=sortrows(ranking,[1 2 3]); selected=ix(order(1:take));
    mask(selected)=true;
    diagnostics(k).cpIndex=cpOrder(k);
    diagnostics(k).availableSpots=numel(ix);
    diagnostics(k).negativeGradientSpots=nnz(negative);
    diagnostics(k).selectedSpots=take;
    diagnostics(k).selectedNegativeGradientSpots=nnz(g(selected)<0);
    diagnostics(k).coverageFillSpots=nnz(g(selected)>=0);
    diagnostics(k).minimumSelectedViolation=min(max(-g(selected),0));
end
end
