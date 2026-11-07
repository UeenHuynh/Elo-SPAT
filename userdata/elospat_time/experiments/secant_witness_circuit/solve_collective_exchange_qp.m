function [h,info] = solve_collective_exchange_qp(C,w,deathLocalIx,rankTolerance)
%SOLVE_COLLECTIVE_EXCHANGE_QP Exact multi-coordinate boundary event.
% min 0.5||W h||^2 s.t. C h=0, h_I=-w_I, h>=-w.
narginchk(4,4);
w=reshape(double(w),[],1); C=double(C); n=numel(w);
deathLocalIx=unique(reshape(double(deathLocalIx),[],1),'stable');
assert(size(C,2)==n && ~isempty(deathLocalIx));
assert(all(deathLocalIx>=1 & deathLocalIx<=n & w(deathLocalIx)>0));

t=tic;
if isempty(C)
    rowBasis=zeros(0,n); s=zeros(0,1); rankC=0;
else
    [~,S,V]=svd(C,'econ'); s=diag(S);
    if isempty(s), scale=1; else, scale=max(1,s(1)); end
    rankC=sum(s>rankTolerance*scale); rowBasis=V(:,1:rankC)';
end
info.svdSeconds=toc(t); info.rank=rankC; info.nullity=n-rankC;
deathEq=sparse(1:numel(deathLocalIx),deathLocalIx,1,numel(deathLocalIx),n);
Aeq=[rowBasis;deathEq]; beq=[zeros(rankC,1);-w(deathLocalIx)];
positive=w(w>0); floorScale=max(median(positive),eps);
relativeScale=max(w,0.05*floorScale);
H=spdiags(1./relativeScale.^2,0,n,n);
opts=optimoptions('quadprog','Display','off','Algorithm','interior-point-convex', ...
    'ConstraintTolerance',1e-10,'OptimalityTolerance',1e-10,'MaxIterations',600);
t=tic;
[h,~,exitflag,output]=quadprog(H,zeros(n,1),[],[],Aeq,beq,-w,[],[],opts);
info.qpSeconds=toc(t); info.exitflag=exitflag; info.output=output;
info.success=~isempty(h) && exitflag>0 && all(isfinite(h));
if ~info.success
    h=zeros(n,1); info.kernelResidual=Inf; info.deathResidual=Inf; return;
end
h=reshape(double(h),[],1);
info.kernelResidual=norm(C*h)/max(1,norm(C,'fro')*norm(h));
info.deathResidual=norm(h(deathLocalIx)+w(deathLocalIx),Inf)/max(1,norm(w(deathLocalIx),Inf));
info.minimumEndpoint=min(w+h);
end
