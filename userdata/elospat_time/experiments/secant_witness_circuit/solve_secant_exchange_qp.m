function [h,info] = solve_secant_exchange_qp(C,w,anchor,rankTolerance)
%SOLVE_SECANT_EXCHANGE_QP Minimum-relative-change boundary exchange.
% Solves min 0.5*||W h||^2 s.t. C*h=0, h(anchor)=-w(anchor), h>=-w.

narginchk(4,4);
w=reshape(double(w),[],1); C=double(C); n=numel(w);
assert(size(C,2)==n && anchor>=1 && anchor<=n && w(anchor)>0);

tSvd=tic;
if isempty(C)
    rowBasis=zeros(0,n); singularValues=zeros(0,1); rankC=0;
else
    [~,S,V]=svd(C,'econ');
    singularValues=diag(S);
    if isempty(singularValues), scale=1; else, scale=max(1,singularValues(1)); end
    rankC=sum(singularValues > rankTolerance*scale);
    rowBasis=V(:,1:rankC)';
end
info.svdSeconds=toc(tSvd);
info.rank=rankC; info.nullity=n-rankC;
info.singularValues=singularValues;

Aeq=[rowBasis; sparse(1,anchor,1,1,n)];
beq=[zeros(rankC,1);-w(anchor)];
positive=w(w>0); floorScale=max(median(positive),eps);
relativeScale=max(w,0.05*floorScale);
H=spdiags(1./(relativeScale.^2),0,n,n);
opts=optimoptions('quadprog','Display','off','Algorithm','interior-point-convex', ...
    'ConstraintTolerance',1e-10,'OptimalityTolerance',1e-10, ...
    'MaxIterations',400);
tQp=tic;
[h,~,exitflag,output]=quadprog(H,zeros(n,1),[],[],Aeq,beq,-w,[],[],opts);
info.qpSeconds=toc(tQp); info.exitflag=exitflag; info.output=output;
info.success=~isempty(h) && exitflag>0 && all(isfinite(h));
if ~info.success
    h=zeros(n,1); info.kernelResidual=Inf; info.anchorResidual=Inf;
    info.minimumEndpoint=NaN; return;
end
h=reshape(double(h),[],1);
info.kernelResidual=norm(C*h)/max(1,norm(C,'fro')*norm(h));
info.anchorResidual=abs(h(anchor)+w(anchor))/max(1,w(anchor));
info.minimumEndpoint=min(w+h);
end
