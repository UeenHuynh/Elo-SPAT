function out=solve_sparse_nnls_fista(A,b,x0,settings)
%SOLVE_SPARSE_NNLS_FISTA Matrix-free-style projected FISTA for sparse A.
narginchk(4,4);x=max(0,reshape(double(x0),[],1));y=x;t=1;
maxIter=value(settings,'maxIterations',200);tol=value(settings,'stationarityTolerance',1e-5);
window=value(settings,'stabilityWindow',8);objTol=value(settings,'objectiveTolerance',1e-7);
L=estimate_lipschitz(A,value(settings,'powerIterations',12));
history=repmat(struct('iteration',0,'objective',NaN,'relativeObjective',NaN, ...
    'normalizedStationarity',NaN,'stepSize',NaN),maxIter,1);
previous=objective(A,b,x);stable=0;timer=tic;reason='maximum_iterations';
for k=1:maxIter
    gy=A'*(A*y-b);accepted=false;
    for bt=1:20
        candidate=max(0,y-gy/L);f=objective(A,b,candidate);
        delta=candidate-y;major=objective(A,b,y)+gy'*delta+0.5*L*(delta'*delta);
        if f<=major+1e-10*max(1,abs(major)),accepted=true;break;end
        L=2*L;
    end
    if ~accepted,reason='line_search_failed';break;end
    mapping=L*(y-candidate);stationarity=norm(mapping)/max(1,norm(candidate));
    rel=abs(f-previous)/max(1,abs(previous));
    history(k)=struct('iteration',k,'objective',f,'relativeObjective',rel, ...
        'normalizedStationarity',stationarity,'stepSize',1/L);
    if rel<=objTol,stable=stable+1;else,stable=0;end
    if stationarity<=tol&&stable>=window
        x=candidate;reason='stationary_and_objective_stable';break;
    end
    tNew=0.5*(1+sqrt(1+4*t*t));yNew=candidate+((t-1)/tNew)*(candidate-x);
    if (candidate-x)'*(yNew-candidate)>0,yNew=candidate;tNew=1;end
    x=candidate;y=yNew;t=tNew;previous=f;
end
used=find([history.iteration]>0,1,'last');if isempty(used),used=0;end
out=struct('x',x,'history',history(1:used),'iterations',used,'reason',reason, ...
    'objective',objective(A,b,x),'wallSeconds',toc(timer),'lipschitz',L, ...
    'objectiveEvaluations',used+1,'gradientEvaluations',used);
end
function f=objective(A,b,x),r=A*x-b;f=0.5*(r'*r);end
function L=estimate_lipschitz(A,n)
v=ones(size(A,2),1);v=v/max(norm(v),eps);
for k=1:n,z=A'*(A*v);nz=norm(z);if nz==0,L=1;return;end;v=z/nz;end
L=max(v'*(A'*(A*v)),eps)*1.05;
end
function x=value(s,n,d),if isfield(s,n),x=s.(n);else,x=d;end,end
