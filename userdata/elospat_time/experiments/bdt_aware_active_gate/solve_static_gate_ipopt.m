function result=solve_static_gate_ipopt(model,w0,A,clGate,cuGate,settings)
% Native matRad objective with a fixed sparse linear clinical-gate matrix.
w0=reshape(double(w0),[],1);A=sparse(double(A));clGate=reshape(double(clGate),[],1);cuGate=reshape(double(cuGate),[],1);
assert(size(A,2)==numel(w0)&&size(A,1)==numel(clGate)&&numel(clGate)==numel(cuGate));if nargin<6,settings=struct();end
lb=reshape(double(model.optiProb.lowerBounds(w0)),[],1);ub=reshape(double(model.optiProb.upperBounds(w0)),[],1);[ncl,ncu]=model.optiProb.matRad_getConstraintBounds(model.cst);ncl=reshape(double(ncl),[],1);ncu=reshape(double(ncu),[],1);cl=[ncl;clGate];cu=[ncu;cuGate];
base=matRad_OptimizerIPOPT();opts=base.options;opts.constr_viol_tol=1e-7;opts.acceptable_constr_viol_tol=1e-7;opts.print_level=3;if isfield(settings,'ipopt'),f=fieldnames(settings.ipopt);for i=1:numel(f),opts.(f{i})=settings.ipopt.(f{i});end,end
options=struct('lb',lb,'ub',ub,'cl',cl,'cu',cu,'ipopt',opts);count=struct('objective',0,'gradient',0,'constraints',0,'jacobian',0);funcs.objective=@obj;funcs.gradient=@grad;funcs.constraints=@con;funcs.jacobian=@jac;funcs.jacobianstructure=@jacstr;
t=tic;[w,info]=ipopt(w0,funcs,options);wall=toc(t);w=reshape(double(w),[],1);[fv,gv,diag]=evaluate_matrad_dose_model(model,w,true);[cv,J]=all_constraints(w);lambda=readv(info,'lambda',numel(cl));zl=readv(info,'zl',numel(w));zu=readv(info,'zu',numel(w));stationarity=gv+J'*lambda-zl+zu;primal=max([0;lb-w;w-ub;cl-cv;cv-cu]);finite=isfinite(ub);comp=max(abs((w-lb).*zl));if any(finite),comp=max(comp,max(abs((ub(finite)-w(finite)).*zu(finite))));end
status=double(info.status);result=struct('w',w,'info',info,'status',status,'statusMessage',status_text(status),'iterations',readn(info,'iter',NaN),'wallSeconds',wall,'callbackCounts',count,'objective',fv,'maxConstraintViolation',diag.maxConstraintViolation,'primalInfeasibility',primal,'dualInfeasibility',norm(stationarity,Inf),'complementarity',comp,'numericallyStationary',status==0&&primal<=1e-6&&norm(stationarity,Inf)<=1e-3&&comp<=1e-3,'solverUsable',any(status==[0 1 -1 -4])&&primal<=1e-6);
 function v=obj(x),count.objective=count.objective+1;v=double(model.optiProb.matRad_objectiveFunction(x,model.dij,model.cst));end
 function v=grad(x),count.gradient=count.gradient+1;v=reshape(double(model.optiProb.matRad_objectiveGradient(x,model.dij,model.cst)),[],1);end
 function v=con(x),count.constraints=count.constraints+1;[v,~]=all_constraints(x);end
 function v=jac(x),count.jacobian=count.jacobian+1;[~,v]=all_constraints(x);end %#ok<INUSD>
 function v=jacstr(),native=sparse(0,numel(w0));if ~isempty(ncl),native=sparse(model.optiProb.matRad_getJacobianStructure(w0,model.dij,model.cst));end;v=[native;spones(A)];end
 function [v,J0]=all_constraints(x),if isempty(ncl),nv=zeros(0,1);nJ=sparse(0,numel(x));else,nv=reshape(double(model.optiProb.matRad_constraintFunctions(x,model.dij,model.cst)),[],1);nJ=sparse(model.optiProb.matRad_constraintJacobian(x,model.dij,model.cst));end;v=[nv;A*x];J0=[nJ;A];end
end
function v=readv(s,n,k),if isfield(s,n)&&numel(s.(n))==k,v=reshape(double(s.(n)),[],1);else,v=zeros(k,1);end,end
function v=readn(s,n,d),if isfield(s,n)&&isscalar(s.(n)),v=double(s.(n));else,v=d;end,end
function s=status_text(x),switch x,case 0,s='Solved';case 1,s='Solved To Acceptable Level';case 2,s='Infeasible Problem Detected';case -1,s='Maximum Iterations Exceeded';case -4,s='Maximum CPU Time Exceeded';otherwise,s=sprintf('IPOPT status %g',x);end,end
