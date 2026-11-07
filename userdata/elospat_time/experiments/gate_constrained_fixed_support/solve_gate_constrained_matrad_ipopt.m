function result=solve_gate_constrained_matrad_ipopt(model,w0,gateA,gateUpper,settings)
% Exact native matRad objective with additional linear dose-gate rows.
narginchk(4,5);if nargin<5,settings=struct();end
w0=reshape(double(w0),[],1);gateA=sparse(double(gateA));gateUpper=reshape(double(gateUpper),[],1);
if isscalar(gateUpper),gateUpper=repmat(gateUpper,size(gateA,1),1);end
assert(numel(w0)==model.nWeights&&size(gateA,2)==model.nWeights&&size(gateA,1)==numel(gateUpper));
lb=reshape(double(model.optiProb.lowerBounds(w0)),[],1);ub=reshape(double(model.optiProb.upperBounds(w0)),[],1);
[nativeCl,nativeCu]=model.optiProb.matRad_getConstraintBounds(model.cst);
nativeCl=reshape(double(nativeCl),[],1);nativeCu=reshape(double(nativeCu),[],1);
cl=[nativeCl;-inf(size(gateUpper))];cu=[nativeCu;gateUpper];
base=matRad_OptimizerIPOPT();opts=base.options;
opts.constr_viol_tol=1e-7;opts.acceptable_constr_viol_tol=1e-7;opts.print_level=5;
if isfield(settings,'ipopt'),n=fieldnames(settings.ipopt);for i=1:numel(n),opts.(n{i})=settings.ipopt.(n{i});end,end
options=struct('lb',lb,'ub',ub,'cl',cl,'cu',cu,'ipopt',opts);
counters=struct('objective',0,'gradient',0,'constraints',0,'jacobian',0);
funcs.objective=@objective_callback;funcs.gradient=@gradient_callback;
funcs.constraints=@constraint_callback;funcs.jacobian=@jacobian_callback;
funcs.jacobianstructure=@jacobian_structure_callback;
t=tic;[w,info]=ipopt(w0,funcs,options);wall=toc(t);w=reshape(double(w),[],1);
[f,g,diag]=evaluate_matrad_dose_model(model,w,true);
c=constraint_callback_no_count(w);J=jacobian_callback_no_count(w);
lambda=read_vector(info,'lambda',numel(cl));zl=read_vector(info,'zl',numel(w));zu=read_vector(info,'zu',numel(w));
stationarity=g+J'*lambda-zl+zu;
primal=max([0;lb-w;w-ub;cl-c;c-cu]);
finiteUpper=isfinite(ub);lc=max(abs((w-lb).*zl));uc=0;if any(finiteUpper),uc=max(abs((ub(finiteUpper)-w(finiteUpper)).*zu(finiteUpper)));end
status=double(info.status);accepted=any(status==[0 1])&&primal<=1e-6;
result=struct('w',w,'info',info,'status',status,'statusMessage',status_message(status), ...
 'solverAccepted',accepted,'iterations',read_scalar(info,'iter',NaN),'wallSeconds',wall, ...
 'callbackCounts',counters,'objective',f,'objectiveGradientNorm',norm(g), ...
 'maxConstraintViolation',diag.maxConstraintViolation,'gateValues',gateA*w, ...
 'gateUpper',gateUpper,'gateViolation',max([0;gateA*w-gateUpper]), ...
 'primalInfeasibility',primal,'dualInfeasibility',norm(stationarity,Inf), ...
 'complementarity',max(lc,uc),'settings',settings,'ipoptOptions',opts);

 function v=objective_callback(x),counters.objective=counters.objective+1;v=double(model.optiProb.matRad_objectiveFunction(x,model.dij,model.cst));end
 function g0=gradient_callback(x),counters.gradient=counters.gradient+1;g0=reshape(double(model.optiProb.matRad_objectiveGradient(x,model.dij,model.cst)),[],1);end
 function v=constraint_callback(x),counters.constraints=counters.constraints+1;v=constraint_callback_no_count(x);end
 function v=constraint_callback_no_count(x)
  if isempty(nativeCl),native=zeros(0,1);else,native=reshape(double(model.optiProb.matRad_constraintFunctions(x,model.dij,model.cst)),[],1);end
  v=[native;gateA*x];
 end
 function J0=jacobian_callback(x),counters.jacobian=counters.jacobian+1;J0=jacobian_callback_no_count(x);end
 function J0=jacobian_callback_no_count(x)
  if isempty(nativeCl),native=sparse(0,numel(x));else,native=sparse(model.optiProb.matRad_constraintJacobian(x,model.dij,model.cst));end
  J0=[native;gateA];
 end
 function J0=jacobian_structure_callback()
  if isempty(nativeCl),native=sparse(0,numel(w0));else,native=sparse(model.optiProb.matRad_getJacobianStructure(w0,model.dij,model.cst));end
  J0=[native;spones(gateA)];
 end
end
function v=read_vector(s,n,k),if isstruct(s)&&isfield(s,n)&&numel(s.(n))==k,v=reshape(double(s.(n)),[],1);else,v=zeros(k,1);end,end
function v=read_scalar(s,n,d),if isstruct(s)&&isfield(s,n)&&isscalar(s.(n)),v=double(s.(n));else,v=d;end,end
function s=status_message(x),switch x,case 0,s='Solved';case 1,s='Solved To Acceptable Level';case 2,s='Infeasible Problem Detected';case 3,s='Search Direction Too Small';case -1,s='Maximum Iterations Exceeded';case -4,s='Maximum CPU Time Exceeded';otherwise,s=sprintf('IPOPT status %g',x);end,end
