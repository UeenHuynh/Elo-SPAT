function result = solve_weighted_matrad_ipopt(doseModel, w0, q, settings)
% solve_weighted_matrad_ipopt  Native IPOPT solve of F_matRad(w)+q'*w.
%
% This project-local entry point preserves matRad objective/constraint
% callbacks while adding a linear weight-space term. It intentionally uses
% primal-only warm starts; no cst semantics are modified.

narginchk(3,4);
if nargin<4 || isempty(settings), settings=struct(); end
w0=reshape(double(w0),[],1); q=reshape(double(q),[],1);
if numel(w0)~=doseModel.nWeights || numel(q)~=numel(w0) || ...
        any(~isfinite(w0)) || any(~isfinite(q)) || any(q<0)
    error('solve_weighted_matrad_ipopt:invalidInput', ...
        'w0 and q must be finite vectors matching doseModel.nWeights; q>=0.');
end

[lb,ub,cl,cu]=problem_bounds(doseModel,w0);
defaults=default_settings(); settings=merge_settings(defaults,settings);
baseOptimizer=matRad_OptimizerIPOPT();
ipoptOptions=baseOptimizer.options;
names=fieldnames(settings.ipopt);
for i=1:numel(names), ipoptOptions.(names{i})=settings.ipopt.(names{i}); end

options=struct('lb',lb,'ub',ub,'cl',cl,'cu',cu,'ipopt',ipoptOptions);
counters=struct('objective',0,'gradient',0,'constraints',0,'jacobian',0);
maxRecords=ipoptOptions.max_iter+2;
iterHistory=initialize_iter_history(maxRecords); nRecords=0;
funcs.objective=@objective_callback;
funcs.gradient=@gradient_callback;
if ~isempty(cl)
    funcs.constraints=@constraint_callback;
    funcs.jacobian=@jacobian_callback;
    funcs.jacobianstructure=@jacobian_structure_callback;
end
funcs.iterfunc=@iteration_callback;

timer=tic;
[wRaw,info]=ipopt(w0,funcs,options);
wallSeconds=toc(timer);
iterHistory=trim_iter_history(iterHistory,nRecords);
wRaw=reshape(double(wRaw),[],1);
w=max(0,wRaw);
projectionNorm=norm(w-wRaw);
projectionRelativeNorm=projectionNorm/max(1,norm(wRaw));

[postKKT,postValues]=compute_post_kkt( ...
    doseModel,wRaw,q,lb,ub,cl,cu,info);
stability=objective_stability(iterHistory.objective,settings.stabilityWindow);
status=double(info.status);
desired=status==0 && ...
    postKKT.primalInfeasibility<=settings.desiredPrimalTolerance && ...
    postKKT.dualInfeasibility<=settings.desiredDualTolerance && ...
    postKKT.complementarity<=settings.desiredComplementarityTolerance;
acceptable=any(status==[0 1]) && ...
    postKKT.primalInfeasibility<=settings.acceptablePrimalTolerance && ...
    postKKT.dualInfeasibility<=settings.acceptableDualTolerance && ...
    postKKT.complementarity<=settings.acceptableComplementarityTolerance && ...
    stability.relativeRange<=settings.acceptableObjectiveWindowTolerance && ...
    stability.nValues>=settings.stabilityWindow;

result=struct('w',w,'wRaw',wRaw,'q',q,'info',info, ...
    'status',status,'statusMessage',status_message(status), ...
    'iterations',read_scalar(info,'iter',NaN),'cpuSeconds', ...
    read_scalar(info,'cpu',NaN),'wallSeconds',wallSeconds, ...
    'callbackCounts',counters,'iterationHistory',iterHistory, ...
    'postKKT',postKKT,'postValues',postValues,'stability',stability, ...
    'desiredConverged',desired,'acceptableStabilized',acceptable, ...
    'innerAccepted',desired||acceptable,'settings',settings, ...
    'warmStartMode','primal_only','projectionNorm',projectionNorm, ...
    'projectionRelativeNorm',projectionRelativeNorm);

    function value=objective_callback(x)
        counters.objective=counters.objective+1;
        value=base_objective(doseModel,x)+q'*x;
        value=double(value);
    end
    function gradient=gradient_callback(x)
        counters.gradient=counters.gradient+1;
        gradient=reshape(double(base_gradient(doseModel,x)),[],1)+q;
    end
    function value=constraint_callback(x)
        counters.constraints=counters.constraints+1;
        value=reshape(double(doseModel.optiProb.matRad_constraintFunctions( ...
            x,doseModel.dij,doseModel.cst)),[],1);
    end
    function J=jacobian_callback(x)
        counters.jacobian=counters.jacobian+1;
        J=sparse(doseModel.optiProb.matRad_constraintJacobian( ...
            x,doseModel.dij,doseModel.cst));
    end
    function J=jacobian_structure_callback()
        J=sparse(doseModel.optiProb.matRad_getJacobianStructure( ...
            w0,doseModel.dij,doseModel.cst));
    end
    function keepGoing=iteration_callback(iter,objective,parameter)
        nRecords=nRecords+1;
        if nRecords>numel(iterHistory.iteration)
            iterHistory=grow_iter_history(iterHistory,maxRecords);
        end
        iterHistory.iteration(nRecords)=iter;
        iterHistory.objective(nRecords)=objective;
        iterHistory.infPr(nRecords)=field_or_nan(parameter,'inf_pr');
        iterHistory.infDu(nRecords)=field_or_nan(parameter,'inf_du');
        iterHistory.mu(nRecords)=field_or_nan(parameter,'mu');
        iterHistory.stepNorm(nRecords)=field_or_nan(parameter,'d_norm');
        iterHistory.alphaPr(nRecords)=field_or_nan(parameter,'alpha_pr');
        iterHistory.alphaDu(nRecords)=field_or_nan(parameter,'alpha_du');
        iterHistory.lineSearchTrials(nRecords)=field_or_nan(parameter,'ls_trials');
        keepGoing=true;
    end
end

function settings=default_settings()
settings=struct( ...
    'desiredPrimalTolerance',1e-6, ...
    'desiredDualTolerance',1e-3, ...
    'desiredComplementarityTolerance',1e-3, ...
    'acceptablePrimalTolerance',1e-5, ...
    'acceptableDualTolerance',1e-2, ...
    'acceptableComplementarityTolerance',1e-2, ...
    'acceptableObjectiveWindowTolerance',1e-6, ...
    'stabilityWindow',8);
settings.ipopt=struct('tol',1e-6,'dual_inf_tol',1e-3, ...
    'constr_viol_tol',1e-6,'compl_inf_tol',1e-3, ...
    'acceptable_iter',8,'acceptable_tol',1e-4, ...
    'acceptable_constr_viol_tol',1e-5, ...
    'acceptable_dual_inf_tol',1e-2, ...
    'acceptable_compl_inf_tol',1e-2, ...
    'acceptable_obj_change_tol',1e-6, ...
    'max_iter',300,'max_cpu_time',21600,'print_level',5, ...
    'hessian_approximation','limited-memory', ...
    'limited_memory_max_history',50,'mu_strategy','adaptive');
end

function merged=merge_settings(base,override)
merged=base; names=fieldnames(override);
for i=1:numel(names)
    key=names{i};
    if isfield(merged,key)&&isstruct(merged.(key))&&isstruct(override.(key))
        merged.(key)=merge_settings(merged.(key),override.(key));
    else
        merged.(key)=override.(key);
    end
end
end

function [lb,ub,cl,cu]=problem_bounds(model,w)
if isfield(model,'optiProb')
    lb=reshape(double(model.optiProb.lowerBounds(w)),[],1);
    ub=reshape(double(model.optiProb.upperBounds(w)),[],1);
    [cl,cu]=model.optiProb.matRad_getConstraintBounds(model.cst);
    cl=reshape(double(cl),[],1); cu=reshape(double(cu),[],1);
else
    lb=zeros(size(w)); ub=Inf(size(w)); cl=zeros(0,1); cu=zeros(0,1);
end
end
function value=base_objective(model,w)
if isfield(model,'objectiveFunction')
    value=model.objectiveFunction(w);
else
    value=model.optiProb.matRad_objectiveFunction(w,model.dij,model.cst);
end
end
function gradient=base_gradient(model,w)
if isfield(model,'objectiveGradient')
    gradient=model.objectiveGradient(w);
else
    gradient=model.optiProb.matRad_objectiveGradient(w,model.dij,model.cst);
end
end
function [kkt,values]=compute_post_kkt(model,w,q,lb,ub,cl,cu,info)
f=double(base_objective(model,w)); grad=reshape(double(base_gradient(model,w)),[],1)+q;
if isempty(cl)
    c=zeros(0,1); J=sparse(0,numel(w)); lambda=zeros(0,1);
else
    c=reshape(double(model.optiProb.matRad_constraintFunctions( ...
        w,model.dij,model.cst)),[],1);
    J=sparse(model.optiProb.matRad_constraintJacobian(w,model.dij,model.cst));
    lambda=read_vector(info,'lambda',numel(cl));
end
zl=read_vector(info,'zl',numel(w)); zu=read_vector(info,'zu',numel(w));
stationarity=grad+J'*lambda-zl+zu;
boundViolation=max([0;lb-w;w-ub]);
constraintViolation=max([0;cl-c;c-cu]);
primal=max(boundViolation,constraintViolation);
lowerComp=max(abs((w-lb).*zl));
finiteUpper=isfinite(ub);
if any(finiteUpper), upperComp=max(abs((ub(finiteUpper)-w(finiteUpper)).*zu(finiteUpper))); else, upperComp=0; end
kkt=struct('primalInfeasibility',primal, ...
    'dualInfeasibility',norm(stationarity,Inf), ...
    'complementarity',max(lowerComp,upperComp), ...
    'boundViolation',boundViolation,'constraintViolation',constraintViolation, ...
    'stationarityTwoNorm',norm(stationarity));
values=struct('doseObjective',f,'linearTerm',q'*w, ...
    'innerObjective',f+q'*w,'gradientNorm',norm(grad));
end
function v=read_vector(s,name,n)
if isstruct(s)&&isfield(s,name)&&numel(s.(name))==n
    v=reshape(double(s.(name)),[],1);
else
    v=zeros(n,1);
end
end
function value=read_scalar(s,name,fallback)
if isstruct(s)&&isfield(s,name)&&isnumeric(s.(name))&&isscalar(s.(name))
    value=double(s.(name));
else
    value=fallback;
end
end
function value=field_or_nan(s,name)
if isstruct(s)&&isfield(s,name), value=double(s.(name)); else, value=NaN; end
end
function h=initialize_iter_history(n)
names={'iteration','objective','infPr','infDu','mu','stepNorm', ...
    'alphaPr','alphaDu','lineSearchTrials'};
for i=1:numel(names), h.(names{i})=NaN(n,1); end
end
function h=grow_iter_history(h,n)
names=fieldnames(h);
for i=1:numel(names), h.(names{i})=[h.(names{i});NaN(n,1)]; end
end
function h=trim_iter_history(h,n)
names=fieldnames(h);
for i=1:numel(names), h.(names{i})=h.(names{i})(1:n); end
end
function s=objective_stability(values,window)
values=values(isfinite(values)); n=min(numel(values),window);
if n==0
    s=struct('nValues',0,'relativeRange',Inf,'relativeLastChange',Inf); return;
end
v=values(end-n+1:end);
s=struct('nValues',n,'relativeRange',(max(v)-min(v))/max(1,abs(v(end))), ...
    'relativeLastChange',Inf);
if n>=2, s.relativeLastChange=abs(v(end)-v(end-1))/max(1,abs(v(end))); end
end
function text=status_message(status)
switch status
    case 0, text='Solved';
    case 1, text='Solved To Acceptable Level';
    case 2, text='Infeasible Problem Detected';
    case 3, text='Search Direction Too Small';
    case -1, text='Maximum Iterations Exceeded';
    case -4, text='Maximum CPU Time Exceeded';
    otherwise, text=sprintf('IPOPT status %g',status);
end
end
