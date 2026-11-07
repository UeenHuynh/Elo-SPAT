function result=run_matrad_working_set(problem,settings)
%RUN_MATRAD_WORKING_SET Exact restricted matRad solves plus full pricing.
%
% Required problem fields:
%   dij, cst, pln, fullModel, groups
% Optional:
%   bdtFunction(w), initialW

narginchk(2,2); settings=defaults(settings); validate_problem(problem);
n=double(problem.dij.totalNumOfBixels);
active=false(n,1); w=zeros(n,1);
if isfield(problem,'initialW') && ~isempty(problem.initialW)
    w=reshape(double(problem.initialW),[],1); active=w>0;
end
history=repmat(empty_history(),settings.maxOuterIterations,1);
supportHashes=cell(settings.maxOuterIterations,1);
objectiveEvaluations=0; gradientEvaluations=0; constraintEvaluations=0;
timer=tic; stopReason='maximum_outer_iterations'; outerUsed=0;
previousObjective=Inf;

% Empty-support initialization is priced directly from the full objective.
[g,pricingMode]=full_pricing_gradient(problem.fullModel,w,[]);
gradientEvaluations=gradientEvaluations+1;
pricing=price_working_set_candidates(g,active,problem.groups,settings);
[addMask,selection]=select_batch(pricing,g,active,problem.groups,settings);
if ~any(addMask)
    stopReason='no_improving_column_at_zero';
    result=finish(); return;
end
active=active|addMask;

for k=1:settings.maxOuterIterations
    outerUsed=k; activeBefore=active; wBefore=w;
    [dijActive,activeToFull]=reduce_dij_to_support(problem.dij,active);
    modelActive=build_matrad_dose_model(dijActive,problem.cst,problem.pln,0);
    warm=w(activeToFull);
    innerSettings=settings.innerSettings;
    inner=solve_weighted_matrad_ipopt(modelActive,warm,zeros(nnz(active),1),innerSettings);
    objectiveEvaluations=objectiveEvaluations+inner.callbackCounts.objective;
    gradientEvaluations=gradientEvaluations+inner.callbackCounts.gradient;
    constraintEvaluations=constraintEvaluations+inner.callbackCounts.constraints;
    w(:)=0; w(activeToFull)=inner.w;

    [g,pricingMode]=full_pricing_gradient(problem.fullModel,w,inner);
    gradientEvaluations=gradientEvaluations+1;
    pricing=price_working_set_candidates(g,active,problem.groups,settings);

    removed=false(n,1);
    if settings.removeZeroNonbasic
        removable=active & w<=settings.zeroTolerance & ...
            g>=-pricing.tolerance;
        active(removable)=false; w(removable)=0; removed=removable;
        pricing=price_working_set_candidates(g,active,problem.groups,settings);
    end
    [addMask,selection]=select_batch(pricing,g,active,problem.groups,settings);

    [objective,diagnostics]=full_objective(problem.fullModel,w);
    objectiveEvaluations=objectiveEvaluations+1;
    bdt=NaN;
    if isfield(problem,'bdtFunction') && ~isempty(problem.bdtFunction)
        bdt=problem.bdtFunction(w);
    end
    h=empty_history(); h.outerIteration=k; h.activeBefore=nnz(activeBefore);
    h.activeAfterSolve=nnz(w>settings.zeroTolerance);
    h.workingSetAfterSolve=nnz(active); h.addedSpots=nnz(addMask);
    h.removedZeroSpots=nnz(removed); h.activeLayers=count_layers(w,problem.groups,settings.zeroTolerance);
    h.objective=objective; h.doseObjective=diagnostics.doseObjective;
    h.objectiveChange=objective-previousObjective;
    h.masterObjectiveMonotone=~isfinite(previousObjective)|| ...
        objective<=previousObjective+settings.masterMonotonicityTolerance*max(1,abs(previousObjective));
    h.maxConstraintViolation=diagnostics.maxConstraintViolation;
    h.maximumInactiveSpotViolation=pricing.maximumInactiveSpotViolation;
    h.maximumInactiveLayerViolation=pricing.maximumInactiveLayerViolation;
    h.pricingTolerance=pricing.tolerance; h.pricingMode=pricingMode;
    h.batchMode=settings.pricingMode; h.selectedBatchScore=selection.score;
    h.ipoptStatus=inner.statusMessage; h.ipoptIterations=inner.iterations;
    h.ipoptAccepted=inner.innerAccepted; h.ipoptPrimal=inner.postKKT.primalInfeasibility;
    h.ipoptDual=inner.postKKT.dualInfeasibility;
    h.ipoptComplementarity=inner.postKKT.complementarity;
    h.objectiveEvaluations=objectiveEvaluations;
    h.gradientEvaluations=gradientEvaluations;
    h.constraintEvaluations=constraintEvaluations;
    h.provisionalBDT=bdt; h.elapsedSeconds=toc(timer);
    h.supportHash=support_hash(active); history(k)=h; supportHashes{k}=h.supportHash;
    previousObjective=objective;
    if settings.verbose
        fprintf(['WS %s k=%d set=%d positive=%d add=%d remove=%d ' ...
            'obj=%.9e rcSpot=%.3e rcLayer=%.3e IPOPT=%s\n'], ...
            settings.pricingMode,k,nnz(active),nnz(w>settings.zeroTolerance), ...
            nnz(addMask),nnz(removed),objective,pricing.maximumInactiveSpotViolation, ...
            pricing.maximumInactiveLayerViolation,inner.statusMessage);
    end
    if ~inner.innerAccepted
        stopReason='restricted_solve_not_accepted'; break;
    end
    if ~any(addMask)
        stopReason='no_inactive_reduced_cost_violation'; break;
    end
    active=active|addMask;
    w(active & w==0)=settings.newColumnInitialWeight;
    if any(strcmp(supportHashes(1:k-1),support_hash(active)))
        stopReason='working_set_cycle'; break;
    end
    if norm(w-wBefore)==0 && isequal(active,activeBefore)
        stopReason='numerical_stagnation'; break;
    end
end
result=finish();

    function out=finish()
        historyTrim=history(1:outerUsed);
        out=struct('w',w,'workingSet',active,'history',historyTrim, ...
            'outerIterations',outerUsed,'stopReason',stopReason, ...
            'wallSeconds',toc(timer),'pricingMode',settings.pricingMode, ...
            'settings',settings,'objectiveEvaluations',objectiveEvaluations, ...
            'gradientEvaluations',gradientEvaluations, ...
            'constraintEvaluations',constraintEvaluations);
    end
end

function [g,mode]=full_pricing_gradient(model,w,inner)
g=reshape(double(model.optiProb.matRad_objectiveGradient(w,model.dij,model.cst)),[],1);
mode='objective_gradient';
if model.hasConstraints && ~isempty(inner) && isfield(inner,'info') && ...
        isstruct(inner.info) && isfield(inner.info,'lambda')
    lambda=reshape(double(inner.info.lambda),[],1);
    if numel(lambda)==numel(model.constraintLower)
        J=sparse(model.optiProb.matRad_constraintJacobian(w,model.dij,model.cst));
        g=g+J'*lambda; mode='lagrangian_gradient';
    end
end
end
function [value,d]=full_objective(model,w)
[value,~,d]=evaluate_matrad_dose_model(model,w,false);
end
function [addMask,selection]=select_batch(p,g,active,groups,s)
addMask=false(size(active)); selection=struct('score',NaN);
switch s.pricingMode
    case 'spot'
        ix=find(p.spotEligible); [scores,order]=sort(p.spotViolation(ix),'descend');
        take=min(s.spotBatchSize,numel(order));
        if take>0, addMask(ix(order(1:take)))=true; selection.score=scores(take); end
    case 'layer'
        % Once a layer is active, its still-inactive spots remain ordinary
        % nonnegative columns and must also satisfy the spot KKT condition.
        inActiveLayer=false(size(active));
        for q0=1:numel(groups)
            ix0=reshape(double(groups(q0).ix),[],1);
            if any(active(ix0)),inActiveLayer(ix0)=true;end
        end
        spotIx=find(p.spotEligible & inActiveLayer);
        [~,spotOrder]=sort(p.spotViolation(spotIx),'descend');
        spotTake=min(s.spotBatchSize,numel(spotOrder));
        if spotTake>0,addMask(spotIx(spotOrder(1:spotTake)))=true;end
        q=find(p.layerEligible); [scores,order]=sort(p.layerViolation(q),'descend');
        take=min(s.layerBatchSize,numel(order));
        for j=1:take
            ix=reshape(double(groups(q(order(j))).ix),[],1);
            addMask(ix(~active(ix)&g(ix)<-p.tolerance))=true;
        end
        if take>0, selection.score=scores(take); end
    otherwise
        error('run_matrad_working_set:pricingMode','Use spot or layer pricing.');
end
end
function n=count_layers(w,groups,tol)
n=0; for q=1:numel(groups),n=n+any(w(groups(q).ix)>tol);end
end
function s=defaults(s)
d=struct('pricingMode','spot','spotBatchSize',50,'layerBatchSize',10, ...
    'absoluteTolerance',1e-7,'relativeTolerance',1e-6, ...
    'masterMonotonicityTolerance',1e-8, ...
    'maxOuterIterations',30,'zeroTolerance',1e-10, ...
    'removeZeroNonbasic',true,'newColumnInitialWeight',0, ...
    'verbose',true,'innerSettings',struct());
n=fieldnames(d); for i=1:numel(n),if ~isfield(s,n{i}),s.(n{i})=d.(n{i});end,end
s.pricingMode=lower(char(s.pricingMode));
end
function validate_problem(p)
for n={'dij','cst','pln','fullModel','groups'}
    if ~isfield(p,n{1}),error('run_matrad_working_set:problem','Missing %s.',n{1});end
end
end
function h=empty_history()
h=struct('outerIteration',NaN,'activeBefore',NaN,'activeAfterSolve',NaN, ...
    'workingSetAfterSolve',NaN,'addedSpots',NaN,'removedZeroSpots',NaN, ...
    'activeLayers',NaN,'objective',NaN,'objectiveChange',NaN, ...
    'masterObjectiveMonotone',false,'doseObjective',NaN, ...
    'maxConstraintViolation',NaN,'maximumInactiveSpotViolation',NaN, ...
    'maximumInactiveLayerViolation',NaN,'pricingTolerance',NaN, ...
    'pricingMode','','batchMode','','selectedBatchScore',NaN, ...
    'ipoptStatus','','ipoptIterations',NaN,'ipoptAccepted',false, ...
    'ipoptPrimal',NaN,'ipoptDual',NaN,'ipoptComplementarity',NaN, ...
    'objectiveEvaluations',NaN,'gradientEvaluations',NaN, ...
    'constraintEvaluations',NaN,'provisionalBDT',NaN, ...
    'elapsedSeconds',NaN,'supportHash','');
end
function value=support_hash(mask)
md=java.security.MessageDigest.getInstance('SHA-256');md.update(uint8(mask(:)));
bytes=typecast(md.digest(),'uint8');value=lower(reshape(dec2hex(bytes,2).',1,[]));
end
