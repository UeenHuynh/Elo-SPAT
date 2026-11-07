function result=run_full_spot_working_set_solver(problem,settings)
%RUN_FULL_SPOT_WORKING_SET_SOLVER Exact restricted solves and full spot pricing.
% Convergence requires both an accepted restricted solve and no inactive
% nonnegative reduced-cost violation. State is checkpointed after every solve.
narginchk(2,2); settings=defaults(settings);
n=double(problem.dij.totalNumOfBixels); w=reshape(double(problem.initialW),[],1);
assert(numel(w)==n); active=w>0;
history=repmat(empty_history(),settings.maxOuterIterations,1);
objectiveEvaluations=0;gradientEvaluations=settings.initialGradientEvaluations;
constraintEvaluations=0;doseMetricEvaluations=0;
timer=tic;stopReason='maximum_outer_iterations';outerUsed=0;previousObjective=Inf;

for k=1:settings.maxOuterIterations
    outerUsed=k;activeBefore=active;
    [dijActive,activeToFull]=reduce_dij_to_support(problem.dij,active);
    modelActive=build_matrad_dose_model(dijActive,problem.cst,problem.pln,0);
    inner=solve_weighted_matrad_ipopt(modelActive,w(activeToFull), ...
        zeros(nnz(active),1),settings.innerSettings);
    objectiveEvaluations=objectiveEvaluations+inner.callbackCounts.objective;
    gradientEvaluations=gradientEvaluations+inner.callbackCounts.gradient;
    constraintEvaluations=constraintEvaluations+inner.callbackCounts.constraints;
    w(:)=0;w(activeToFull)=inner.w;

    [g,pricingMode]=full_pricing_gradient(problem.fullModel,w,inner);
    gradientEvaluations=gradientEvaluations+1;
    pricing=price_working_set_candidates(g,active,problem.groups,settings);
    removable=active&w<=settings.zeroTolerance&g>=-pricing.tolerance;
    if settings.removeZeroNonbasic,active(removable)=false;w(removable)=0;end
    pricing=price_working_set_candidates(g,active,problem.groups,settings);
    eligible=find(pricing.spotEligible);
    [score,order]=sort(pricing.spotViolation(eligible),'descend');
    take=min(settings.spotBatchSize,numel(order));addMask=false(n,1);
    if take>0,addMask(eligible(order(1:take)))=true;cutoff=score(take);else,cutoff=NaN;end

    [objective,~,diagnostics]=evaluate_matrad_dose_model(problem.fullModel,w,false);
    objectiveEvaluations=objectiveEvaluations+1;
    dose=problem.doseFunction(w);doseMetricEvaluations=doseMetricEvaluations+1;
    bdt=problem.bdtFunction(w); meaningful=problem.meaningfulFunction(dose,diagnostics,objective);
    h=empty_history();h.outerIteration=k;h.activeBefore=nnz(activeBefore);
    h.workingSetAfterSolve=nnz(active);h.positiveSpots=nnz(w>settings.zeroTolerance);
    h.activeLayers=count_layers(w,problem.groups,settings.zeroTolerance);
    h.addedSpots=nnz(addMask);h.removableSpots=nnz(removable);
    h.eligibleInactiveSpots=numel(eligible);h.rejectedEligibleSpots=numel(eligible)-take;
    h.objective=objective;h.objectiveChange=objective-previousObjective;
    h.objectiveMonotone=~isfinite(previousObjective)||objective<=previousObjective+ ...
        settings.masterMonotonicityTolerance*max(1,abs(previousObjective));
    h.maximumInactiveSpotViolation=pricing.maximumInactiveSpotViolation;
    h.pricingTolerance=pricing.tolerance;h.admissionCutoff=cutoff;
    h.pricingMode=pricingMode;h.maxConstraintViolation=diagnostics.maxConstraintViolation;
    h.ipoptStatus=inner.statusMessage;h.ipoptIterations=inner.iterations;
    h.restrictedSolveAccepted=inner.innerAccepted;
    h.ipoptPrimal=inner.postKKT.primalInfeasibility;
    h.ipoptDual=inner.postKKT.dualInfeasibility;
    h.ipoptComplementarity=inner.postKKT.complementarity;
    h.objectiveEvaluations=objectiveEvaluations;h.gradientEvaluations=gradientEvaluations;
    h.constraintEvaluations=constraintEvaluations;h.doseMetricEvaluations=doseMetricEvaluations;
    h.targetD95=dose.summary.targetD95Min;h.targetD98=dose.summary.targetD98Min;
    h.oarDmax=dose.summary.oarDmaxMax;h.clinicallyMeaningful=meaningful;
    h.provisionalBDT=bdt.T_total;h.elapsedSeconds=toc(timer);
    h.supportHash=support_hash(active);history(k)=h;previousObjective=objective;
    save_state(k,inner);
    fprintf(['FULL-WS k=%d set=%d positive=%d layers=%d add=%d removable=%d ' ...
        'eligible=%d obj=%.9e rc=%.3e D95=%.4f D98=%.4f BDT=%.2f IPOPT=%s\n'], ...
        k,nnz(active),h.positiveSpots,h.activeLayers,h.addedSpots,h.removableSpots, ...
        h.eligibleInactiveSpots,objective,h.maximumInactiveSpotViolation, ...
        h.targetD95,h.targetD98,h.provisionalBDT,h.ipoptStatus);
    if ~inner.innerAccepted,stopReason='restricted_solve_not_accepted';break;end
    if ~any(addMask),stopReason='restricted_and_pricing_converged';break;end
    active=active|addMask;w(addMask)=settings.newColumnInitialWeight;
end
result=make_state([]);save(fullfile(settings.outputDir,'full_ws_selector_final.mat'),'result','-v7.3');

    function save_state(iteration,innerResult)
        state=make_state(innerResult); %#ok<NASGU>
        save(fullfile(settings.outputDir,sprintf('full_ws_iter_%03d.mat',iteration)), ...
            'state','-v7.3');
        writetable(struct2table(state.history),fullfile(settings.outputDir,'full_ws_trajectory.csv'));
    end
    function state=make_state(innerResult)
        state=struct('w',w,'workingSet',active,'history',history(1:outerUsed), ...
            'outerIterations',outerUsed,'stopReason',stopReason,'wallSeconds',toc(timer), ...
            'objectiveEvaluations',objectiveEvaluations,'gradientEvaluations',gradientEvaluations, ...
            'constraintEvaluations',constraintEvaluations, ...
            'doseMetricEvaluations',doseMetricEvaluations,'lastInner',innerResult, ...
            'settings',settings);
    end
end

function [g,mode]=full_pricing_gradient(model,w,inner)
g=reshape(double(model.optiProb.matRad_objectiveGradient(w,model.dij,model.cst)),[],1);
mode='objective_gradient';
if model.hasConstraints&&isstruct(inner.info)&&isfield(inner.info,'lambda')
    lambda=reshape(double(inner.info.lambda),[],1);
    if numel(lambda)==numel(model.constraintLower)
        J=sparse(model.optiProb.matRad_constraintJacobian(w,model.dij,model.cst));
        g=g+J'*lambda;mode='lagrangian_gradient';
    end
end
end
function n=count_layers(w,groups,tol),n=0;for q=1:numel(groups),n=n+any(w(groups(q).ix)>tol);end,end
function s=defaults(s)
d=struct('spotBatchSize',1000,'maxOuterIterations',20,'absoluteTolerance',1e-7, ...
    'relativeTolerance',1e-6,'zeroTolerance',1e-10,'removeZeroNonbasic',true, ...
    'newColumnInitialWeight',0,'masterMonotonicityTolerance',1e-8, ...
    'initialGradientEvaluations',1,'innerSettings',struct(),'outputDir','');
f=fieldnames(d);for i=1:numel(f),if ~isfield(s,f{i}),s.(f{i})=d.(f{i});end,end
assert(~isempty(s.outputDir));
end
function h=empty_history()
h=struct('outerIteration',NaN,'activeBefore',NaN,'workingSetAfterSolve',NaN, ...
    'positiveSpots',NaN,'activeLayers',NaN,'addedSpots',NaN,'removableSpots',NaN, ...
    'eligibleInactiveSpots',NaN,'rejectedEligibleSpots',NaN,'objective',NaN, ...
    'objectiveChange',NaN,'objectiveMonotone',false,'maximumInactiveSpotViolation',NaN, ...
    'pricingTolerance',NaN,'admissionCutoff',NaN,'pricingMode','', ...
    'maxConstraintViolation',NaN,'ipoptStatus','','ipoptIterations',NaN, ...
    'restrictedSolveAccepted',false,'ipoptPrimal',NaN,'ipoptDual',NaN, ...
    'ipoptComplementarity',NaN,'objectiveEvaluations',NaN,'gradientEvaluations',NaN, ...
    'constraintEvaluations',NaN,'doseMetricEvaluations',NaN,'targetD95',NaN, ...
    'targetD98',NaN,'oarDmax',NaN,'clinicallyMeaningful',false, ...
    'provisionalBDT',NaN,'elapsedSeconds',NaN,'supportHash','');
end
function value=support_hash(mask)
md=java.security.MessageDigest.getInstance('SHA-256');md.update(uint8(mask(:)));
bytes=typecast(md.digest(),'uint8');value=lower(reshape(dec2hex(bytes,2).',1,[]));
end
