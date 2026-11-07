function result=run_stage1_dose_coreset(problem,settings)
%RUN_STAGE1_DOSE_CORESET Greedy active-column weighted dose reconstruction.
n=size(problem.D,2);active=reshape(logical(problem.initialMask),[],1);
w=zeros(n,1);w(active)=problem.referenceW(active);timer=tic;
history=repmat(empty_history(),settings.maxOuterIterations,1);
objEvals=0;gradEvals=0;stopReason='maximum_outer_iterations';used=0;
for k=1:settings.maxOuterIterations
    used=k;before=active;As=problem.weightedD(:,active);
    nnls=solve_sparse_nnls_fista(As,problem.weightedReference,w(active),settings.nnls);
    objEvals=objEvals+nnls.objectiveEvaluations;gradEvals=gradEvals+nnls.gradientEvaluations;
    w(:)=0;w(active)=nnls.x;
    removable=active&w<=settings.zeroTolerance;active(removable)=false;w(removable)=0;
    residual=problem.weightedD*w-problem.weightedReference;
    grad=reshape(double(problem.weightedD'*residual),[],1);gradEvals=gradEvals+1;
    score=max(-grad,0).^2./max(2*problem.columnNormSquared,eps);score(active)=-Inf;
    eligible=find(isfinite(score)&score>settings.minimumAdmissionScore);
    [sorted,order]=sort(score(eligible),'descend');
    capacity=max(0,settings.maximumSupportSize-nnz(active));
    take=min([settings.batchSize,numel(order),capacity]);add=false(n,1);
    if take>0,add(eligible(order(1:take)))=true;cutoff=sorted(take);else,cutoff=NaN;end
    weightedRelativeError=norm(residual)/max(norm(problem.weightedReference),eps);
    dose=problem.doseFunction(w);clinicalObjective=problem.objectiveFunction(w);objEvals=objEvals+1;
    gate=problem.gateFunction(weightedRelativeError,dose,clinicalObjective);
    h=empty_history();h.outerIteration=k;h.supportBefore=nnz(before);
    h.activeSpots=nnz(active);h.activeLayers=count_layers(w,problem.groups,settings.zeroTolerance);
    h.addedSpots=nnz(add);h.removedZeroSpots=nnz(removable);
    h.eligibleColumns=numel(eligible);h.rejectedEligibleColumns=numel(eligible)-take;
    h.reconstructionObjective=0.5*(residual'*residual);
    h.weightedRelativeError=weightedRelativeError;h.clinicalObjective=clinicalObjective;
    h.targetD95=dose.summary.targetD95Min;h.targetD98=dose.summary.targetD98Min;
    h.oarDmax=dose.summary.oarDmaxMax;h.representative=gate.pass;
    h.gateReason=gate.reason;h.maximumAdmissionScore=max([0;score(isfinite(score))]);
    h.admissionCutoff=cutoff;h.nnlsIterations=nnls.iterations;h.nnlsReason=nnls.reason;
    h.nnlsStationarity=last_stationarity(nnls);h.objectiveEvaluations=objEvals;
    h.gradientEvaluations=gradEvals;h.elapsedSeconds=toc(timer);h.supportHash=support_hash(active);
    history(k)=h;checkpoint(k,nnls,gate);
    fprintf(['CORESET k=%d spots=%d layers=%d add=%d remove=%d relErr=%.4e ' ...
        'D95=%.4f D98=%.4f clinicalObj=%.6e score=%.3e gate=%d NNLS=%s\n'], ...
        k,h.activeSpots,h.activeLayers,h.addedSpots,h.removedZeroSpots, ...
        h.weightedRelativeError,h.targetD95,h.targetD98,h.clinicalObjective, ...
        h.maximumAdmissionScore,h.representative,h.nnlsReason);
    if gate.pass,stopReason='representativeness_gate_passed';break;end
    if k==settings.maxOuterIterations
        % Do not place a priced-but-unsolved admission batch in the final
        % artifact. The current w/active pair is the last solved NNLS state.
        stopReason='maximum_outer_iterations';break;
    end
    if take==0
        if capacity==0,stopReason='maximum_support_without_representativeness';
        else,stopReason='no_positive_reconstruction_column';end
        break;
    end
    active=active|add;w(add)=problem.referenceW(add);
end
result=state([]);save(fullfile(settings.outputDir,'coreset_final.mat'),'result','-v7.3');
    function checkpoint(iter,nnlsResult,gateResult)
        stateValue=state(nnlsResult);stateValue.lastGate=gateResult; %#ok<NASGU>
        save(fullfile(settings.outputDir,sprintf('coreset_iter_%03d.mat',iter)), ...
            'stateValue','-v7.3');
        writetable(struct2table(history(1:used)),fullfile(settings.outputDir,'coreset_trajectory.csv'));
    end
    function out=state(lastNnls)
        out=struct('w',w,'support',active,'history',history(1:used), ...
            'iterations',used,'stopReason',stopReason,'wallSeconds',toc(timer), ...
            'objectiveEvaluations',objEvals,'gradientEvaluations',gradEvals, ...
            'lastNnls',lastNnls,'settings',settings);
    end
end
function x=last_stationarity(n),if isempty(n.history),x=Inf;else,x=n.history(end).normalizedStationarity;end,end
function n=count_layers(w,g,tol),n=0;for q=1:numel(g),n=n+any(w(g(q).ix)>tol);end,end
function h=empty_history()
h=struct('outerIteration',NaN,'supportBefore',NaN,'activeSpots',NaN, ...
    'activeLayers',NaN,'addedSpots',NaN,'removedZeroSpots',NaN, ...
    'eligibleColumns',NaN,'rejectedEligibleColumns',NaN, ...
    'reconstructionObjective',NaN,'weightedRelativeError',NaN, ...
    'clinicalObjective',NaN,'targetD95',NaN,'targetD98',NaN,'oarDmax',NaN, ...
    'representative',false,'gateReason','','maximumAdmissionScore',NaN, ...
    'admissionCutoff',NaN,'nnlsIterations',NaN,'nnlsReason','', ...
    'nnlsStationarity',NaN,'objectiveEvaluations',NaN,'gradientEvaluations',NaN, ...
    'elapsedSeconds',NaN,'supportHash','');
end
function value=support_hash(mask)
md=java.security.MessageDigest.getInstance('SHA-256');md.update(uint8(mask(:)));
bytes=typecast(md.digest(),'uint8');value=lower(reshape(dec2hex(bytes,2).',1,[]));
end
