function result=solve_bidirectional_repair(fullModel,wDefect,correctionIx,settings)
%SOLVE_BIDIRECTIONAL_REPAIR Exact full-objective solve on bounded corrections.
narginchk(4,4);
wDefect=reshape(double(wDefect),[],1);correctionIx=unique(reshape(double(correctionIx),[],1),'stable');
assert(numel(wDefect)==fullModel.nWeights&&~isempty(correctionIx));
assert(all(correctionIx>=1&correctionIx<=numel(wDefect)));
if fullModel.hasConstraints
    error('solve_bidirectional_repair:constraintsNeedSplit', ...
        'This frozen case unexpectedly contains hard constraints; do not silently drop them.');
end
y0=wDefect(correctionIx);callback=struct('objective',0,'gradient',0);
restricted=struct('nWeights',numel(correctionIx), ...
    'objectiveFunction',@objective,'objectiveGradient',@gradient);
result=solve_weighted_matrad_ipopt(restricted,y0,zeros(size(y0)),settings);
result.wFull=compose(result.w);result.correctionIx=correctionIx;
result.fullCallbackCounts=callback;

    function w=compose(y)
        w=wDefect;w(correctionIx)=reshape(double(y),[],1);
    end
    function f=objective(y)
        callback.objective=callback.objective+1;
        f=evaluate_matrad_dose_model(fullModel,compose(y),false);
    end
    function g=gradient(y)
        callback.gradient=callback.gradient+1;
        [~,gFull]=evaluate_matrad_dose_model(fullModel,compose(y),true);
        g=gFull(correctionIx);
    end
end
