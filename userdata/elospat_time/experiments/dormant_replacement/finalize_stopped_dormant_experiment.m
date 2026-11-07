function finalize_stopped_dormant_experiment(resultDir)
%FINALIZE_STOPPED_DORMANT_EXPERIMENT Read-only report from event-1 checkpoints.
stage1File='/mnt/beegfs/ydang/project1/data/Elo-SPAT/userdata/elospat_time/results/full_integrated_layer_20deg_job22180/alpha_0p3/alpha_0p3_final.mat';
stage1=load(stage1File,'finalDose');
s=load(fullfile(resultDir,'event_01_survivor.mat'),'checkpoint');
a=load(fullfile(resultDir,'event_01_scalar_gradient.mat'),'rescueCheckpoint');
b=load(fullfile(resultDir,'event_01_component_aware.mat'),'rescueCheckpoint');
names={'survivor','scalar_gradient','component_aware'};
plans={s.checkpoint.survivor.plan,a.rescueCheckpoint.solution.plan,b.rescueCheckpoint.solution.plan};
rows=repmat(empty_row(),0,1);
for p=1:numel(plans)
    d=plans{p}.dose;
    for i=1:numel(d.structures)
        x=d.structures(i);r=stage1.finalDose.structures(i);assert(strcmp(x.name,r.name));
        q=empty_row();q.event=1;q.method=names{p};q.structure=x.name;q.type=x.type;
        q.stage1D98=r.D98;q.trialD98=x.D98;q.stage1D95=r.D95;q.trialD95=x.D95;
        q.stage1Dmean=r.Dmean;q.trialDmean=x.Dmean;q.stage1Dmax=r.Dmax;q.trialDmax=x.Dmax;
        rows(end+1,1)=q; %#ok<AGROW>
    end
end
writetable(struct2table(rows,'AsArray',true),fullfile(resultDir,'clinical_plan_metrics.csv'));
status=table({'DORMANT_REPLACEMENT_NOT_SUPPORTED'},{'replacement_search_cost_dominated_and_no_accepted_birth_event'}, ...
    1,0,3,56,677,49,'VariableNames',{'status','stopReason','completedEvents', ...
    'acceptedReplacementEvents','observedHistoricalBirthsInRejectedTrial', ...
    'sameEventReturnsInRejectedTrial','proposedDeaths','initialDormantPool'});
writetable(status,fullfile(resultDir,'dormant_replacement_status.csv'));
fprintf('STOPPED_DORMANT_REPORT_COMPLETE\n');disp(status);
end

function r=empty_row()
r=struct('event',0,'method','','structure','','type','','stage1D98',NaN, ...
    'trialD98',NaN,'stage1D95',NaN,'trialD95',NaN,'stage1Dmean',NaN, ...
    'trialDmean',NaN,'stage1Dmax',NaN,'trialDmax',NaN);
end
