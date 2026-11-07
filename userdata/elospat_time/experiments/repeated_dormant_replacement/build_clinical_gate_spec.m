function spec=build_clinical_gate_spec(D,cst,referenceDose,pctTol,absTol)
% Exact project acceptance gates, separate from the matRad objective.
D=sparse(double(D));spec=repmat(struct('structure','','metric','','sense','', ...
    'limit',NaN,'vox',[],'percentile',NaN,'jacobianStructure',sparse(1,size(D,2))),0,1);
for i=1:numel(referenceDose.structures)
    r=referenceDose.structures(i);ci=find(strcmp(cst(:,2),r.name),1);vox=reshape(double(cst{ci,4}{1}),[],1);
    if r.isTarget
        keys={'D95',5;'D98',2};
        for k=1:size(keys,1),s=empty(r.name,keys{k,1},'lower',(1-pctTol/100)*r.(keys{k,1}),vox,keys{k,2},D);spec(end+1,1)=s;end %#ok<AGROW>
    elseif strcmpi(r.type,'OAR')&&~strcmpi(r.name,'BODY')&&~strcmpi(r.name,'EXTERNAL')
        keys={'Dmean','mean';'Dmax','max'};
        for k=1:size(keys,1),m=keys{k,1};lim=r.(m)+max(absTol,pctTol/100*max(abs(r.(m)),1));s=empty(r.name,m,'upper',lim,vox,keys{k,2},D);spec(end+1,1)=s;end %#ok<AGROW>
    end
end
end
function s=empty(name,metric,sense,limit,vox,kind,D)
p=NaN;if isnumeric(kind),p=kind;end
s=struct('structure',name,'metric',metric,'sense',sense,'limit',limit,'vox',vox, ...
    'percentile',p,'jacobianStructure',spones(sum(D(vox,:),1)));
end
