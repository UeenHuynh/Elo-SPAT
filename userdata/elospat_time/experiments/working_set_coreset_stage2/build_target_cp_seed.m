function [mask,diagnostics]=build_target_cp_seed(D,dRef,cst,map,perRegionCp)
%BUILD_TARGET_CP_SEED Target-region/CP-balanced quadratic-improvement seed.
narginchk(5,5);n=size(D,2);mask=false(n,1);cp=reshape(double([map.cpIndex]),[],1);
targetRows=find(strcmpi(string(cst(:,3)),'TARGET'));assert(~isempty(targetRows));
% Higher-dose targets are made exclusive first; lower-dose targets retain
% their surrounding region. If subtraction empties a target, use it whole.
prescription=zeros(numel(targetRows),1);
for k=1:numel(targetRows)
    fs=cst{targetRows(k),6};vals=[];
    for j=1:numel(fs)
        try
            if iscell(fs),item=fs{j};else,item=fs(j);end
            parameters=item.getDoseParameters();
            vals=[vals;double(parameters(:))]; %#ok<AGROW>
        catch
            % Ordering is only used to make nested targets exclusive from
            % high to low prescription. Missing legacy metadata falls back
            % to stable cst target order without changing candidate access.
        end
    end
    if ~isempty(vals),prescription(k)=max(vals);end
end
[~,order]=sort(prescription,'descend');claimed=[];
records=repmat(struct('target','','targetRegionVoxels',0,'cpIndex',0, ...
    'selectedSpots',0,'minimumScore',NaN),0,1);
for oo=1:numel(order)
    row=targetRows(order(oo));whole=unique(double(cst{row,4}{1}(:)));
    region=setdiff(whole,claimed);if isempty(region),region=whole;end;claimed=union(claimed,whole);
    Ar=D(region,:);b=dRef(region);norm2=full(sum(Ar.^2,1))';corr=reshape(double(Ar'*b),[],1);
    score=max(corr,0).^2./max(2*norm2,eps);
    cpOrder=unique(cp,'stable');
    for q=1:numel(cpOrder)
        ix=find(cp==cpOrder(q));[~,rank]=sort(score(ix),'descend');take=min(perRegionCp,numel(ix));
        chosen=ix(rank(1:take));mask(chosen)=true;
        records=[records;struct('target',char(string(cst{row,2})), ... %#ok<AGROW>
            'targetRegionVoxels',numel(region),'cpIndex',cpOrder(q), ...
            'selectedSpots',take,'minimumScore',score(chosen(end)))];
    end
end
diagnostics=struct2table(records);
end
