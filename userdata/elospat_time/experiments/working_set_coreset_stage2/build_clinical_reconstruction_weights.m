function [weights,tableOut]=build_clinical_reconstruction_weights(cst,nVoxels)
%BUILD_CLINICAL_RECONSTRUCTION_WEIGHTS Structure-balanced voxel weights.
% Each target contributes aggregate weight 10 and each named OAR contributes
% aggregate weight 1. BODY/EXTERNAL is excluded to avoid swamping specific
% clinical structures. Overlapping structures add their contributions.
narginchk(2,2);weights=zeros(nVoxels,1);
rows=repmat(struct('structure','','type','','nVoxels',0,'aggregateWeight',0),0,1);
for s=1:size(cst,1)
    ix=unique(reshape(double(cst{s,4}{1}),[],1));
    if isempty(ix),continue;end
    name=char(string(cst{s,2}));type=char(string(cst{s,3}));
    if strcmpi(type,'TARGET'),beta=10;
    elseif strcmpi(name,'BODY')||strcmpi(name,'EXTERNAL'),beta=0;
    else,beta=1;end
    if beta>0,weights(ix)=weights(ix)+beta/numel(ix);end
    rows(end+1)=struct('structure',name,'type',type,'nVoxels',numel(ix), ...
        'aggregateWeight',beta); %#ok<AGROW>
end
positive=weights>0;assert(any(positive),'No clinically weighted voxels found.');
weights=weights/mean(weights(positive));tableOut=struct2table(rows);
end
