function [values,J,rows]=evaluate_clinical_gates(D,w,spec,needJacobian)
% Exact metric values and a valid local gradient/subgradient row per gate.
if nargin<4,needJacobian=nargout>=2;end
d=D*w;m=numel(spec);values=zeros(m,1);J=sparse(m,size(D,2));
rows=repmat(struct('structure','','metric','','sense','','value',NaN,'limit',NaN,'margin',NaN,'pass',false),m,1);
for i=1:m
    vox=spec(i).vox;x=reshape(double(d(vox)),[],1);
    if strcmp(spec(i).metric,'Dmax')
        values(i)=max(x);if needJacobian,tol=1e-10*max(1,abs(values(i)));tie=vox(abs(x-values(i))<=tol);J(i,:)=mean(D(tie,:),1);end
    elseif strcmp(spec(i).metric,'Dmean')
        values(i)=mean(x);if needJacobian,J(i,:)=mean(D(vox,:),1);end
    else
        [values(i),lo,hi,a]=percentile_active(x,spec(i).percentile);if needJacobian,J(i,:)=(1-a)*tie_row(D,vox,x,x(lo))+a*tie_row(D,vox,x,x(hi));end
    end
    if strcmp(spec(i).sense,'lower'),margin=values(i)-spec(i).limit;else,margin=spec(i).limit-values(i);end
    rows(i)=struct('structure',spec(i).structure,'metric',spec(i).metric,'sense',spec(i).sense, ...
        'value',values(i),'limit',spec(i).limit,'margin',margin,'pass',margin>=-1e-7);
end
end
function r=tie_row(D,vox,x,value),tol=1e-10*max(1,abs(value));ix=vox(abs(x-value)<=tol);r=mean(D(ix,:),1);end
function [p,lo,hi,a]=percentile_active(x,q)
[x,~]=sort(x,'ascend');n=numel(x);pos=1+(n-1)*q/100;lo=floor(pos);hi=ceil(pos);a=pos-lo;p=(1-a)*x(lo)+a*x(hi);
end
