function report=validate_working_set_callbacks(dij,cst,pln,wReference,outputFile)
%VALIDATE_WORKING_SET_CALLBACKS Restricted/full objective-gradient checks.
narginchk(5,5); n=double(dij.totalNumOfBixels);
w=reshape(double(wReference),[],1);assert(numel(w)==n&&all(w>=0));
fullModel=build_matrad_dose_model(dij,cst,pln,0);

% Deterministic representative support: strongest 200 reference weights.
[~,order]=sort(w,'descend');mask=false(n,1);mask(order(1:min(200,n)))=true;
[dijRestricted,restrictedToFull]=reduce_dij_to_support(dij,mask);
restrictedModel=build_matrad_dose_model(dijRestricted,cst,pln,0);
wRestricted=w(restrictedToFull);wEmbedded=zeros(n,1);wEmbedded(restrictedToFull)=wRestricted;
fFull=fullModel.optiProb.matRad_objectiveFunction(wEmbedded,fullModel.dij,fullModel.cst);
fRestricted=restrictedModel.optiProb.matRad_objectiveFunction( ...
    wRestricted,restrictedModel.dij,restrictedModel.cst);
gFull=reshape(double(fullModel.optiProb.matRad_objectiveGradient( ...
    wEmbedded,fullModel.dij,fullModel.cst)),[],1);
gRestricted=reshape(double(restrictedModel.optiProb.matRad_objectiveGradient( ...
    wRestricted,restrictedModel.dij,restrictedModel.cst)),[],1);
objectiveRelativeError=abs(fFull-fRestricted)/max(1,abs(fFull));
gradientRelativeError=norm(gFull(restrictedToFull)-gRestricted)/max(1,norm(gRestricted));
gReference=reshape(double(fullModel.optiProb.matRad_objectiveGradient( ...
    w,fullModel.dij,fullModel.cst)),[],1);

% Directional derivative checks on fixed, spread coordinates.
indices=unique(round(linspace(1,n,min(10,n))))';
rows=repmat(struct('index',NaN,'analytic',NaN,'finiteDifference',NaN, ...
    'relativeError',NaN),numel(indices),1);
for q=1:numel(indices)
    j=indices(q);h=1e-6*max(1,abs(w(j)));wp=w;wp(j)=wp(j)+h;
    fp=fullModel.optiProb.matRad_objectiveFunction(wp,fullModel.dij,fullModel.cst);
    if w(j)>h
        wm=w;wm(j)=wm(j)-h;
        fm=fullModel.optiProb.matRad_objectiveFunction(wm,fullModel.dij,fullModel.cst);
        fd=(fp-fm)/(2*h);
    else
        f0=fullModel.optiProb.matRad_objectiveFunction(w,fullModel.dij,fullModel.cst);
        fd=(fp-f0)/h;
    end
    rows(q).index=j;rows(q).analytic=gReference(j);rows(q).finiteDifference=fd;
    rows(q).relativeError=abs(fd-gReference(j))/max([1,abs(fd),abs(gReference(j))]);
end
writetable(struct2table(rows),outputFile);
report=struct('objectiveRelativeError',objectiveRelativeError, ...
    'gradientRelativeError',gradientRelativeError, ...
    'maximumFiniteDifferenceRelativeError',max([rows.relativeError]), ...
    'nFiniteDifferences',numel(rows));
assert(objectiveRelativeError<=1e-11,'Restricted/full objective mismatch.');
assert(gradientRelativeError<=1e-10,'Restricted/full gradient mismatch.');
assert(report.maximumFiniteDifferenceRelativeError<=1e-3, ...
    'Full matRad gradient finite-difference check failed.');
end
