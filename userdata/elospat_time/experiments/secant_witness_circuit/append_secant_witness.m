function [Qnew,info] = append_secant_witness(Q,a0,v,deltaH,rankTolerance)
%APPEND_SECANT_WITNESS Construct, verify and orthogonalize a finite secant.
narginchk(5,5);
a0=reshape(double(a0),[],1); v=reshape(double(v),[],1);
assert(size(Q,1)==numel(v) && numel(a0)==numel(v));
vv=v'*v;
if vv<=eps(max(1,norm(v)^2))
    error('append_secant_witness:zeroTransition','Dose transition is zero.');
end
abar=a0+((double(deltaH)-a0'*v)/vv)*v;
secantError=abs(abar'*v-double(deltaH))/max([1,abs(deltaH),abs(abar'*v)]);
if isempty(Q), qPerp=abar; else, qPerp=abar-Q*(Q'*abar); end
perpNorm=norm(qPerp); threshold=rankTolerance*max(1,norm(abar));
independent=perpNorm>threshold;
if independent, Qnew=[Q,qPerp/perpNorm]; else, Qnew=Q; end %#ok<AGROW>
info=struct('secantError',secantError,'independent',independent, ...
    'perpendicularNorm',perpNorm,'independenceThreshold',threshold, ...
    'abar',abar,'newWitness',qPerp/max(perpNorm,eps));
end
