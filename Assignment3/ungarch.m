function [eout,sigvec]=ungarch(e,c0,c1,d1)
sigvec=zeros(length(e),1); e2=e.^2;
denom=1-c1-d1; if denom>0.001, sinit=c0/denom; else sinit=mean(e2); end
einit=sinit;
% do the recursion in sigvec^delta because it is faster
sigvec(1)=c0+c1*einit+d1*sinit;
for t=2:length(e), sigvec(t)=c0 + c1 *e2(t-1) + d1*sigvec(t-1); end
sigvec=sigvec.^(1/2); eout=e./sigvec;