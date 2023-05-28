function [param,stderr,loglik,zvec] = babygarch(y)
% normal-GARCH(1,1) with power=2. y is vector of log percentage returns
initvec= [0.04 0.05 0.8 10]; % c_0 c_1 d_1 dof
bound.lo =    [0   0 0 0];
bound.hi =    [0.5 1 1 500];
bound.which = [1   1 1 1];
opt=optimset('Display','None', 'Maxiter',500, 'TolFun',1e-6, ...
    'TolX',1e-6,'LargeScale','off');
init=einschrk(initvec,bound);
[pout,~,~,~,~,hess] = fminunc(@(param) like(param,y,bound),init,opt);
[loglik,zvec]=like(pout,y,bound);
V=pinv(hess)/length(y);
[param,V]=einschrk(pout,bound,V); stderr=sqrt(diag(V));

function [loglik,zvec]=like(param,y,bound)
param=einschrk(real(param),bound,999);
c0=param(1); c1=param(2); d1=param(3); dof=param(4);
[zvec,sigvec]=ungarch(y,c0,c1,d1);
ll = gammaln((dof+1)/2) - log(sqrt(dof * pi)) - gammaln(dof / 2) - log(sigvec) - ((dof + 1)/2) * log(1 + zvec.^2/dof);
loglik = -mean(ll);

function [eout,sigvec]=ungarch(e,c0,c1,d1)
sigvec=zeros(length(e),1); e2=e.^2;
denom=1-c1-d1; if denom>0.001, sinit=c0/denom; else sinit=mean(e2); end
einit=sinit;
% do the recursion in sigvec^delta because it is faster
sigvec(1)=c0+c1*einit+d1*sinit;
for t=2:length(e), sigvec(t)=c0 + c1 *e2(t-1) + d1*sigvec(t-1); end
sigvec=sigvec.^(1/2); eout=e./sigvec;