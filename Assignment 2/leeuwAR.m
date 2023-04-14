function [Vi,detVi]=leeuwAR(a,T);
% a is [a1 a2 ... ap] of a stationary AR(p) model
p=length(a); a=-a; if nargin < 2, T=p+1; end
firrowP = [1 zeros(1,T-1)]; fircolP = [1 a zeros(1,T-p-1)];
P = toeplitz(fircolP,firrowP); P1 = P(1:p,1:p);
firrowQ1 = a(p:-1:1); fircolQ1 = [a(p) zeros(1,p-1)];
Q1 = toeplitz(fircolQ1,firrowQ1); Q = [Q1; zeros(T-p,p)];
Vi = P'*P - Q*Q';
if nargout>=2, detVi = det ( P1'*P1 - Q1*Q1' ); end