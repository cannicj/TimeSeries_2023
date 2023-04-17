function y=ma2sim(nobs,sigma,b1,b2)
 u=sigma*randn(nobs+2,1); 
 y=zeros(nobs,1); 

 for i=3:nobs+2
     y(i-2) = u(i) + b1*u(i-1) + b2*u(i-2);
 end
