bvec=-0.9:0.1:0.9; sig=1; T=100; sim=100;
b1=zeros(sim,length(bvec)); b2=b1; b3=b1; true=kron(ones(sim,1),bvec);
Mdl = arima(0,0,1);

for bloop=1:length(bvec), b=bvec(bloop);
  for i=1:sim
    y=ma1sim(T,sig,b); 
    if mod(i,100)==0, disp([i, b]); end 
    param=ma1(y,2); b2(i,bloop)=param(1); % Marc Paolella code
    b1(i,bloop)=vertcat(estimate(Mdl,y,'Display','off').MA{:}); % built in code
    b3(i,bloop)=DurbinMA1959(y,1); % DurbinMA1959
  end 
end
%%
figure, set(gca,'fontsize',16)
plot(bvec,mean(b1)-bvec,'r-',bvec,mean(b2)-bvec,'g--',bvec,mean(b3)-bvec,'b--','linewidth',3)
title('Bias'), legend('built in function','Marc Paolella code', 'DurbinMA1959')
figure, set(gca,'fontsize',16)
plot(bvec, mean((b1-true).^2), 'r-', ...
     bvec, mean((b2-true).^2), 'g--', ...
     bvec, mean((b3-true).^2), 'b--', 'linewidth',3)
title('MSE'), legend('built in function','Marc Paolella code','DurbinMA1959')
