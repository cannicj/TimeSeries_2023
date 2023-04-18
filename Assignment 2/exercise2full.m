% exercise 2, loop for results
%% Define ARMA(3,2) model parameters
a = [0.2, -0.3, 0.5]'; % AR coefficients
b = [0.4, -0.6]; % MA coefficients
sigma = 1; % noise standard deviation
Tvec = [100 1000]; % length of time series, need 100 and 1000
rep = 100;
pmax = 10;
qmax = 10;
pvecAIC=zeros(rep,2);
pvecBIC=zeros(rep,2);
qvecAIC=zeros(rep,2);
qvecBIC=zeros(rep,2);
%%
tic
for t=1:length(Tvec)
    T = Tvec(t);
    for i=1:rep
        % Generate ARMA(3,2) data
        Y = armasim(T, sigma, a, b);
        %Part b: Estimate optimal AR(p) using AIC and BIC
        AICminAR=1e20; BICminAR=1e20;
        for p=1:pmax
            ar_est = exactarp(Y,p);
            %param=armareg(Y,X,p,0,1);
            lsig2=log(ar_est(end)^2);
            z=p; AIC=lsig2+(2*z)/T; BIC=lsig2 + z*log(T)/T;
            if AIC<AICminAR, AICminAR=AIC; pvecAIC(i,t)=p; end
            if BIC<BICminAR, BICminAR=BIC; pvecBIC(i,t)=p; end
        end

        % Part c: Estimate optimal MA(q) using AIC and BIC
        AICminMA=1e20; BICminMA=1e20;
        for q = 1:qmax
            %estimate MA(q) with sigma
            ma_est = maq(Y, 2, q);
            lsig2=log(ma_est(end)^2);
            z=q; AIC=lsig2+(2*z)/T; BIC=lsig2 + z*log(T)/T;
            if AIC<AICminMA, AICminMA=AIC; qvecAIC(i,t)=q; end
            if BIC<BICminMA, BICminMA=BIC; qvecBIC(i,t)=q; end
        end
    end
end
toc
