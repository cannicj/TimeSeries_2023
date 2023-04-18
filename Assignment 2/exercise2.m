% exercise 2, one execution
%% Define ARMA(3,2) model parameters
a = [0.2, -0.3, 0.5]'; % AR coefficients
b = [0.4, -0.6]; % MA coefficients
sigma = 1; % noise standard deviation
T = 100; % length of time series, need 100 and 1000
rep = 1e4;

% Generate ARMA(3,2) data
Y = armasim(T, sigma, a, b);

%% Part b: Estimate optimal AR(p) using AIC and BIC
AICminAR=1e20; BICminAR=1e20;
pmax = 10;
%pvecAIC=zeros(rep,1);
%pvecBIC=zeros(rep,1);
for p=1:pmax
    ar_est = exactarp(Y,p);
    %param=armareg(Y,X,p,0,1);
    lsig2=log(ar_est(end)^2);
    z=p; AIC=lsig2+(2*z)/T; BIC=lsig2 + z*log(T)/T;
    disp(["AIC", num2str(AIC)])
    if AIC<AICminAR, AICminAR=AIC; chosenpAIC = p; end %pvecAIC(p)=p; end
    if BIC<BICminAR, BICminAR=BIC; chosenpBIC = p;end %pvecBIC(p)=p; end
end

%%
%for p = 1:pmax
%    Mdl = arima(p,0,0);
%    [est,~,logL] = estimate(Mdl,Y,'display','off');
%    aic_p(p) = -2*logL + 2*p;
%    bic_p(p) = -2*logL + log(N)*p;
%end
fprintf('AIC selected AR(p) = %d\n', opt_aic_p);
fprintf('BIC selected AR(p) = %d\n', opt_bic_p);

%% Part c: Estimate optimal MA(q) using AIC and BIC
AICminMA=1e20; BICminMA=1e20;
qmax = 10;
%qvecAIC=zeros(rep,1);
%qvecBIC=zeros(rep,1);
for q = 1:qmax
    %estimate MA(q) with sigma
    ma_est = maq(Y, 2, q);
    lsig2=log(ma_est(end)^2);
    z=q; AIC=lsig2+(2*z)/T; BIC=lsig2 + z*log(T)/T;
    if AIC<AICminMA, AICminMA=AIC; chosenqAIC = q; end %qvecAIC(i)=q; end
    if BIC<BICminMA, BICminMA=BIC; chosenqBIC = q; end %qvecBIC(i)=q; end
end
