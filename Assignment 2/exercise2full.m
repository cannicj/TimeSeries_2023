% exercise 2, loop for results
%% Define ARMA(3,2) model parameters
a = [0.4, -0.5, -0.2]'; % AR coefficients
b = [-0.5, -0.24]; % MA coefficients
sigma = 1; % noise standard deviation
Tvec = [100 1000]; % length of time series, need 100 and 1000
rep = 10000;
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
        if mod(i,50)==0, disp([num2str(i/100) "% done."]); end
        % Generate ARMA(3,2) data
        Y = armasim(T, sigma, a, b);
        %Part b: Estimate optimal AR(p) using AIC and BIC
        AICminAR=1e20; BICminAR=1e20;
        for p=1:pmax
            [ar_est, logL] = exactarp(Y,p);
            AIC=2*(p+1) - 2*(-logL); BIC = log(T)*(p+1)-2*(-logL);
            if AIC<AICminAR, AICminAR=AIC; pvecAIC(i,t)=p; end
            if BIC<BICminAR, BICminAR=BIC; pvecBIC(i,t)=p; end
        end

        % Part c: Estimate optimal MA(q) using AIC and BIC
        AICminMA=1e20; BICminMA=1e20;
        for q = 1:qmax
            %[ma_est, logL] = maq(Y, 2, q);
            ma_est = DurbinMA1959(Y,q);
            logL = condmaq([ma_est' 1]', Y);
            AIC=2*q - 2*(-logL); BIC = log(T)*q-2*(-logL);
            if AIC<AICminMA, AICminMA=AIC; qvecAIC(i,t)=q; end
            if BIC<BICminMA, BICminMA=BIC; qvecBIC(i,t)=q; end
        end
    end
end
toc
%%
qvecAIC10000 = array2table(qvecAIC);
qvecBIC10000 = array2table(qvecBIC);
pvecAIC10000 = array2table(pvecAIC);
pvecBIC10000 = array2table(pvecBIC);
writetable(qvecAIC10000,'qvecAIC10000.csv', 'WriteVariableNames', true);
writetable(qvecBIC10000,'qvecBIC10000.csv', 'WriteVariableNames', true);
writetable(pvecAIC10000,'pvecAIC10000.csv', 'WriteVariableNames', true);
writetable(pvecBIC10000,'pvecBIC10000.csv', 'WriteVariableNames', true);
%%

qvecAIC10000 = readtable('qvecAIC10000.csv');
qvecBIC10000 = readtable('qvecBIC10000.csv');
pvecAIC10000 = readtable('pvecAIC10000.csv');
pvecBIC10000 = readtable('pvecBIC10000.csv');
%%
qvecAIC10000 = qvecAIC10000{:,:};
qvecBIC10000 = qvecBIC10000{:,:};
pvecAIC10000 = pvecAIC10000{:,:};
pvecBIC10000 = pvecBIC10000{:,:};
