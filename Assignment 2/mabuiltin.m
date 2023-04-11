%Test to estimate parameters using built in estimate function

%Specify MA(1) model
DGP = arima('Constant',0.1,'MA',{0.7},...
    'MALags',[1],'Variance',0.15);

%Simulate MA(1) model
rng(5); % For reproducibility
T = 500;
y = simulate(DGP,T);

%Specify model structure we want to estimate
Mdl = arima(0,0,1);

%estimate parameters
EstMdl = estimate(Mdl,y)