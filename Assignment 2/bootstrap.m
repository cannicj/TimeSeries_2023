% Exercise 1) Part 3: Creating Confidence Intervals based on 
% the method of single parametric bootstrap
% initialize parameters
true_b = .2; sig=1; T=1000; sim=100;
B = 400;
b1=zeros(sim); b2=b1; b3=b1;
Mdl = arima(0,0,1);
alpha = .1;
%initialize parameters for actual coverage
bool_Durbin=zeros(sim, 1);
bool_MLE=zeros(sim, 1);
bool_matlab=zeros(sim, 1);
%initialize parameters for length of confidence intervals
length_Durbin=zeros(sim, 1);
length_MLE=zeros(sim, 1);
length_matlab=zeros(sim, 1);

for i=1:sim %loop to get average CI coverage
    y=ma1sim(T,sig,true_b); % simulation of the process
    b1=DurbinMA1959(y,1); % estimation by DurbinMA1959 method
    param=ma1(y,2); b2=param(1); % estimation by MP MLE method
    b3=vertcat(estimate(Mdl,y,'Display','off').MA{:}); % estimation by builtin method
    %initialize parameters simulation of bootstrap data
    vec_timeseries_Durbin = zeros(T, B);
    vec_timeseries_approxMLE = zeros(T, B);
    vec_timeseries_matlab = zeros(T, B);
    %initialize parameters for length of confidence intervals
    boot_MA1est_Durbin = zeros(B, 1);
    boot_MA1est_approxMLE = zeros(B, 1);
    boot_MA1est_matlab = zeros(B, 1);
    if mod(i,10)==0, disp(['Iteration:' num2str(i)]); end
    for j = 1:B
        % simulation of the MA(1)-process with the estimated b based on all
        % three methods
        vec_timeseries_Durbin(:,j) = ma1sim(T, sig, b1); %simulate bootstrap data with durbin estimation param
        vec_timeseries_approxMLE(:,j) = ma1sim(T, sig, b2); %simulate bootstrap data with MP MLE estimation param
        vec_timeseries_matlab(:,j) = ma1sim(T, sig, b3); %simulate bootstrap data with builtin estimation param
        % parameter estimation based from bootstrap samples
        boot_MA1est_Durbin(j) = DurbinMA1959(vec_timeseries_Durbin(:,j), 1);
        param=ma1(vec_timeseries_approxMLE(:,j),2); boot_MA1est_approxMLE(j) = param(1);
        boot_MA1est_matlab(j) = vertcat(estimate(Mdl,vec_timeseries_matlab(:,j),'Display','off').MA{:});
    end
    % calculate CI
    ci_Durbin = quantile(boot_MA1est_Durbin, [alpha/2 1-alpha/2]);
    ci_approxMLE = quantile(boot_MA1est_approxMLE, [alpha/2 1-alpha/2]);
    ci_matlab = quantile(boot_MA1est_matlab, [alpha/2 1-alpha/2]);
    % store length of the confidence intervals
    length_Durbin(i)=ci_Durbin(2) -ci_Durbin(1);
    length_MLE(i)=ci_approxMLE(2) - ci_approxMLE(1);
    length_matlab(i)=ci_matlab(2) - ci_matlab(1);
    % store actual coverage
    bool_Durbin(i)=(true_b>ci_Durbin(1)) & (true_b<ci_Durbin(2));
    bool_MLE(i)=(true_b>ci_approxMLE(1)) & (true_b<ci_approxMLE(2));
    bool_matlab(i)=(true_b>ci_matlab(1)) & (true_b<ci_matlab(2));
end
%%
disp(["Actual coverage matlab method:", num2str(mean(bool_matlab))])
disp(["Actual coverage MLE method:", num2str(mean(bool_MLE))])
disp(["Actual coverage Durbin method:", num2str(mean(bool_Durbin))])

disp(["CI Length matlab method:", num2str(mean(length_matlab))])
disp(["CI Length MLE method:", num2str(mean(length_MLE))])
disp(["CI Length Durbin method:", num2str(mean(length_Durbin))])

%%
% calculate the mean of it to get 




n=20; p=0.3; B=1e4; alpha =0.05; sim=1e5; bool=zeros(sim, 1);
for i = 1:sim
    phat0 = binornd(n, p, [1, 1])/n; % the estimate of p from Bin(n, p) data
    phatvec = binornd(n, phat0, [B, 1])/n; % B samples of S/n , S~Bin(n, phat0)
    ci = quantile(phatvec, [alpha/2 1-alpha/2]); low = ci(1); high = ci(2);
    bool(i) = (p>low) & (p<high) ; % is the true p in the interval?
end
actualcoverage = mean(bool)

%% aaron codes haha
% prep work
b_true = -.5; % freely assumed
T=100;
n_obs = T;% + p;
B = 400; % check if can be smaller or should be larger for reliable results
alpha = .1;

% initialize variables
vec_timeseries_Durbin = zeros(n_obs, B);
vec_timeseries_approxMLE = zeros(n_obs, B);
vec_timeseries_matlab = zeros(n_obs, B);
boot_MA1est_Durbin = zeros(B, 1);
boot_MA1est_approxMLE = zeros(B, 1);
boot_MA1est_matlab = zeros(B, 1);

% simulate MA(1) process and then estimate the parameter
% % for reproducibility
rng(42); seed = rng;

% % simulate MA(1) process with true parameter
timeseries = armasim(n_obs, 1, 0, b_true);

% % estimate MA parameter using the three methods
% %% Durbin
MA1est_Durbin = DurbinMA1959(timeseries, 1);
% %% approx MLE
MA1_temp = ma1(timeseries, 0); % 0 for approx MLE
MA1est_approxMLE = MA1_temp(1);
% %% matlab
MA1_temp = armax(timeseries, [0, 1]);
MA1est_matlab = MA1_temp.c(2:end);

for i = 1:B
    rng(seed.Seed + i);
    % simulation of MA process with ESTIMATED b
    vec_timeseries_Durbin(:,i) = armasim(n_obs, 1, 0, MA1est_Durbin);
    vec_timeseries_approxMLE(:,i) = armasim(n_obs, 1, 0, MA1est_approxMLE);
    vec_timeseries_matlab(:,i) = armasim(n_obs, 1, 0, MA1est_matlab);

    % estimation of bootstrapped timeseries
    boot_MA1est_Durbin(i) = DurbinMA1959(vec_timeseries_Durbin(:,i), 1);
    MA1_temp = ma1(vec_timeseries_approxMLE(:,i), 0);
    boot_MA1est_approxMLE(i) = MA1_temp(1);
    MA1_temp = armax(vec_timeseries_matlab(:,i), [0, 1]);
    boot_MA1est_matlab(i) = MA1_temp.c(2:end);
end

% calculate CI
ci_Durbin = quantile(boot_MA1est_Durbin, [alpha/2 1-alpha/2]);
ci_approxMLE = quantile(boot_MA1est_approxMLE, [alpha/2 1-alpha/2]);
ci_matlab = quantile(boot_MA1est_matlab, [alpha/2 1-alpha/2]);