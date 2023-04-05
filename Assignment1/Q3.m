% Q3 - AR(1) disturbances

% Parameters
samp_sizes = [20, 50];
p = 4; % Number of regressors
sim = 1e4; % Number of simulations
X = [ones(samp_sizes(1)+1, 1), randn(samp_sizes(1)+1, p)]; % Design matrix
beta0 = [0; zeros(p, 1)]; % True beta vector under null hypothesis
%beta1 = [1; ones(p, 1)]; % True beta vector under alternative hypothesis
sigmas = 1; % Residual variances
tol = 1e-5;
a_values = [0, 0.2, 0.4, 0.6, 0.8, 0.99]; % results should get worse with increasing true a


%% Tests for just one value
% simulation
e = zeros(samp_sizes(1)+1,1); % initial vector for regression error terms
U = normrnd(0,sqrt(sigmas),samp_sizes(1)+1,1); % Residual vector of the AR(1) process
e(1) = U(1);
%create AR(1) process for the error terms
for i=2:(samp_sizes(1)+1)
    e(i) = a_values(6)*e(i-1) + U(i);
end
y = X * beta0 + e; % get response vector

% estimate parameters
betaHat = inv((X' * X))*X'*y; % OLS regression for first step
yHat = X * betaHat; % get estimated y with betaHat
resid = y - yHat; % calculate the residuals from true values
% now calculate aHat based on resid (formula 4.14 in the book)
resid_t = resid(2:end); % vector of length T
resid_tminus1 = resid(1:end-1); % vector of length T
aHat = (resid_tminus1' * resid_t)/(resid_tminus1' * resid_tminus1); %estimate aHat
aHatnew=0;
betaHatnew = [0, 0, 0, 0, 0];
counter = 0;

% algorithm until convergence from book page 225
while abs(aHatnew - aHat) > tol | any(abs(betaHatnew - betaHat) > tol) 
    counter = counter + 1;
    aHat = aHatnew;
    betaHat = betaHatnew;
    a_mat = (1/(1-aHat^2))*(toeplitz([aHat.^(0:samp_sizes(1))])); %formula 4.2
    betaHatnew = inv(X'*inv(a_mat)*X)*X'*inv(a_mat)*y; %gls
    yHat = X * betaHatnew;
    resid = y - yHat;
    % calculate a_hat based on resid
    resid_t = resid(2:end); %vector of length T
    resid_tminus1 = resid(1:end-1); %vector of length T
    aHatnew = (resid_tminus1' * resid_t)/(resid_tminus1' * resid_tminus1); %estimate a hat
    if abs(aHatnew - aHat) < tol
        disp(["a converges after: " num2str(counter) "runs"])
    end
    if all(abs(betaHatnew - betaHat) < tol)
        disp(["the betas converge after: " num2str(counter) "runs"])
    end
end

%% now connect everything with all parameters and save results for beta and a
a_accuracy = zeros(sim, length(a_values));
beta_accuracy = zeros(sim, length(a_values));
for t=1:length(samp_sizes)
    X = [ones(samp_sizes(t)+1, 1), randn(samp_sizes(t)+1, p)]; % Design matrix
    for a=1:length(a_values)
        for n=1:sim
            e = zeros(samp_sizes(t)+1,1); % initial vector for regression error terms
            U = normrnd(0,sqrt(sigmas),samp_sizes(t)+1,1); % Residual vector of the AR(1) process
            e(1) = U(1);
            %create AR(1) process for the error terms
            for i=2:(samp_sizes(1)+1)
                e(i) = a_values(a)*e(i-1) + U(i);
            end
            y = X * beta0 + e; % get response vector
            % estimate the parameters
            betaHat = inv((X' * X))*X'*y; % OLS regression for first step
            yHat = X * betaHat; % get estimated y with betaHat
            resid = y - yHat; % calculate the residuals from true values
            % now calculate aHat based on resid (formula 4.14 in the book)
            resid_t = resid(2:end); % vector of length T
            resid_tminus1 = resid(1:end-1); % vector of length T
            aHat = (resid_tminus1' * resid_t)/(resid_tminus1' * resid_tminus1); %estimate aHat
            aHatnew=0;
            betaHatnew = [0, 0, 0, 0, 0];

            % algorithm until convergence from book page 225
            while abs(aHatnew - aHat) > tol | any(abs(betaHatnew - betaHat) > tol)
                aHat = aHatnew;
                betaHat = betaHatnew;
                a_mat = (1/(1-aHat^2))*(toeplitz([aHat.^(0:samp_sizes(t))])); %formula 4.2
                betaHatnew = inv(X'*inv(a_mat)*X)*X'*inv(a_mat)*yHat; %gls
                yHat = X * betaHatnew;
                resid = y - yHat;
                % calculate a_hat based on resid
                resid_t = resid(2:end); %vector of length T
                resid_tminus1 = resid(1:end-1); %vector of length T
                aHatnew = (resid_tminus1' * resid_t)/(resid_tminus1' * resid_tminus1); %estimate a hat
            end
            a_accuracy(n,a) = aHatnew - a_values(a);
            beta_accuracy(n,a) = mean(abs(betaHat));
        end
    end
    a_accuracy1 = array2table(a_accuracy);
    beta_accuracy1 = array2table(beta_accuracy);
    writetable(a_accuracy1,['a_accuracy_T' num2str(samp_sizes(t)) '.csv'], 'WriteVariableNames', true);
    writetable(beta_accuracy1,['beta_accuracy_T' num2str(samp_sizes(t)) '.csv'], 'WriteVariableNames', true);
end


%% graphs / tables
a_accuracy_T20 = readtable('a_accuracy_T20.csv');
a_accuracy_T50 = readtable('a_accuracy_T50.csv');
beta_accuracy_T20 = readtable('beta_accuracy_T20.csv');
beta_accuracy_T50 = readtable('beta_accuracy_T50.csv');
a_accuracy_T20 = a_accuracy_T20{:,:};
a_accuracy_T50 = a_accuracy_T50{:,:};
beta_accuracy_T20 = beta_accuracy_T20{:,:};
beta_accuracy_T50 = beta_accuracy_T50{:,:};
%%

figure
boxplot(beta_accuracy_T20, 'Labels', string(a_values))
xlabel('a')
ylabel('Mean distance to true \beta')
title('Mean distance to true \beta values for different a, T=20')

figure
boxplot(beta_accuracy_T50, 'Labels', string(a_values))
xlabel('a')
ylabel('Mean distance to true \beta')
title('Mean distance to true \beta values for different a, T=50')

figure
boxplot(a_accuracy_T20, 'Labels', string(a_values))
xlabel('a')
ylabel('Distance to true a')
title('Distance to true a-values for different a, T=20')

figure
boxplot(a_accuracy_T50, 'Labels', string(a_values))
xlabel('a')
ylabel('Distance to true a')
title('Distance to true a-values for different a, T=50')

mean(beta_accuracy_T20)
mean(beta_accuracy_T50)




%% visualization and tables
mean(median(a_accuracy_T50)-median(a_accuracy_T20))

