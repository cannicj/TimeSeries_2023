% Parameters
n = 50; % Sample size
p = 3; % Number of regressors
sim = 1e5; % Number of simulations
X = [ones(n, 1), randn(n, p)]; % Design matrix
beta0 = [0; zeros(p, 1)]; % True beta vector under null hypothesis
beta1 = [1; ones(p, 1)]; % True beta vector under alternative hypothesis
sigmas = [1, 2, 4]; % Residual variances
%% 

% Simulate for test one
R2_test_one = zeros(sim, 3);
F_test_one = zeros(sim, 3);
for j = 1:length(sigmas)
    sigma = sigmas(j);
    for i = 1:sim
        e = normrnd(0,sqrt(sigma),n,1); % Residual vector
        y = X * beta0 + e; % Response vector
        % OLS regression
        [b, ~, r, ~, stats] = regress(y, X);
        % R^2 and F statistics
        R2_test_one(i, j) = stats(1);
        F_test_one(i, j) = stats(2);
    end
end


%%
% Simulate for test two
R2_test_two = zeros(sim, 3);
F_test_two = zeros(sim, 3);
for j = 1:length(sigmas)
    sigma = sigmas(j);
    for i = 1:sim
        e = normrnd(0,sqrt(sigma),n,1); % Residual vector
        y = X * beta1 + e; % Response vector
        % OLS regression
        [b, ~, r, ~, stats] = regress(y, X);
        % R^2 and F statistics
        R2_test_two(i, j) = stats(1);
        F_test_two(i, j) = stats(2);
    end
end


%%
 % Plot kernel density estimate of simulated distribution
 for s=1:3
    figure
    hold on
    x = linspace(0,1,1000);
    yR2 = ksdensity(R2_test_one(:,s),x);
    yF = ksdensity(F_test_one(:,s),x);
    plot(x,yR2,'LineWidth',2)
    plot(x,yF,'LineWidth',2)
    % Plot theoretical null distributions
    plot(x,betapdf(x,p,n-p),'--','LineWidth',2)
    plot(x,fpdf(x,p,n-p),'--','LineWidth',2)
    % Add labels and legend
    xlabel('Statistic Value')
    ylabel('Density')
    title(['Simulated vs. Theoretical Distributions for \sigma^2 = ' num2str(sigmas(s))])
    legend('R-squared','F','Null R-squared','Null F','Location','NorthWest')
    hold off
 end

 %%
 % Plot kernel density estimate of simulated distribution
 for s=1:3
    figure
    hold on
    x = linspace(0,1,1000);
    yR2 = ksdensity(R2_test_two(:,s),x);
    yF = ksdensity(F_test_two(:,s),x);
    plot(x,yR2,'LineWidth',2)
    plot(x,yF,'LineWidth',2)
    % Add labels and legend
    xlabel('Statistic Value')
    ylabel('Density')
    title(['Simulated Distributions for \sigma^2 = ' num2str(sigmas(s))])
    legend('R-squared','F','Location','NorthWest')
    hold off
end