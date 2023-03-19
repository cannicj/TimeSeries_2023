% Parameters
n = 50; % Sample size
p = 3; % Number of regressors
sim = 1e5; % Number of simulations
X = [ones(n, 1), randn(n, p)]; % Design matrix
beta0 = [0;0;0;0]; % True beta vector under null hypothesis
beta1 = [1;1;1;1]; % True beta vector under alternative hypothesis
sigmas = [1, 2, 4]; % Residual variances


%%
% Simulate for test one
R2_test_one = zeros(sim, 3);
F_test_one = zeros(sim, 3);
for j = 1:length(sigmas)
    sigma = sigmas(j);
    e = normrnd(0,sqrt(sigma),n,sim); % Residual vector
    for i = 1:sim
        y = X*beta0 + e(:,i); % Response vector
        % OLS regression
        betaHat = inv((X' * X))*X'*y;
        yHat = X * betaHat;
        RSS = sum((y - yHat).^2);
        TSS = sum((y - mean(y)).^2);
        R2 = 1 - RSS/TSS; % save here for clearer use later on
        R2_test_one(i,j) = R2;
        F_test_one(i,j) = ((n-p+1)/p)*(R2/(1-R2));
    end
end

%%
% Simulate for test two
R2_test_two = zeros(sim, 3);
F_test_two = zeros(sim, 3);
for j = 1:length(sigmas)
    sigma = sigmas(j);
    e = normrnd(0,sqrt(sigma),n,sim); % Residual vector
    for i = 1:sim
        y = X*beta1 + e(:,i); % Response vector
        % OLS regression
        betaHat = inv((X' * X))*X'*y;
        yHat = X * betaHat;
        RSS = sum((y - yHat).^2);
        TSS = sum((y - mean(y)).^2);
        R2 = 1 - RSS/TSS; % save here for clearer use later on
        R2_test_two(i,j) = R2;
        F_test_two(i,j) = ((n-p+1)/p)*(R2/(1-R2));
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
    plot(x,betapdf(x,p/2,(n-p+1)/2),'--','LineWidth',2)
    plot(x,fpdf(x,p,n-p+1),'--','LineWidth',2)
    % Add labels and legend
    xlabel('Statistic Value')
    ylabel('Density')
    title(['Simulated vs. Theoretical Distributions for \sigma^2 = ' num2str(sigmas(s))])
    legend('R-squared','F','Null R-squared','Null F','Location','NorthEast')
    hold off
 end

 %%
 % Plot kernel density estimate of simulated distribution
 for c=1:2
     figure
     %hold on
     if c==1
         x = linspace(0,1,1000);
         for s=1:3
             hold on
             yR2 = ksdensity(R2_test_two(:,s),x);
             plot(x,yR2,'LineWidth',2)
             hold off
         end
         % Add labels and legend
         xlabel('Statistic Value')
         ylabel('Density')
         title('Simulated R-squared distributions')
         legend(['\sigma^2 = ' num2str(sigmas(1))],['\sigma^2 = ' num2str(sigmas(2))],['\sigma^2 = ' num2str(sigmas(3))], 'Location','NorthWest')
     else
         x = linspace(0,100,1000);
         for s=1:3
             hold on
             yF = ksdensity(F_test_two(:,s),x);
             plot(x,yF,'LineWidth',2)
             hold off
         end
         % Add labels and legend
         xlabel('Statistic Value')
         ylabel('Density')
         title('Simulated F distributions')
         legend(['\sigma^2 = ' num2str(sigmas(1))],['\sigma^2 = ' num2str(sigmas(2))],['\sigma^2 = ' num2str(sigmas(3))], 'Location','NorthWest')
     end
    %hold off
end
