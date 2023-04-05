T = 100; % Sample size
sigma2 = 2; % Variance of error term
k_vals = 2:10; % Number of regressors
sim = 1e4; % Number of simulations

% Preallocate array for results
R2_vals = zeros(sim, length(k_vals));

% Loop over k values
for a = 1:length(k_vals)
    k = k_vals(a); % Current k value
    % Loop over simulations
    for i = 1:sim
        % Generate design matrix X
        X = [ones(T, 1), normrnd(0, 1, T, k)];
        % Generate error term
        e = normrnd(0, sqrt(sigma2), T, 1);
        % Generate response variable
        y = X * zeros(k+1, 1) + e;
        % Fit linear regression model
        mdl = fitlm(X(:, 2:end), y);
        % Compute R-squared statistic
        R2_vals(i, a) = mdl.Rsquared.Ordinary;
    end
end

%%

% Plot boxplots of R-squared values for each k value
figure
boxplot(R2_vals, 'Labels', string(k_vals))
xlabel('Number of Regressors')
ylabel('R-squared Value')
title('Boxplots of R-squared Values for Different Numbers of Regressors')


%%

%The resulting plot shows nine boxplots, one for each k value, showing 
% the distribution of R-squared values for each simulation. From the plot, 
% we can see that the median R-squared value decreases as the number of 
% bogus regressors increases. This is because as more regressors are added,
% some of them are likely to be non-significant and add noise to the model. 
% This noise reduces the fit of the model, and as a result, the R-squared 
% value decreases. We can also see that the spread of the R-squared 
% values increases as the number of bogus regressors increases. This is
% because with more regressors, there is more variability in the model 
% fit due to the added noise from non-significant regressors. This 
% increased variability is reflected in the larger interquartile range
% and longer whiskers of the boxplots. Overall, this plot shows the 
% importance of careful selection of regressors in linear regression models. 
% Adding too many non-significant regressors can decrease the fit of the model
% and increase its variability, leading to inaccurate conclusions and predictions.
