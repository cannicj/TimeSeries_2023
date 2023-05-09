% mixed normal garch model by chatGPT
% COMPLETELY WRONG ATM
% Mixed Normal GARCH(1,1) Parameters
omega = 0.1;  % Constant term
alpha = [0.1, 0.1];  % ARCH coefficients
beta = [0.2, 0.3];   % GARCH coefficients
T = 1000;     % Number of observations

% Preallocate arrays
epsilon = zeros(T, 1);   % Random shocks
sigma = zeros(T, 1);     % Conditional volatility

% Generate random shocks from mixed normal distribution
mu = [0.0, 0.0];       % Mean of the mixed normal distribution
sigma_mixed = [0.2, 0.5];  % Standard deviations of the mixed normal distribution
rho = 0.5;             % Correlation between the two components

rng(0); % Set the random seed for reproducibility
z1 = randn(T, 1);
z2 = randn(T, 1);
epsilon = sqrt(1-rho^2) * (sigma_mixed(1)*z1 + sigma_mixed(2)*z2) + rho * z1;

% Calculate conditional volatility
sigma(1) = sqrt(omega / (1 - sum(alpha) - sum(beta)));
for t = 2:T
    sigma(t) = sqrt(omega + sum(alpha .* epsilon(t-1).^2) + sum(beta .* sigma(t-1).^2));
end

% Plot Mixed Normal GARCH(1,1) process
figure;
plot(sigma, 'LineWidth', 1.5);
xlabel('Time');
ylabel('Conditional Volatility');
title('Simulated Mixed Normal GARCH(1,1) Process');
