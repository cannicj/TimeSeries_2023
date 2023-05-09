% Simulation of GARCH and respective parameter estimation
% own simple program that does this
% initialize parameters for a GARCH(1,1)
c0=0.05;  % Constant term
c1=0.2;  % ARCH coefficient
d1=0.7;  % GARCH coefficient
T=1000;  % Number of observations
% Preallocate arrays
epsilon = zeros(T, 1);   % Random shocks
sigma = zeros(T, 1);     % Conditional volatility
% Generate random shocks
sigma(1) = c0/(1-c1-d1); % we use vol, not variance as in the book
epsilon(1) = normrnd(0,sigma(1));
% simulate the rest of the process
for t = 2:T
    sigma(t) = c0 + c1 * epsilon(t-1) + d1 * sigma(t-1); % vol
    epsilon(t) = normrnd(0,sigma(t));
end
% Plot GARCH(1,1) process
figure;
plot(sigma, 'LineWidth', 1.5);
xlabel('Time');
ylabel('Conditional Volatility');
title('Simulated GARCH(1,1) Process');


% matlab versions using garch, simulate, estimate, forecast


