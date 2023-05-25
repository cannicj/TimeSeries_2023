% 1a) simulate the data according the a mixed-normal DGP

% parameters
T = 1000; % number of periods
k = 2; % number of components

% initialization of all parameters with the given paper page 230 (cite it!)
lambda = [0.82 0.18]; % weights from paper, change variable name! lamba
mu = [0.091 -0.415]; % means from paper
alpha_0 = [0.002 0.075]; % alpha_01, alpha_02, GARCH coefficients from paper
alpha_1 = [0.051 0.512]; % alpha_11, alpha_12, GARCH coefficients from paper
beta = diag([0.92 0.727]); % beta from paper (diag??)

% initialize zero vector for the simulation
sigma2 = zeros(T,k); % conditional variances
y = zeros(T,1); % return series
%%
% ensure zero mean for k = 2
mu(2) = -sum(lambda(1)/lambda(2)*mu(1));

% set initial values
sigma2(1, :) = alpha_0;
y(1) = sum(lambda.*alpha_0);

% simulation of the paths
for i = 2:T
    sigma2(i,:) = alpha_0 + alpha_1*(y(i-1)^2) + (beta * sigma2(i-1,:)')';
    y(i) = sum(lambda.*normrnd(mu, sqrt(sigma2(i,:))));
end
%%
figure;
plot(y, 'LineWidth', 1.5);
%xlabel('Time');
%ylabel('Conditional Volatility');
%title('Simulated GARCH(1,1) Process');

%% 
% parameter estimation with student-t-garch model
% mle extended for dof
% data vector is stored in y vector
babygarch(y)

%%
% 1c: one-step ahead forecast
% Now question 1c, we estimate the student-t several times, and record the
% VaR violations
violations = zero0;
for i = 250:(T-1)
    mle_params = babygarch(y(1:i));
    [zvec,sigvec] = ungarch(y(1:i), mle_params(1), mle_params(2), mle_params(3));
    VaR = tinv(0.01, mle_params(4))*sigvec(end);
    if VaR > y(i+1)
        violations = violations + 1;
    end
end
violations
