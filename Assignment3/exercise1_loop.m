% loop exercise 1
%%
% parameters
T = 1000; % number of periods
k = 2; % number of components

% initialization of all parameters with the given paper page 230 (cite it!)
lambda = [0.82 0.18]; % weights from paper, change variable name! lamba
mu = [0.091 -0.415]; % means from paper
alpha_0 = [0.002 0.075]; % alpha_01, alpha_02, GARCH coefficients from paper
alpha_1 = [0.051 0.512]; % alpha_11, alpha_12, GARCH coefficients from paper
beta = diag([0.92 0.727]); % beta from paper
% ensure zero mean for k = 2
mu(2) = -sum(lambda(1)/lambda(2)*mu(1));
%%
%initialize storage vectors for the parameters of interest
c0_vec = zeros(10,1);
c1_vec = zeros(10,1);
d1_vec = zeros(10,1);
dof_vec = zeros(10,1);
violations_vec = zeros(10,1);

for l=1:10
    % initialize zero vector for the simulation
    sigma2 = zeros(T,k); % conditional variances
    y = zeros(T,1); % return series
    % set initial values
    sigma2(1, :) = alpha_0;
    y(1) = sum(lambda.*alpha_0);
    % simulation
    for i = 2:T
        sigma2(i,:) = alpha_0 + alpha_1*(y(i-1)^2) + (beta * sigma2(i-1,:)')';
        y(i) = sum(lambda.*normrnd(mu, sqrt(sigma2(i,:))));
    end
    params = babygarch(y);
    c0_vec(l) = params(1); c1_vec(l) = params(2); d1_vec(l) = params(3); dof_vec(l) = params(4);
    for i = 250:(T-1)
        mle_params = babygarch(y(1:i));
        [zvec,sigvec] = ungarch(y(1:i), mle_params(1), mle_params(2), mle_params(3));
        VaR = tinv(0.01, mle_params(4))*sigvec(end);
        if VaR > y(i+1)
            violations_vec(l) = violations_vec(l) + 1;
        end
    end
end

