% Simulation
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
c0_vec = zeros(1000,1);
c1_vec = zeros(1000,1);
d1_vec = zeros(1000,1);
dof_vec = zeros(1000,1);
violations_vec = zeros(1000,1);

for l=1:1000
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

%%

% Histogram settings
numBins = 20;   % Number of bins
barWidth = 0.8; % Manual bar width

% Plotting
figure;

% Histogram 1
subplot(2, 2, 1);
histogram(c0_vec, numBins, 'FaceColor', 'blue', 'EdgeColor', 'none', 'Normalization', 'pdf');
hold on;
x1 = linspace(min(c0_vec), max(c0_vec), 100);
density1 = ksdensity(c0_vec, x1);
plot(x1, density1, 'r', 'LineWidth', 2);
mean1 = mean(c0_vec);
line([mean1, mean1], ylim, 'Color', 'black', 'LineWidth', 2);
hold off;
xlim([0, 0.05]);
title('Estimates for c0');

% Histogram 2
subplot(2, 2, 2);
histogram(c1_vec, numBins, 'FaceColor', 'green', 'EdgeColor', 'none', 'Normalization', 'pdf');
hold on;
x2 = linspace(min(c1_vec), max(c1_vec), 100);
density2 = ksdensity(c1_vec, x2);
plot(x2, density2, 'r', 'LineWidth', 2);
mean2 = mean(c1_vec);
line([mean2, mean2], ylim, 'Color', 'black', 'LineWidth', 2);
hold off;
xlim([0, 0.2]);
title('Estimates for c1');

% Histogram 3
subplot(2, 2, 3);
histogram(d1_vec, numBins, 'FaceColor', 'red', 'EdgeColor', 'none', 'Normalization', 'pdf');
hold on;
x3 = linspace(min(d1_vec), max(d1_vec), 100);
density3 = ksdensity(d1_vec, x3);
plot(x3, density3, 'r', 'LineWidth', 2);
mean3 = mean(d1_vec);
line([mean3, mean3], ylim, 'Color', 'black', 'LineWidth', 2);
hold off;
xlim([0, 1]);
title('Estimates for d1');

% Histogram 4
subplot(2, 2, 4);
histogram(dof_vec, numBins, 'FaceColor', 'yellow', 'EdgeColor', 'none', 'Normalization', 'pdf');
hold on;
x4 = linspace(min(dof_vec), max(dof_vec), 100);
density4 = ksdensity(dof_vec, x4);
plot(x4, density4, 'r', 'LineWidth', 2);
mean4 = mean(dof_vec);
line([mean4, mean4], ylim, 'Color', 'black', 'LineWidth', 2);
hold off;
xlim([-0, 500]);
title('Estimates for dof');

%%

% Histogram 
figure;
histogram(violations_vec, 20, 'Normalization', 'pdf');
hold on;
x4 = linspace(min(violations_vec), max(violations_vec), 100);
density4 = ksdensity(violations_vec, x4);
plot(x4, density4, 'r', 'LineWidth', 2);
mean4 = mean(violations_vec);
line([mean4, mean4], ylim, 'Color', 'black', 'LineWidth', 2);
hold off;
xlim([0, 20]);
title('Violations');
