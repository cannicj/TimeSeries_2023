% Exercise 1) Part 3: Creating Confidence Intervals based on 
% the method of single parametric bootstrap
% initialize parameters
true_b = .2; sig=1; T=1000; sim=200;
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
% data for the table
disp(["Actual coverage matlab method:", num2str(mean(bool_matlab))])
disp(["Actual coverage MLE method:", num2str(mean(bool_MLE))])
disp(["Actual coverage Durbin method:", num2str(mean(bool_Durbin))])

disp(["CI Length matlab method:", num2str(mean(length_matlab))])
disp(["CI Length MLE method:", num2str(mean(length_MLE))])
disp(["CI Length Durbin method:", num2str(mean(length_Durbin))])
%%
% plot
plot(sort(length_matlab), 'r', 'LineWidth', 2);
hold on;
plot(sort(length_MLE), 'g', 'LineWidth', 2);
plot(sort(length_Durbin), 'b', 'LineWidth', 2);

% Add legend and title
legend('Matlab', 'MLE', 'Durbin');
title('CI lengths (sorted)');

%%
% Compute the density and mean for each data vector
[f1, xi1] = ksdensity(length_matlab);
mean1 = mean(length_matlab);
[f2, xi2] = ksdensity(length_MLE);
mean2 = mean(length_MLE);
[f3, xi3] = ksdensity(length_Durbin);
mean3 = mean(length_Durbin);

% Create the plot with the density lines and vertical mean lines
hold on;
plot(xi1, f1, 'r', 'LineWidth', 2);
line([mean1 mean1], [0 max(f1)], 'Color', 'r', 'LineStyle', '--');
plot(xi2, f2, 'g', 'LineWidth', 2);
line([mean2 mean2], [0 max(f2)], 'Color', 'g', 'LineStyle', '--');
plot(xi3, f3, 'b', 'LineWidth', 2);
line([mean3 mean3], [0 max(f3)], 'Color', 'b', 'LineStyle', '--');
hold off;

% Add labels and title
xlabel('Value');
ylabel('Density');
title('Density Plot of CI Lengths');

% Add legend
legend('Matlab','mean matlab', 'MLE','mean MLE', 'Durbin','mean Durbin');



