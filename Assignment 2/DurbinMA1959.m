% Estimation method from the Paper from Durbin
% what to do description from Marc P.
% program his method of computing the estimator of the MA(1) parameter
% (equation 7), and the MA(q) parameters (equation 15).
% "a" terms are the AR estimated parameters, obtained from OLS: chapter 4
% For the MA(q), compute program for arbitrary q, provided q < T (check
% this, T = number of observations in time series)



% I suggest you make ONE computer program (DurbinMA1959)
% (input: a vector time series, and a value of q). If q=1, then you 
% have an IF THEN ELSE, and use his eq 7 (special case of eq 15).
% run the codes for the q>1 case but in fact for q=1, you get exactly 
% same output.

% Now for your simulations: For the MA(1) model, use T=100 observations, 
% and a grid of true b values: -0.9, -0.8, ... 0, 0.1, up to 0.9. Use the 
% Durbin method in conjunction with simulation of the MA(1) process (it is
% easy to simulate an MA model, and I have codes in my book) to make nice 
% performance *graphics* (not tables).

function [param, stderr, resid] = DurbinMA1959(y,q)
if length(y) <= q
    disp("q has to be smaller than number of observations!")
    return
end
k = 10;
if q = 1
    % estimation of MA(1): equation 7
    % first, estimate a by OLS, see chapter 4.3.1
    y_t = y(2:end); % vector of length T
    y_tminus1 = y(1:end-1); % vector of length T
    a = (y_tminus1' * y_t)/(y_tminus1' * y_tminus1); %get a
    for i=1:(length(a)-1)
        b_num(i) = a(i)*a(i+1);
    end
    b_num_sum = sum(b_num);
    b_denom = sum(a^2);
    b_estimate = - b_num_sum/b_denom;
else
    % estimation of MA(q): equation 15
    
end
