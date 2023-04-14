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

function b = DurbinMA1959(y,q)
if length(y) <= q
    error("q has to be smaller than number of observations!")
end

% first, estimate a by OLS, see program 6.3
k = round(sqrt(length(y)));
[~, b] = yw(y,k); %get the a vector
a = [1; -b];

if q == 1
    % estimation of MA(1): equation 7
    num = a(1:(end-1))' * a(2:end);
    denom = sum(a.^2);
    b = - num / denom;

else
    LHS = zeros(q, q);
    RHS = zeros(q, 1);
    List = zeros(1,q);
    for j = 1:q
        List(1,j) = a(1:(end-j+1))' * a(j:end);
    end
    %Generate LHS and RHS
    for u=1:q
        LHS(u,:)=circshift(List,[0,u-1]);
        LHS = triu(LHS);
        LHS = LHS + tril(LHS',-1);
        RHS(u) = a(1:(end-u))' * a(u+1:end);
    end
    b = LHS \ (-RHS);
end
