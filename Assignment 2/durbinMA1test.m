y = ma1sim(1000,1,0.5);

k = round(sqrt(length(y)));
[~, b] = yw(y,k); %get the a vector
a = [1; -b];
q=1;


% estimation of MA(1): equation 7
num = a(1:(end-1))' * a(2:end);
denom = sum(a.^2);
test1 = - num / denom


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
test2 = LHS \ (-RHS)
