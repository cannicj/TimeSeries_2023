pAIC100 = pvecAIC10000(:,1);
pAIC1000 = pvecAIC10000(:,2);
pBIC100 = pvecBIC10000(:,1);
pBIC1000 = pvecBIC10000(:,2);
qAIC100 = qvecAIC10000(:,1);
qAIC1000 = qvecAIC10000(:,2);
qBIC100 = qvecBIC10000(:,1);
qBIC1000 = qvecBIC10000(:,2);
%%
est_matrix_p = [pvecAIC10000 pvecBIC10000]; % real model order is 3
est_matrix_q = [qvecAIC10000 qvecBIC10000]; % real model order is 2
occurances_p = zeros(4,10);
occurances_q = zeros(4,10);
for i=1:4
    eval_p = est_matrix_p(:,i);
    eval_q = est_matrix_q(:,i);
    for j=1:10
        occurances_p(i,j) = numel(find(eval_p==j));
        occurances_q(i,j) = numel(find(eval_q==j));
    end
end
%%
occurances_p/10000
occurances_q/10000
%%

tabulate(est_matrix_p(:,1))
tabulate(est_matrix_p(:,2))
tabulate(est_matrix_p(:,3))
tabulate(est_matrix_p(:,4))


%% Histograms for p
figure;
subplot(2,2,1);
histogram(pAIC100, 'BinMethod', 'integers','Normalization', 'probability');
hold on;
title('AIC: T=100');
xlabel('p');
ylabel('Probability');

subplot(2,2,2);
histogram(pAIC1000, 'BinMethod', 'integers','Normalization', 'probability');
hold on;
title('AIC: T=1000');
xlabel('p');
ylabel('Probability');

subplot(2,2,3);
histogram(pBIC100, 'BinMethod', 'integers','Normalization', 'probability');
hold on;
title('BIC: T=100');
xlabel('p');
ylabel('Probability');

subplot(2,2,4);
histogram(pBIC1000, 'BinMethod', 'integers','Normalization', 'probability');
hold on;
title('BIC: T=1000');
xlabel('p');
ylabel('Probability');

%% Histograms
figure;
subplot(2,2,1);
histogram(qAIC100, 'BinMethod', 'integers','Normalization', 'probability');
%h.Parent.XTick = h.BinEdges;
hold on;
title('AIC: T=100');
xlabel('q');
ylabel('Probability');

subplot(2,2,2);
histogram(qAIC1000, 'BinMethod', 'integers','Normalization', 'probability');
hold on;
title('AIC: T=1000');
xlabel('q');
ylabel('Probability');

subplot(2,2,3);
histogram(qBIC100, 'BinMethod', 'integers','Normalization', 'probability');
hold on;
title('BIC: T=100');
xlabel('q');
ylabel('Probability');

subplot(2,2,4);
histogram(qBIC1000, 'BinMethod', 'integers','Normalization', 'probability');
hold on;
title('BIC: T=1000');
xlabel('q');
ylabel('Probability');

