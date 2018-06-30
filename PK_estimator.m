function [pka_est, pk_est] = PK_estimator
% estimate pka (psychophysical kernel amplitude) & pk (psychophysical
% kernel)

% generate stimulus
ntrsig = 100;
nframe = 7;
trst = dynamic_stm_generator(ntrsig, nframe);

% decision variable
threshold = 1;
pk = [-0.1 -0.3 -0.5 -0.25 0.01 0.24 0.45 0.33 0.05];
pka = ones(1, nframe);
for n = 2:nframe
    pka(n) = pka(n-1)*0.8;
end
beh = weight_on_stm(trst, pk, pka, threshold);

% estimate pka by logistic regression
zerotr = trst.C == 0;
beta = glmfit(trst.stm(zerotr,:), beh.ch(zerotr), 'binomial', 'link', 'logit', 'constant', 'on');
pka_est = beta(2:end)';

% estimate pk by optimization procedure
pk0 = mean(trst.hdx_prob, 1).*sign(trst.hdx);
c = @(p)cost(p, trst, beh, pka_est);
options = optimset('MaxFunEvals',10000,'maxiter',10000);
pk_est = fminsearch(c, pk0, options);

% compare results
close all;
figure;

% pk
subplot(1,2,1)
plot(trst.hdx, pk/(max(pk) - min(pk)), '-k')
hold on; 
plot(trst.hdx, pk_est/(max(pk_est) - min(pk_est)), '-r')
xlabel('hdx')
ylabel({'normalized','PK'})
set(gca, 'box', 'off'); set(gca, 'TickDir', 'out');

% pka
subplot(1,2,2)
plot(1:nframe, pka/mean(pka), '-k')
hold on; 
plot(1:nframe, pka_est/mean(pka_est), '-r')
xlabel('time bin')
ylabel({'normalized', 'PKA'})
set(gca, 'box', 'off'); set(gca, 'TickDir', 'out');

% subfunction
function c = cost(p, trst, beh, pka_est)
% simulation of choice behavior
beh_pred = weight_on_stm(trst, p, pka_est, 0);
% choice prediction accuracy
predacc = sum(beh.ch == beh_pred.ch)/size(beh.ch, 1);
% print parameters
disp(['pk = ' num2str(p)])
disp(['pred acc = ' num2str(100*predacc) ' %'])
% cost
c = 1 - predacc;