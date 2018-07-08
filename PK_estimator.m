function PK_estimator
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

% Nienborg & Cumming, 2007, 2009
[pka_est_hn, pk_est_hn] = PKA_hn(trst.C, trst.stm, beh.ch, nframe);

% just use logistic regression
beta = glmfit(zscore(trst.stm, 1), beh.ch, 'binomial', 'link', 'logit', 'constant', 'on');
pka_est_logreg = beta(2:end)';

% prior 
pka0 = mean(trst.stm(trst.C == 0 & beh.ch==1, :), 1) ...
    - mean(trst.stm(trst.C == 0 & beh.ch==0, :), 1);
pk0 = mean(trst.hdx_prob, 1).*sign(trst.hdx);

% parameter estimations
options = optimset('MaxFunEvals',10000,'maxiter',10000);
p0 = [pka0, pk0];
c = @(p)cost_params(p, trst, beh, options);
p = fminsearch(c, p0, options);
pka_est = p(1:nframe);
pk_est = p(nframe+1:end);

% visualize results
pk = {pk, pk_est_hn, pk_est};
pka = {pka, pka_est_hn, pka_est, pka0, pka_est_logreg};
labels = {'original', 'hn', 'me', 'ic', 'logreg'};
visualize(trst, beh, pk, pka, labels)

% subfunction
function trst = dynamic_stm_generator(ntrsig, nframe, hdx, sig, stmdist)
% generate dynamic stimulus sequence 

% deal with inputs
if nargin < 1; ntrsig = 100; end
if nargin < 2; nframe = 50; end
if nargin < 3; hdx = 0.3*[-1:0.25:1]; end
if nargin < 4
    sig = 0.5*[-1 -0.5 -0.25 -0.125 -0.0625 0 0.0625 0.125 0.25 0.5 1];
end
if nargin < 6; stmdist = 'gaussian'; end

% seed
rng(19891220);
% hdx (stimulus feature)
nhdx = length(hdx);
% percent signal
nsig = length(sig);
% % stimulus distribution
% stmMean = mean(sig);
% stmSD = std(sig);
% switch lower(stmdist)
%     case 'uniform' 
%         pc1 = ones(1, nsig);
%         pc2 = ones(1, nsig);
%     case 'gaussian' 
%         pc1 = gaussmf(dc,[stmSD stmMean]);
%         pc2 = gaussmf(-dc,[stmSD -stmMean]);
% end
% generate dynamic stimulus sequence
begin = 1;
pd = zeros(nsig, nhdx);
stm = zeros(ntrsig*nsig, nframe);
signal_hdx_idx = 2;
zeropos = find(hdx==0);
C = zeros(ntrsig*nsig, 1);
for n = 1:nsig
    % signal strength
    C(begin:begin+ntrsig-1) = sig(n);
    % probability distribution of hdx in each signal
    restper = 1 - abs(sig(n));
    pd(n,:) = (restper/nhdx)*ones(1,nhdx);
    sigpos = zeropos + sign(sig(n))*signal_hdx_idx;
    pd(n, sigpos) = abs(sig(n)) + pd(n, sigpos);
    % generate dynamic stimulus sequence in each signal
    stm_temp = datasample(hdx, ntrsig*nframe, 'Replace', true, 'Weights', pd(n,:));
    stm(begin:begin+ntrsig-1, :) = reshape(stm_temp, [ntrsig, nframe]);
    begin = begin + ntrsig;
end
% randomize order
rtr = randperm(size(stm,1));
C = C(rtr);
stm = stm(rtr, :);
% into the output
trst.hdx = hdx;
trst.sig = sig;
trst.hdx_prob = pd;
trst.C = C;
trst.stm = stm;

function beh = weight_on_stm(trst, pk, pka, threshold)
% assign time-weight (pka) and dv-weight (pk) on the stimulus sequence and
% determine the choice
% 
% PK (Psychophysical kernel) ... stm to decision variable (template)
% PKA (Psychophysical Kernel Amplitude) ... weight on time
% threshold ... psychophysical threshold to determine the noise level
%

% stimulus to decision variable
[ntr, nframe] = size(trst.stm);
idv = stm2dv(trst.stm, pk);

% noise on stm to dv
idv = idv + normrnd(0, threshold, ntr, nframe);

% weights on time
idv = idv.*repmat(pka, ntr, 1);

% integrate instantaneous decision variable
dv = cumsum(idv, 2);

% choice
ch = sign(dv(:,end));
ch(ch==0) = datasample([-1 1], sum(ch==0));
ch(ch==-1) = 0;

% accuracy
acc = datasample([0 1], ntr)';
acc(sign(ch) == sign(trst.C)) = 1;

% into output
beh.idv = idv;
beh.dv = dv;
beh.ch = ch;
beh.acc = acc;

function idv = stm2dv(stm, pk)
% assign PK weight
hdx = unique(stm);
lenhdx = length(hdx);
idv = stm;
for n = 1:lenhdx
    idv(stm == hdx(n)) = pk(n);
end

function c = cost_pk(pk, trst, beh, pka_est)
% simulation of choice behavior
beh_pred = weight_on_stm(trst, pk, pka_est, 0);
% choice prediction accuracy
predacc = sum(beh.ch == beh_pred.ch)/size(beh.ch, 1);
% print parameters
disp(['pred acc = ' num2str(100*predacc) ' %'])
% cost
c = 1 - predacc;

function c = cost_pka(p, trst, beh, pk)
% negative log likelihood of logistic regression
X = zscore(stm2dv(trst.stm, pk), 1);
X = [ones(size(X,1), 1) X];
p = [1 p];
c = sum(log(1 + exp(X*p'))) - (beh.ch)'*X*p';

function stmbin = binstm(stm, nbin)
% bin the stimulus metrix along with time
nframe = size(stm, 2);
begin = 1;
frameperbin = floor(nframe/nbin);
stmbin = nan(size(stm, 1), nbin);
for a = 1:nbin
    stmbin(:, a) = mean(stm(:,begin:begin+frameperbin-1), 2);
    begin = begin + frameperbin;
end

function [pk, pkvar] = getKernel(stm, disval, ch)
% trial averaged stimulus distributions split by choice
nd = length(disval);
ntr = length(ch);
svmat = zeros(ntr, nd);
for r = 1:ntr
    for d = 1:nd
        svmat(r,d) = sum(stm(r,:)==disval(d));
    end
end
% compute PK for 0% stimulus
pk = mean(svmat(ch==1,:)) - mean(svmat(ch==0,:));
pkvar = var(svmat(ch==1, :), [], 1) + var(svmat(ch==0, :), [], 1);

function [pka, pk] = PKA_hn(ss, stm, ch, nbin)
% only 0% signal trials
stm = stm(ss==0,:);
ch = ch(ss==0);
% time-averaged PK in each time bin
disval = unique(stm);
nd = length(disval);
% discretize stimuli, if too many unique values
if nd > 25
    [~,~,stm] = histcounts(stm, 11);
    stm = stm - mean(mean(stm));
end
disval = unique(stm);
nd = length(disval);
tkernel = nan(nd, nbin);
begin = 1;
nframe = size(stm, 2);
frameperbin = floor(nframe/nbin);
for a = 1:nbin
    pk0 = getKernel(stm(:,begin:begin+frameperbin-1), disval, ch);
    tkernel(:,a) = pk0';
    begin = begin + frameperbin;
end
% PKA
pk = mean(tkernel,2);
pka = nan(1,nbin);
for a = 1:nbin
    pka(a) = dot(tkernel(:,a), pk);
end

function c = cost_params(p, trst, beh, options)
% estimate pka by logistic regression
nframe = size(trst.stm, 2);
c_pka = @(x)cost_pka(x, trst, beh, p(nframe+1:end));
pka_est = fminsearch(c_pka, p(1:nframe), options);

% estimate pk by optimization procedure
c = cost_pk(p(nframe+1:end), trst, beh, pka_est);

function [x, y] = getPM(trst, beh)
% psychometric function
x = unique(trst.sig);
lenx = length(x);
y = zeros(1, lenx);
for i = 1:lenx
    y(i) = sum(beh.ch==1 & trst.C==x(i))/sum(trst.C==x(i));
end

function visualize(trst, beh, pk, pka, labels)
% compare results --- visualization
close all;
figure;

nframe = size(trst.stm, 2);

% hdx probability in each signal strength
subplot(2,5,1)
imagesc(trst.hdx, trst.sig, trst.hdx_prob)
xlabel('hdx')
ylabel('signal')
title('probability of hdx')
set(gca, 'box', 'off'); set(gca, 'TickDir', 'out');

% dynamic stimulus sequence
subplot(2,5,2)
imagesc(trst.stm)
xlabel('time (frames)')
ylabel('trials')
title('stimulus sequence')
set(gca, 'box', 'off'); set(gca, 'TickDir', 'out');

% instantaneous decision variable
subplot(2,5,3)
imagesc(beh.idv)
xlabel('time (frames)')
ylabel('trials')
title({'instantaneous', 'decision variable'})
set(gca, 'box', 'off'); set(gca, 'TickDir', 'out');

% dynamic stimulus sequence
subplot(2,5,4)
imagesc(beh.dv)
xlabel('time (frames)')
ylabel('trials')
title({'integrated', 'decision variable'})
set(gca, 'box', 'off'); set(gca, 'TickDir', 'out');

% psychometric function
subplot(2,5,5)
[x, y] = getPM(trst, beh);
plot(x,y,'-ok')
xlabel('% signal')
ylabel('P(ch = 1)')
title({'psychometric', 'function'})
set(gca, 'box', 'off'); set(gca, 'TickDir', 'out');

len_pk = length(pk);
len_pka = length(pka);
cols = hsv(max([len_pk, len_pka]));

% pk
for i = 1:length(pk)
    subplot(2,5,[6 7])
    hold on;
    plot(trst.hdx, pk{i}/(max(pk{i}) - min(pk{i})), '-', 'color', cols(i,:))
end
xlabel('hdx')
ylabel({'normalized','PK'})
set(gca, 'box', 'off'); set(gca, 'TickDir', 'out');

% pka
subplot(2,5,[9 10])
for i = 1:length(pka)
    hold on;
    plot(1:nframe, pka{i}/mean(pka{i}), '-', 'color', cols(i,:))
end
xx = [0.5 nframe+0.5];
xlim(xx)
yy = get(gca, 'YLim');

% labels
for i = 1:length(pka)
    text(xx(1)+0.1*(xx(2)-xx(1)), yy(1)+i*0.1*(yy(2)-yy(1)), ...
        labels{i}, 'color', cols(i,:))
end
xlabel('time bin')
ylabel({'normalized', 'PKA'})
set(gca, 'box', 'off'); set(gca, 'TickDir', 'out');