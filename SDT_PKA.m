function para = SDT_PKA(varargin)
%% 
% simulation of Signal detection theory (SDT) based models to compute
% psychophysical kernel amplitude (PKA) in a 2AFC task
%
% INPUT:
% 'db' ... float; decision boundary (for Integration-to-Bound model)
% 'dt' ... float; contribution of decision time to confidence
% 'race' ... race model (2 integrators)
% 'acceleration' ... float; acceleration parameter 
%  (negative: leaky-integration, 0: perfect integration, positive: attractor)
% 'link' ... link function for acceleration parameter: 'linear' or
% 'sigmoid'
% 'ntr' ... int; the number of trials (default is 10^6). The half is
%            automatically assigned as 0% signal trials. 
% 'nbin' ... int; the number of time bins to compute the time-resolved PKA 
% 'noise' ... float; pooling noise (internal noise). 22.8 is default. 
% 'cfnoise' ... float; noise on confidence judgement. 0 is default. 
% 'weights' ... vector with the same length of nframe (20 in default)
% 'pkmethod' ... method to compute psychophysical kernel amplitude: 
%              0, weights as occurences of frames (Nienborg &
%              Cumming, 2009). In default.
%              1, image classification (stm for ch1 - stm for ch2)
%              2, logistic regression
% 'stmdist' ... type of stimulus distribution:  'uniform' or 'Gaussian'
% 'conftype' ... confidence type: 'sdt' (Hangya et al., 2016) or 'Bayes' (Adler & Ma, 2017)
% 'overlap'... overlap (0 or >1) between P(s|C=-1) and P(s|C=1). Default is 0.
% 'repeat' ... the number of repeats for resampling (bootstrap)
% 'plot' ... plot PKA
%
% OUTPUT:
% matlab structure including relevant information
%
% EXAMPLE:
% para = SDT_PKA('db',140, 'plot')
%
% NOTE:
% Unit of decision bound here is arbitrary. Formally it needs to be normalized (divided) by
% the standard deviation of instantaneous decision variables (para.noiseidv)
%
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++

% pre-set parameters
ntr = 10^5;
nframe = 20;

% race model
race_flag = 0;

% pooling noise (internal noise)
noise = 22.8;

% confidence noise
cfnoise = 0;

% weights on sensory time-course
weights = ones(1, nframe);

% decision boundry
db = 10000; % this is essentially infinite

% taking into acount decision time
dtw = 0;

% acceleration parameter
acceleration = 0;

% link function for acceleration parameter
link = 'linear';

% number of bin
nbin = 4;

% way to compute PKA
pkmethod = 0;

% stimulus distribution
stmdist = 'uniform';

% confidence
conftype = 'sdt';

% overlap between stimulus distributions
overlap = 0;

% resampling
repeat = 0;

% figure
plot_flag = 0;

j = 1;              
while  j<= length(varargin)
    switch varargin{j}
        case 'ntr'
            ntr = varargin{j+1};
            j = j + 2;
        case 'noise'
            noise = varargin{j+1};
            j = j + 2;
        case 'cfnoise'
            cfnoise = varargin{j+1};
            j = j + 2;
        case 'race'
            race_flag = 1;
            j = j + 1;
        case 'weights' 
            weights = varargin{j+1};
            j = j + 2;
        case 'db'
            db = varargin{j+1};
            j = j + 2;
        case 'acceleration'
            acceleration = varargin{j+1};
            j = j + 2;
        case 'link'
            link = varargin{j+1};
            j = j + 2;
        case 'dt'
            dtw = varargin{j+1};
            j = j + 2;
        case 'nbin'
            nbin = varargin{j+1};            
             j = j + 2;
        case 'pkmethod'
            pkmethod = varargin{j+1};
            j = j + 2;
        case 'stmdist'
            stmdist = varargin{j+1};
            j = j + 2;
        case 'conftype'
            conftype = varargin{j+1};
            j = j + 2;
        case 'overlap'
            overlap = varargin{j+1};
            j = j + 2;
        case 'repeat'
            repeat = varargin{j+1};
            j = j + 2;
        case 'plot'
            plot_flag = 1;
            j = j + 1;
    end
end

% evidence strength
dc = -overlap:5:50;

% stimulus distribution
stmMean = 25;
stmSD = 15;
switch lower(stmdist)
    case 'uniform' 
        lendc = length(dc);
        pc1 = ones(1, lendc);
        pc2 = ones(1, lendc);
    case 'gaussian' 
        pc1 = gaussmf(dc,[stmSD stmMean]);
        pc2 = gaussmf(-dc,[stmSD -stmMean]);
end

% true categories
disp(['simulating ' num2str(ntr) ' trials...'])
C = datasample([-1 1], ntr)';

% signed stimulus
stm = C;
stm(C==1) = datasample(dc, sum(C==1), 'Weights', pc1);
stm(C==-1) = datasample(-dc, sum(C==-1), 'Weights', pc2);

% half of them are 0% signal trials
stm(randi(ntr, round(ntr/2), 1)) = 0;
ss = stm;

% trials x frames
stm = repmat(stm, 1, nframe);
stm = arrayfun(@(x) normrnd(x, 2*stmSD), stm);

% instantaneous noisy measurements (decision variable)
idv = arrayfun(@(x) normrnd(x, noise), stm);

% sensory weighting
idv = idv.*repmat(weights, ntr, 1);

% noise levels (experiment & internal)
noisestm = std(stm(:));
noiseidv = std(idv(:));

% race model
dt = nframe*ones(ntr,1);
dbreach = zeros(ntr, 1);
if race_flag == 1
    % two accumulators
    idv1 = idv; idv2 = idv; 
    idv1(idv1 < 0) = 0;
    idv2(idv2 > 0) = 0;

    % evidence integration
    dv1 = idv1; dv2 = idv2;
    for f = 2:nframe
        dv1(:,f) = dv1(:,f-1) + acceleration*linkf(dv1(:,f-1), link) + dv1(:,f);
        dv2(:,f) = dv2(:,f-1) + acceleration*linkf(dv2(:,f-1), link) + dv2(:,f);
    end

    % integration-to-bound and decision time
    for n = 1:ntr   
        idx1 = find(abs(dv1(n,:)) >= db, 1, 'first');
        idx2 = find(abs(dv2(n,:)) >= db, 1, 'first');
        if ~isempty(idx1) || ~isempty(idx2)
            idx = min([idx1, idx2]);
            dv1(n, idx:end) = dv1(n, idx);
            dv2(n, idx:end) = dv2(n, idx);
            dt(n) = idx;
            dbreach(n) = 1;
        end
    end
    dv = dv1 + dv2;
else % one integrator
    % evidence integration
    dv = idv;
    for f = 2:nframe
        dv(:,f) = dv(:,f-1) + acceleration*linkf(dv(:,f-1), link) + dv(:,f);
    end

    % integration-to-bound and decision time
    dt = nframe*ones(ntr,1);
    dbreach = zeros(ntr, 1);
    for n = 1:ntr   
        idx = find(abs(dv(n,:)) >= db, 1, 'first');
        if ~isempty(idx)
            dv(n, idx:end) = dv(n, idx);
            dt(n) = idx;
            dbreach(n) = 1;
        end
    end
end

nreach0 = 100*sum(dbreach==1)/ntr;
disp(['The % trials reaching the DB: ' num2str(nreach0)])

% confidence & choice
[conf, ch] = compute_confidence(dv(:,end), (nframe/2)*noise, conftype, stmdist,...
    dc, stmMean, stmSD);

% influence of decision time on confidence
if dtw > 0
    conf = 0.5 + (1/pi)*atan(2*(conf-0.5)./(dtw*dt/nframe));
end

% noise on confidence judgement
conf = conf + normrnd(0, cfnoise, size(conf));

% accuracy
acc = 1*(ch==C);

% median split
med = median(conf);
nreach1 = 100*sum(dbreach==1 & conf < med)/sum(conf < med);
nreach2 = 100*sum(dbreach==1 & conf > med)/sum(conf > med);
disp([num2str(nreach2) ...
    '% trials reached DB in high confidence, '...
    num2str(nreach1)...
    '% trials reached DB in low confidence '])
ch(ch==-1) = 0;
disp([num2str(sum(ch==0)) ' near-choices, ' num2str(sum(ch==1)) ' far-choices'])
disp('----------------------------------------------------------------')

% psychophysical kernel amplitude (PKA)
[pka_all, pka_hc, pka_lc] = getPKA(ss, stm, ch, conf, nbin, pkmethod);
if mean(isnan(pka_hc))
    pka_hc = 2*pka_all - pka_lc;
elseif mean(isnan(pka_lc))
    pka_lc = 2*pka_all - pka_hc;
end

% output argumant
para = struct('category', C, 'assigned_stm', ss, 'stm', stm, 'choice', ch, ...
    'accuracy', acc, 'confidence', conf, 'decisiontime',dt,'pka_method', pkmethod, 'pka', pka_all, ...
    'pka_highconf', pka_hc, 'pka_lowconf', pka_lc,...
    'choice_bias', sum(ch==0)/sum(ch==1), ...
    'nreach',nreach0,'nreach_highconf',nreach2,'nreach_lowconf',nreach1,...
    'noisestm',noisestm,'noiseidv',noiseidv);
if race_flag == 1
    para.dv = {dv1, dv2, dv};
else
    para.dv = dv;
end

% visualization
if plot_flag==1
    % yellow and green
    y = [0.9576    0.7285    0.2285];
    g = [0.1059    0.4706    0.2157];
    
    close all;
    h = figure;
    subplot(1,2,1)
    nom = mean(pka_all);
    if repeat > 0
        [err, errl, errh] = resamplePK(ss, stm, ch, conf, nbin, repeat, pkmethod);
        errorbar(1:nbin, pka_all/nom, err/nom, '-', 'color', [0 0 0], 'linewidth', 2)
    else
        plot(1:nbin, pka_all/nom, '-', 'color', [0 0 0], 'linewidth', 2)
    end    
    xlim([0.5 nbin + 0.5])
    xlabel('time bin')
    ylabel({'normalized', 'PKA'})
    set(gca, 'XTick', 1:nbin)
    set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
    
    subplot(1,2,2)
%     nom = mean([pka_hc, pka_lc]);
    nom = max(pka_all);
    l = zeros(1,2);
    if repeat > 0
        l(2) = errorbar(1:nbin, pka_lc/nom, errl/nom, '-', 'color', g, 'linewidth', 2);
        hold on;     
        l(1) = errorbar(1:nbin, pka_hc/nom, errh/nom, '-', 'color', y, 'linewidth', 2);
    else
        l(2) = plot(1:nbin, pka_lc/nom, '-', 'color', g, 'linewidth', 2);
        hold on;
        l(1) = plot(1:nbin, pka_hc/nom, '-', 'color', y, 'linewidth', 2);
    end    
    legend(l, 'high confidence', 'low confidence', 'location', 'best')
    legend('boxoff')
    xlim([0.5 nbin + 0.5])
    xlabel('time bin')
    set(gca, 'XTick', 1:nbin)
    set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
end
    
%%
%subfunctions
function y = linkf(x, link)
% link function for acceleration parameter
switch link
    case 'linear'
        y = x;
    case 'sigmoid'
        y = tanh(x);
end

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

function [pka_all, pka_hc, pka_lc] = getPKA(ss, stm, ch, cf, nbin, pkmethod)
% compute PKA in a specified method
switch pkmethod
    case 0 % Nienborg & Cumming, 2009
        [pka_all, pka_hc, pka_lc] = PKA_hn(ss, stm, ch, cf, nbin);
    case 1 % image classification
        [pka_all, pka_hc, pka_lc] = PKA_ic(ss, stm, ch, cf, nbin);
    case 2 % logistic regression
        [pka_all, pka_hc, pka_lc] = PKA_logreg(stm, ch, cf, nbin);
end

function pk = getKernel(stm, disval, ch)
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

function [pka_all, pka_hc, pka_lc] = PKA_hn(ss, stm, ch, cf, nbin)
% only 0% signal trials
stm = stm(ss==0,:);
ch = ch(ss==0);
cf = cf(ss==0);
% time-averaged PK in each time bin
stm = round(stm);
disval = unique(stm);
nd = length(disval);
med = median(cf);
tkernel = nan(nd, nbin);
tkernel_h = nan(nd, nbin);
tkernel_l = nan(nd, nbin);
begin = 1;
nframe = size(stm, 2);
frameperbin = floor(nframe/nbin);
for a = 1:nbin
    pk0 = getKernel(stm(:,begin:begin+frameperbin-1), disval, ch);
    pkh = getKernel(stm(cf > med, begin:begin+frameperbin-1), disval, ch(cf > med));
    pkl = getKernel(stm(cf < med, begin:begin+frameperbin-1), disval, ch(cf < med));
    tkernel(:,a) = pk0';
    tkernel_h(:,a) = pkh';
    tkernel_l(:,a) = pkl';
    begin = begin + frameperbin;
end
% PKA
pk = getKernel(stm, disval, ch);
pka_all = nan(1,nbin);
pka_hc = nan(1,nbin);
pka_lc = nan(1,nbin);
for a = 1:nbin
    pka_all(a) = dot(tkernel(:,a), pk);
    pka_hc(a) = dot(tkernel_h(:,a), pk);
    pka_lc(a) = dot(tkernel_l(:,a), pk);
end

function [pka_all, pka_hc, pka_lc] = PKA_ic(ss, stm, ch, cf, nbin)
% only 0% signal trials
stm = stm(ss==0,:);
ch = ch(ss==0);
cf = cf(ss==0);
% image classification to compute PKA
pka_all = mean(stm(ch==1,:), 1) - mean(stm(ch==0,:), 1);
med = median(cf);
pka_hc = mean(stm(ch==1 & cf > med,:), 1) - mean(stm(ch==0 & cf > med,:), 1);
pka_lc = mean(stm(ch==1 & cf < med,:), 1) - mean(stm(ch==0 & cf < med,:), 1);
pka_all = binstm(pka_all, nbin);
pka_hc = binstm(pka_hc, nbin);
pka_lc = binstm(pka_lc, nbin);

function [pka_all, pka_hc, pka_lc] = PKA_logreg(stm, ch, cf, nbin)
% logistic regression (intercept included to capture a bias)
stm = zscore(binstm(stm, nbin));
pka_all = glmfit(stm, ch, ...
    'binomial', 'link', 'logit', 'constant', 'on');
med = median(cf);
pka_hc = glmfit(stm(cf > med, :), ch(cf > med), ...
    'binomial', 'link', 'logit', 'constant', 'on');
pka_lc = glmfit(stm(cf < med, :), ch(cf < med), ...
    'binomial', 'link', 'logit', 'constant', 'on');
pka_all = pka_all(2:end); pka_hc = pka_hc(2:end); pka_lc = pka_lc(2:end);

function [err, err0, err1] = resamplePK(ss, stm, ch, cf, nbin, repeat, pkmethod)
% resampling procedure
ampr = nan(repeat, nbin);
ampr0 = nan(repeat, nbin);
ampr1 = nan(repeat, nbin);
ntr = length(ch);
for r = 1:repeat
    rtr = randi(ntr, ntr, 1);
    [ampr(r,:), ampr1(r,:), ampr0(r,:)] = getPKA(ss(rtr), stm(rtr,:), ch(rtr), cf(rtr), nbin, pkmethod);
end
err = std(ampr, [], 1);
err0 = std(ampr0, [], 1);
err1 = std(ampr1, [], 1);