function para = SDT_PKA(varargin)
%% 
% simulation of a Signal detection theory (SDT) based model to compute
% psychophysical kernel amplitude (PKA) in a 2AFC task
%
% INPUT:
% 'db' ... float; decision boundary (for Integration-to-Bound model)
% 'dt' ... float; contribution of decision time to confidence
% 'dc' ... float vector; range of signal strength. It can start from
% negative like [-5:5:50] to represent overlapped stimulus distributions
% 'race' ... race model (2 integrators). Set 'db' to be 2-element vector
% (e.g. [200 180])
% 'acceleration' ... float; acceleration parameter 
%  (negative: leaky-integration, 0: perfect integration, positive: attractor)
% 'link' ... link function for acceleration parameter: 'linear' or
% 'sigmoid'
% 'ntr' ... int; the number of trials (default is 10^6). The half is
%            automatically assigned as 0% signal trials. 
% 'nbin' ... int; the number of time bins to compute the time-resolved PKA 
% 'noise' ... float; pooling noise (internal noise). 22.8 is default. 
% 'cfnoise' ... float; noise on confidence judgement. 0 is default. 
% 'weights' ... vector with the same length of nframe 
% 'pkmethod' ... method to compute psychophysical kernel amplitude: 
%              0, weights as occurences of frames (Nienborg &
%              Cumming, 2009). In default.
%              1, image classification (stm for ch1 - stm for ch2)
%              2, logistic regression
% 'stmdist' ... type of stimulus distribution:  'uniform' or 'Gaussian'
% 'conftype' ... confidence type: 'sdt' (Hangya et al., 2016) or 'Bayes' (Adler & Ma, 2017)
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
% the stimulus noise (para.stmSD)
%
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++

% pre-set parameters
ntr = 10^5;
nframe = 20;

% evidence strength
dc = 0:5:50;

% race model
race_flag = 0;

% pooling noise (internal noise)
noise = 22.8;

% confidence noise
cfnoise = 0;

% weights on sensory time-course
weights = ones(1, nframe);

% decision boundry
db = inf; % this is essentially infinite

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

% resampling
repeat = 0;

% figure
plot_flag = 0;

j = 1;              
while  j <= length(varargin)
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
        case 'dc'
            dc = varargin{j+1};
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
        case 'repeat'
            repeat = varargin{j+1};
            j = j + 2;
        case 'plot'
            plot_flag = 1;
            j = j + 1;
    end
end

% stimulus distribution
stmMean = mean(dc);
stmSD = std(dc);
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

% trials x frames (dynamic stimuli)
stm = repmat(stm, 1, nframe);
stm = arrayfun(@(x) normrnd(x, 2*stmSD), stm);

% instantaneous noisy measurements (decision variable)
idv = arrayfun(@(x) normrnd(x, noise), stm);

% race model
dt = nframe*ones(ntr,1);
dbreach = zeros(ntr, 1);
if race_flag == 1 % 2 integrators
    if length(db)==1
        db = [db, db];
    end
    % instantaneous noisy measurements    
    idv1 = stm; idv2 = stm;
    sigma = [1 -sqrt(0.5); -sqrt(0.5) 1]; % van den Berg et al., 2016
    parfor n = 1:ntr
        for f = 1:nframe
            R = stm(n,f) + mvnrnd((1+noise/100)*[abs(stm(n,f)), -abs(stm(n,f))], sigma, 1);
            idv1(n,f) = R(1);
            idv2(n,f) = R(2);
        end
    end
    
    % sensory weighting
    idv1 = idv1.*repmat(weights, ntr, 1);
    idv2 = idv2.*repmat(weights, ntr, 1);
    
    % evidence integration
    dv1 = idv1; dv2 = idv2;
    for f = 2:nframe
        dv1(:,f) = dv1(:,f-1) + acceleration*linkf(dv1(:,f-1), link) + dv1(:,f);
        dv2(:,f) = dv2(:,f-1) + acceleration*linkf(dv2(:,f-1), link) + dv2(:,f);
    end

    % integration-to-bound and decision time
    for n = 1:ntr   
        idx1 = find(abs(dv1(n,:)) >= db(1), 1, 'first');
        idx2 = find(abs(dv2(n,:)) >= db(2), 1, 'first');
        if isempty(idx1)
            idx1 = nframe;
        end
        if isempty(idx2)
            idx2 = nframe;
        end       
        idx = min([idx1, idx2]);
        if idx1 > idx2
            dv1(n, idx:end) = dv1(n, idx);
            dv2(n, idx:end) = sign(dv2(n, idx))*db(2);
        elseif idx1 < idx2
            dv1(n, idx:end) = sign(dv1(n, idx))*db(1);
            dv2(n, idx:end) = dv2(n, idx);
        elseif idx1==idx2
            dv1(n, idx:end) = sign(dv1(n, idx))*db(1);
            dv2(n, idx:end) = sign(dv2(n, idx))*db(2);
        end
        dt(n) = idx;
        dbreach(n) = 1;
    end
    dv = dv1 + dv2;
    idv = [idv1, idv2];
else % one integrator
    % sensory weighting
    idv = idv.*repmat(weights, ntr, 1);

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
            dv(n, idx:end) = sign(dv(n, idx))*db;
            dt(n) = idx;
            dbreach(n) = 1;
        end
    end
end
nreach0 = 100*sum(dbreach==1)/ntr;
disp(['The % trials reaching the DB: ' num2str(nreach0)])

% noise levels (experiment & internal)
noisestm = std(stm(:));
noiseidv = std(idv(:));

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
[cf0, cf1] = median_split(conf);
% med = median(conf);
nreach1 = 100*sum(dbreach==1 & ismember(conf, cf0))/sum(ismember(conf, cf0));
nreach2 = 100*sum(dbreach==1 & ismember(conf, cf1))/sum(ismember(conf, cf1));
disp([num2str(nreach2) ...
    '% trials reached DB in high confidence, '...
    num2str(nreach1)...
    '% trials reached DB in low confidence '])
ch(ch==-1) = 0;
disp([num2str(sum(ch==0)) ' near-choices, ' num2str(sum(ch==1)) ' far-choices'])
disp('----------------------------------------------------------------')

% psychophysical kernel amplitude (PKA)
[pka_all, pka_hc, pka_lc] = getPKA(ss, stm, ch, conf, nbin, pkmethod, repeat);
if mean(isnan(pka_hc))
    pka_hc = 2*pka_all - pka_lc;
elseif mean(isnan(pka_lc))
    pka_lc = 2*pka_all - pka_hc;
end

% output argumant
para = struct('category', C, 'assigned_stm', ss, 'stm', stm, 'choice', ch, ...
    'accuracy', acc, 'confidence', conf, 'decisiontime',dt,'pka_method', pkmethod, 'pka', pka_all(1,:), ...
    'pka_highconf', pka_hc(1,:), 'pka_lowconf', pka_lc(1,:),...
    'choice_bias', sum(ch==0)/sum(ch==1), ...
    'nreach',nreach0,'nreach_highconf',nreach2,'nreach_lowconf',nreach1,...
    'noisestm',noisestm,'noiseidv',noiseidv, 'stmMean', stmMean, 'stmSD', stmSD);
if race_flag == 1
    para.dv = {dv1, dv2, dv};
else
    para.dv = dv;
end
if repeat > 0
    para.pka_error = pka_all(2,:);
    para.pka_highconf_error = pka_hc(2,:);
    para.pka_lowconf_error = pka_lc(2,:);
end

% visualization
if plot_flag==1
    % yellow and green
    y = [0.9576    0.7285    0.2285];
    g = [0.1059    0.4706    0.2157];
    
    close all;
    h = figure;
    subplot(1,2,1)
    nom = mean(pka_all(1,:));
    if repeat > 0
        err = pka_all(2,:); errl = pka_lc(2,:); errh = pka_hc(2,:); 
        errorbar(1:nbin, pka_all(1,:)/nom, err/nom, '-', 'color', [0 0 0], 'linewidth', 2)
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
    nom = max(pka_all(1,:));
    l = zeros(1,2);
    if repeat > 0
        l(2) = errorbar(1:nbin, pka_lc(1,:)/nom, errl/nom, '-', 'color', g, 'linewidth', 2);
        hold on;     
        l(1) = errorbar(1:nbin, pka_hc(1,:)/nom, errh/nom, '-', 'color', y, 'linewidth', 2);
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

%subfunctions
function y = linkf(x, link)
% link function for acceleration parameter
switch link
    case 'linear'
        y = x;
    case 'sigmoid'
        y = tanh(x);
end

function [cf0, cf1] = median_split(cf)
% median split
[~,idx] = sort(cf);
n = floor(length(cf)/2);
cf0 = cf(idx(1:n));
cf1 = cf(idx(n+1:end));