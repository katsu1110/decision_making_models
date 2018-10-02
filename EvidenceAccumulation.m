function [pka, dv] = EvidenceAccumulation(acc, e, ntime, ntrial, stim, model, cfnoise, plot_flag)
%%
% simulate & visualize behavior of the 'Evidence Accumulation' model 
% given by Ralf Haefner @ Rochester Uni.
%
% simulate decision-variable for simulations of zero-signal case
% INPUT:
% - acc ... acceleration paramters: e.g. [-0.4 0 0.4] 
% - e ... average unsigned evidence per frame
% - ntime ... the number of frames: 100
% - ntrial ... the number of trials: 1000
% - stim ... 'binary' or 'normal'
% - model ... 'linear' or 'sigmoid'
% - cfnoise ... noise on confidence judgement
% - plot_flag ... 0 (no figure) or 1 (plot)
%
% OUTPUT: 
% pka (psychophysical kernel amplitude split by confidence)
% dv (decision variables)
%
% EXAMPLE: 
% [pka, dv] = EvidenceAccumulation([-0.4 0 0.4], 1, 4, 10000, 'normal', 'sigmoid', 1);
% +++++++++++++++++++++++++++++++++++++++++++++++++

% preset parameters
if nargin<1, acc=[0 0.01 -0.01]; end % [perfect, accelerate, deccelerate]
if nargin<2, e=1; end % average unsigned evidence per frame
if nargin<3, ntime=100; end % 100 time steps
if nargin<4, ntrial=1000; end % number of simulated trials
if nargin<5, stim='normal'; end % stimulus to be simulated (see below)
if nargin<6, model='linear'; end % model to be simulated (see below)
if nargin<7, cfnoise = 0; end % noise on confidence judgement
if nargin<8, plot_flag=1; end % visualize or just store data

% Ralf's function (in subfunction)
[dv, stm] = Simulate_Evidence_Accumulation(acc,e,ntime,ntrial,stim,model);
stm = stm';
nrow = size(dv, 2);
ncol = 5;
binsize = ntime;    
pka = nan(3*nrow, ntime);
c = 1;
for n = 1:nrow
    [ch, conf] = getCh(dv, n);
    conf = conf + normrnd(0, cfnoise*std(conf), size(conf));
    pka(c:c+2, :) = getPK(stm, ch, conf, binsize);
    c = c + 3;
end

% visualize
if plot_flag==1
    % yellow and green
    y = [0.9576    0.7285    0.2285];
    g = [0.1059    0.4706    0.2157];

    close all;
    h = figure;
    rtr = randi(size(dv, 3), 100, 1);
    c = 1;
    for n = 1:nrow
        % trial-by-trial decision variable
        subplot(nrow, ncol, n*5 -4)
        imagesc(squeeze(dv(:,n,:))')
        ylabel({['acc = ' num2str(acc(n))], 'trials'})
        if n == 1
            title('decision variable')
        elseif n==nrow
            xlabel('time')
        end
        % randomly selected traces of decision variables
        subplot(nrow, ncol, n*5 -3)
        plot(squeeze(dv(:,n,rtr)))
        xlabel('frames')
        ylabel('decision variable')
        subplot(nrow, ncol, n*5 -2)
        histogram(squeeze(dv(end, n, :)))
        xlabel('last decision variable')
        ylabel('trials')
        subplot(nrow, ncol, n*5 -1)
        plot(pka(c,:), '-r')
        xlim([0.5 binsize+0.5])
        ylabel('PKA')
        xlabel('frames')
        subplot(nrow, ncol, n*5)
        plot(pka(c+1,:),'-','color',g)
        hold on;
        plot(pka(c+2,:),'-','color',y)
        xlim([0.5 binsize+0.5])
        ylabel('PKA')
        xlabel('frames')
        c = c + 3;
    end
    set(h, 'Name', ['Ralf-Evidence-Accumulation: acc = ' num2str(acc(n))], 'NumberTitle','off')
end

% subfunctions
function [ch, conf] = getCh(dv, acc)
% trial-by-trial choice & confidence
ch = sign(squeeze(dv(end, acc, :)));
ch(ch==0) = datasample([-1 1], sum(ch==0));
ch(ch==-1) = 0;
conf = abs(squeeze(dv(end, acc, :)));

function pka = getPK(stm, ch, cf, binsize)
% get psychophysical kernel amplitude
frameperbin = floor(size(stm, 2)/binsize);
begin = 1;
binmat = nan(size(stm, 1), binsize);
med = median(cf);
for b = 1:binsize
    binmat(:,b) = mean(squeeze(stm(:, begin:begin+frameperbin-1)),2);
    begin = begin + frameperbin;
end
pka = nan(3, binsize);
pka(1,:) = mean(binmat(ch==1,:),1) - mean(binmat(ch==0,:), 1);
pka(2,:) = mean(binmat(ch==1 & cf < med,:),1) - mean(binmat(ch==0 & cf < med,:), 1);
pka(3,:) = mean(binmat(ch==1 & cf > med,:),1) - mean(binmat(ch==0 & cf > med,:), 1);
    
function [dv, evidence] = Simulate_Evidence_Accumulation(acc,e,ntime,ntrial,stim,model)
% returns decision-variable for simulations of zero-signal case
% dv[number of time steps, number of acceleration parameters, number of
% trials]
% Allows one to obtain reaction time and choice for a particular bound B as
%      reaction time for first trial for first acc parameter: 
%               rt = find(dv(:,idx_acc,idx_trial)>2,1)
%               if isempty(rt), rt=size(dv,1); end
%      corresponding choice: 
%               choice = sign(dv(rt,idx_acc,idx_trial));

dv=zeros(ntime, length(acc),ntrial); % decision variable
switch stim
  case 'normal'
    evidence=randn(ntime,ntrial);
  case 'binary'
    evidence=-0.5+[ones(ntime/2,ntrial); zeros(ntime/2,ntrial)];
    for i=1:ntrial
      evidence(:,i)=evidence(randperm(ntime),i);
    end
  otherwise
    error(['invaid stim: ' stim]);
end
evidence=e*evidence;
evidence_upscale=ones(length(acc),1);
acc=acc'*ones(1,ntrial);
dv(1,:,:)=evidence_upscale*evidence(1,:); % starting point
for i=2:ntime
  x=reshape(dv(i-1,:,:),size(acc,1),ntrial);
  switch model
    case 'linear'
      x=x+acc.*x;       % dv=dv+acc*dv
    case 'sigmoid'
      x=x+acc.*tanh(x); % dv=dv+acc*tanh(dv)
    otherwise
      error(['invalid model: ' model]);
  end
  dv(i,:,:)=x+evidence_upscale*evidence(i,:);
end