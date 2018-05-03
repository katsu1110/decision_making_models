function [pka, dv] = EvidenceAccumulation(acc, e, ntime, ntrial, stim, model, plot_flag, logreg_flag)
%%
% simulate & visualize behavior of the 'Evidence Accumulation' model given by Ralf
% Haefner @ Rochester Univ.
%
% simulate decision-variable for simulations of zero-signal case
% INPUT:
% - acc ... acceleration paramters: [0 0.05 -0.05] 
% - e ... average unsigned evidence per frame
% - ntime ... the number of frames: 100
% - ntrial ... the number of trials: 1000
% - stim ... 'binary' or 'normal'
% - model ... 'linear' or 'sigmoid'
% - plot_flag ... 0 (no figure) or 1 (plot)
% - logreg_flag ... 0 (image classification) or 1 (logistic regression)
% +++++++++++++++++++++++++++++++++++++++++++++++++

% preset parameters
if nargin<1, acc=[0 0.01 -0.01]; end % [perfect, accelerate, deccelerate]
if nargin<2, e=1; end % average unsigned evidence per frame
if nargin<3, ntime=100; end % 100 time steps
if nargin<4, ntrial=1000; end % number of simulated trials
if nargin<5, stim='normal'; end % stimulus to be simulated (see below)
if nargin<6, model='linear'; end % model to be simulated (see below)
if nargin<7, plot_flag=1; end % visualize or just store data
if nargin<8, logreg_flag=1; end % logistic regression to compute PKA

% Ralf's function (subfunction)
[dv, stm] = Simulate_Evidence_Accumulation(acc,e,ntime,ntrial,stim,model);
stm = stm';

nrow = size(dv, 2);
ncol = 5;
binsize = ntime;    
pka = nan(3*nrow, ntime);
c = 1;
%%
% store data
for n = 1:nrow
    ch = getCh(dv, n);
    pka(c, :) = getPK(stm, ch, binsize, logreg_flag);
    idx_conf = conf_split(dv, n);
    pka(c+1, :) = getPK(stm(idx_conf==0, :), ch(idx_conf==0), binsize, logreg_flag);
    pka(c+2, :) = getPK(stm(idx_conf==1, :), ch(idx_conf==1), binsize, logreg_flag);
    c = c + 3;
end

%%
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
    %     h = figure;
    %     subplot(1,ncol,1)
        subplot(nrow, ncol, n*5 -4)
        imagesc(squeeze(dv(:,n,:))')
        ylabel({['acc = ' num2str(acc(n))], 'trials'})
        if n == 1
            title('decision variable')
        elseif n==nrow
            xlabel('time')
        end
    %     subplot(1,ncol,2)
        subplot(nrow, ncol, n*5 -3)
        plot(squeeze(dv(:,n,rtr)))
    %     subplot(1,ncol,3)
        subplot(nrow, ncol, n*5 -2)
        histogram(squeeze(dv(end, n, :)))
    %     xlim([-20 20])
    %     histogram(squeeze(dv(end, n, :)), [-2000:10:2000])
    %     subplot(1,ncol,4)
        subplot(nrow, ncol, n*5 -1)
        plot(pka(c,:), '-r')
        xlim([0.5 binsize+0.5])
        if n == 1
            ylabel('PKA')
        elseif n==nrow
            xlabel('time bin')
        end
    %     subplot(1,ncol,5)
        subplot(nrow, ncol, n*5)
        plot(pka(c+1,:),'-','color',g)
        hold on;
        plot(pka(c+2,:),'-','color',y)
        xlim([0.5 binsize+0.5])
        c = c + 3;

        set(h, 'Name', ['Ralf-Evidence-Accumulation: acc = ' num2str(acc(n))], 'NumberTitle','off')

    %     figure(200);
    %     delta = (pka2 - pka1)/mean([pka2, pka1]);
    %     scatter(acc(n), mean(delta(end - floor(ntime/4) + 1:end)), 100, 'filled', 'markerfacecolor','r','markerfacealpha',0.4)
    %     hold on;
    end
end
% figure(200);
% xx = get(gca, 'XLim');
% plot(xx, [0 0], ':k')
% xlabel('accleralation')
% ylabel(['\Delta PKA (high - low conf) at time ' num2str(4)])


function ch = getCh(dv, acc)
ch = sign(squeeze(dv(end, acc, :)));
ch(ch==0) = datasample([-1 1], sum(ch==0));
ch(ch==-1) = 0;

function pka = getPK(stm, ch, binsize, logreg_flag)
frameperbin = floor(size(stm, 2)/binsize);
begin = 1;
binmat = nan(size(stm, 1), binsize);
for b = 1:binsize
    binmat(:,b) = mean(squeeze(stm(:, begin:begin+frameperbin-1)),2);
    begin = begin + frameperbin;
end
if logreg_flag==0
    % image classification
    pka = mean(binmat(ch==1,:),1) - mean(binmat(ch==0,:), 1);
else
    % logistic regression
%     pka = glmfit(binmat, ch, 'binomial', 'link', 'logit', 'constant', 'on');
%     pka = pka(2:end)';
    B = lassoglm(binmat, ch, 'binomial');
    pka = B(:, 1)';
end

function [idx_conf] = conf_split(dv, acc)
conf = abs(squeeze(dv(end, acc, :)));
med = median(conf);
idx_conf = 1*(conf > med);
    
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


    