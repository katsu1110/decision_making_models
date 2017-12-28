function [pka] = plot_EvidenceAccumulation(acc, e, ntime, ntrial, stim,model)
%%
% visualize behavior of the 'Evidence Accumulation' model given by Ralf
%
% simulate decision-variable for simulations of zero-signal case
% INPUT:
% - acc ... acceleration paramters: [0 0.05 -0.05] 
% - e ... average unsigned evidence per frame
% - ntime ... the number of frames: 100
% - ntrial ... the number of trials: 1000
% - stim ... 'binary' or 'normal'
% - model ... 'linear' or 'sigmoid'
% +++++++++++++++++++++++++++++++++++++++++++++++++

% preset parameters
if nargin<1, acc=[0 0.01 -0.01]; end % [paaerfect, accelerate, deccelerate]
if nargin<2, e=1; end % average unsigned evidence per frame
if nargin<3, ntime=100; end % 100 time steps
if nargin<4, ntrial=1000; end % number of simulated trials
if nargin<5, stim='normal'; end % stimulus to be simulated (see below)
if nargin<6, model='linear'; end % model to be simulated (see below)

% Ralf's function
%  [dv, stm] = Evidence_Accumulation(acc,ntime,ntrial,noise);
 [dv, stm] = Evidence_Accumulation2(acc,e,ntime,ntrial,stim,model);
stm = stm';

% visualization
nrow = size(dv, 2);
ncol = 5;
binsize = ntime;
% yellow and green
y = [0.9576    0.7285    0.2285];
g = [0.1059    0.4706    0.2157];
    
pka = nan(3*nrow, ntime);
c = 1;
close all;
h = figure;
rtr = randi(size(dv, 3), 100, 1);
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
    ch = getCh(dv, n);
    pka(c, :) = getPK(stm, ch, binsize);
    plot(pka(c,:), '-r')
    xlim([0.5 binsize+0.5])
    if n == 1
        ylabel('PKA')
    elseif n==nrow
        xlabel('time bin')
    end
%     subplot(1,ncol,5)
    subplot(nrow, ncol, n*5)
    idx_conf = conf_split(dv, n);
    pka1 = getPK(stm(idx_conf==0, :), ch(idx_conf==0), binsize);
    plot(pka1,'-','color',g)
    hold on;
    pka2 = getPK(stm(idx_conf==1, :), ch(idx_conf==1), binsize);
    plot(pka2,'-','color',y)
    xlim([0.5 binsize+0.5])
    pka(c+1, :) = pka1;
    pka(c+2, :) = pka2;
    c = c + 3;
    
    set(h, 'Name', ['Ralf-Evidence-Accumulation: acc = ' num2str(acc(n))], 'NumberTitle','off')
    
%     figure(200);
%     delta = (pka2 - pka1)/mean([pka2, pka1]);
%     scatter(acc(n), mean(delta(end - floor(ntime/4) + 1:end)), 100, 'filled', 'markerfacecolor','r','markerfacealpha',0.4)
%     hold on;
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

function pka = getPK(stm, ch, binsize)
frameperbin = floor(size(stm, 2)/binsize);
begin = 1;
binmat = nan(size(stm, 1), binsize);
for b = 1:binsize
    binmat(:,b) = mean(squeeze(stm(:, begin:begin+frameperbin-1)),2);
    begin = begin + frameperbin;
end
pka = mean(binmat(ch==1,:),1) - mean(binmat(ch==0,:), 1);

function [idx_conf] = conf_split(dv, acc)
conf = abs(squeeze(dv(end, acc, :)));
med = median(conf);
idx_conf = 1*(conf > med);
    
    