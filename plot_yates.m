function plot_yates(res)
% plot simulation results from 'Yate_simulation.m'
% INPUT: result structure from 'Yate_simulation.m'

close all;

% stm/co kernels
figure(1);
subplot(3,4,1)
t = [1:length(res.kernel_stm)] - 1;
plot(t, res.kernel_stm, '-b', 'linewidth',2)
hold on;
plot(t, -res.kernel_stm, '-r', 'linewidth',2)
xlim([t(1) t(end)])
xlabel('time (ms)')
title('stimulus kernel')
set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')

subplot(3,4,2)
plot(t, res.kernel_co, '-k', 'linewidth',2)
xlim([t(1) t(end)])
xlabel('time (ms)')
title('contrast kernel')
set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')

subplot(3,4,3)
if isfield(res, 'kernel_ht')
    plot(res.kernel_ht, '-k', 'linewidth',2)
    set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
end
xlabel('time (ms)')
title('history kernel')    

% example PSTH
subplot(3,4,4)
nneuron = 1;
offset = 100;
time = [1:length(res.neuron(1).spk1(1,:))] - offset;
imagesc(time, 1:2*nneuron, [res.neuron(1).spk1(1,:); res.neuron(1).spk2(1,:)])
colormap(gca, flipud(bone))
hold on;
plot(time, 1.5*ones(1, length(time)), '-m')
hold on;
plot(time, 1.5+ 3*res.stm(1,:), '-g', 'linewidth',2)
hold on;
xlim([time(1) time(end)])
ylim([0 2*nneuron+1])   
xlabel('time (ms)')
ylabel('neurons')
title('example: trial 1')
set(gca,'box','off'); set(gca, 'TickDir', 'out')

% overall stimulus sign
sumstm = sum(res.stm,2);
stmsign = sign(sumstm);
med_p = median(sumstm(stmsign > 0));
med_n = median(sumstm(stmsign < 0));
stmsign_p2 = stmsign > 0 & sumstm > med_p;
stmsign_n2 = stmsign < 0 & sumstm < med_n;

lent = length(time);
ras1 = zeros(nneuron, lent);
ras2 = zeros(nneuron, lent);
spkch1 = zeros(1, lent);
spkch2 = zeros(1, lent);
spkstm1 = zeros(1, lent);
spkstm2 = zeros(1, lent);
for n = 1:nneuron
    ras1(n,:) = 1*(res.neuron(n).spk1(1,:) > 0);
    ras2(n,:) = 1*(res.neuron(n).spk2(1,:) > 0);
    spkch1 = spkch1 + mean(res.neuron(n).spk1(res.ch==1, 1:lent), 1);
    spkch2 = spkch2 + mean(res.neuron(n).spk1(res.ch==0, 1:lent), 1);
    spkstm1 = spkstm1 + mean(res.neuron(n).spk1(stmsign_p2, 1:lent), 1);
    spkstm2 = spkstm2 + mean(res.neuron(n).spk1(stmsign_n2, 1:lent), 1);
end
len_frame = length(res.neuron(1).spk1(1,:)) - 2*offset;
nbin = length(res.pka);
frameperbin = len_frame/nbin;

% PSTH (preferred or null stimulus)
subplot(3,4,5)
plot(time, spkstm1/nneuron, '-b', 'linewidth',1)
hold on;
plot(time, spkstm2/nneuron, '-r', 'linewidth',1)
hold on;
yy = get(gca, 'YLim');
begin = 1;
for i = 1:nbin+1
    hold on;
    plot(begin*[1 1],yy, ':k')
    begin = begin + frameperbin;
end
xlim([-offset len_frame+offset])
ylim(yy)    
xlabel('time (ms)')
ylabel('psth (by stm)')
set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')

% PSTH (preferred or null choice)
subplot(3,4,6)
plot(time, spkch1/nneuron, '-b', 'linewidth',1)
hold on;
plot(time, spkch2/nneuron, '-r', 'linewidth',1)
hold on;
yy = get(gca, 'YLim');
begin = 1;
for i = 1:nbin+1
    hold on;
    plot(begin*[1 1],yy, ':k')
    begin = begin + frameperbin;
end
xlim([-offset len_frame+offset])
ylim(yy)    
xlabel('time (ms)')
ylabel('psth (by choice)')
set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')

% PTA
subplot(3,4,7);
col = copper(nbin);
PTA = cell(1, nbin);
for n = 1:nneuron
    PTA1 = PTA_easy(res.neuron(n).spk1, res.stm, max(res.stm(:)),nbin, offset, frameperbin);
    PTA2 = PTA_easy(res.neuron(n).spk2, res.stm, min(res.stm(:)),nbin, offset, frameperbin);
    for b = 1:nbin
        if n==1
            PTA{b} = (PTA1{b} + PTA2{b})/2;
        else
            PTA{b} = PTA{b} + (PTA1{b} + PTA2{b})/2;
        end
    end
end
for n = 1:nbin
    plot(1:len_frame, PTA{n}/nneuron,...
        'color',col(n,:),'linewidth',2)
    hold on;
end
xlabel('time (ms)')
ylabel('PTA')
set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')

% confidence
subplot(3,4,8)
histogram(res.confidence)
hold on;
me = median(res.confidence);
yy = get(gca, 'YLim');
plot(me*[1 1], yy, '-r')
text(me*1.1, yy(1)+(yy(2)-yy(1))*0.8, ['median = ' num2str(me)], 'color', 'r')
xlabel('confidence')
ylabel('trials')
set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')

% PKA
subplot(3,4,9)
nom = mean(res.pka);
if isfield(res, 'pka_err')
    errorbar(1:nbin, res.pka/nom, res.pka_err/nom, '-k', ...
        'linewidth', 2, 'CapSize',0)
else
    plot(1:nbin, res.pka/nom, '-k', 'linewidth', 2)
end
xlim([0.75 nbin+0.25])
xlabel('time bin')
ylabel('PKA')
set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')

% yellow and green
y = [0.9576    0.7285    0.2285];
g = [0.1059    0.4706    0.2157];

subplot(3,4,[10 11])
nom = max(res.pka);
l = nan(1,2);
if isfield(res, 'pka_err')
    l(1) = errorbar(1:nbin, res.pka_highconf/nom,res.pka_highconf_err/nom, ...
        '-', 'color', y, 'linewidth', 2, 'CapSize', 0);
    hold on;
    l(2) = errorbar(1:nbin, res.pka_lowconf/nom, res.pka_lowconf_err/nom, ...
        '-', 'color', g, 'linewidth', 2, 'CapSize', 0);
else
    l(1) = plot(1:nbin, res.pka_highconf/nom, '-', 'color', y, 'linewidth', 2);
    hold on;
    l(2) = plot(1:nbin, res.pka_lowconf/nom, '-', 'color', g, 'linewidth', 2);
end
xlim([0.75 nbin+0.25])
xlabel('time bin')
ylabel('PKA')
set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
legend(l, 'high confidence', 'low confidence', 'location', 'eastoutside')
legend('boxoff')

function PTA = PTA_easy(spk, stm, pref, nbin, offset, frameperbin)
% compute the pulse triggered averaging quickly
PTA = cell(1, nbin);
begin = 1 + offset;
for n = 1:nbin
    PTA{n} = mean(spk(stm(:,begin)==pref,1+offset:nbin*frameperbin+offset),1)...
        - mean(spk(:,1+offset:nbin*frameperbin+offset),1);
    begin = begin + frameperbin;
end