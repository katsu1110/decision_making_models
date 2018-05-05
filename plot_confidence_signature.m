function plot_confidence_signature(confidence, accuracy, stimulus)
% plot signatures of decision confidence:
% 1) accuracy vs confidence
% 2) psychometric functions split by confidence (median split)
% 3) confidence as a function of signal strength (correct and incorrect)
%
% each input should have trials x 1.

close all;
h = figure;
ntr = length(confidence);

% yellow and green
y = [0.9576    0.7285    0.2285];
g = [0.1059    0.4706    0.2157];

% 1) accuracy vs confidence
binsize = 20;
frameperbin = floor(ntr/binsize);
[sorted_cf, si] = sort(confidence);
sorted_acc = accuracy(si);
cfx = zeros(1, binsize);
acx = cfx;
begin = 1;
for b = 1:binsize
    cfx(b) = mean(sorted_cf(begin:begin+frameperbin-1));
    v = sorted_acc(begin:begin+frameperbin-1);
    acx(b) = sum(v==1)/length(v);
    begin = begin + frameperbin;
end
subplot(1,3,1)
plot(cfx, acx, '-k', 'linewidth', 2)
set(gca, 'box', 'off')
xlabel('confidence')
ylabel('accuracy (%)')
set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')

% 2) psychometric functions split by confidence
dc = unique(abs(stimulus));
lendc = length(dc);
pm0 = zeros(1, lendc);
pm1 = zeros(1, lendc);
cfc = zeros(1, lendc);
cfe = cfc;
for d = 1:lendc
    med = median(confidence(stimulus==dc(d)));
    acc0 = accuracy(stimulus==dc(d) & confidence < med);
    acc1 = accuracy(stimulus==dc(d) & confidence > med);
    pm0(d) = 100*sum(acc0==1)/length(acc0);
    pm1(d) = 100*sum(acc1==1)/length(acc1);
    cfc(d) = mean(confidence(accuracy==1 & stimulus==dc(d)));
    cfe(d) = mean(confidence(accuracy==0 & stimulus==dc(d)));
end
subplot(1,3,2)
plot(dc, pm1, '-', 'color', y, 'linewidth', 2)
hold on;
plot(dc, pm0, '-', 'color', g, 'linewidth', 2)
set(gca, 'box', 'off')
xlim([dc(1) dc(end)])
xlabel('signal strength')
ylabel('accuracy (%)')
legend('high', 'low', 'location', 'northwest')
legend('boxoff')
set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')

% 3) confidence as a function of signal strength (correct and incorrect)
subplot(1,3,3)
plot(dc, cfc, '-', 'color', zeros(1,3), 'linewidth', 2)
hold on;
plot(dc, cfe, '-', 'color', 0.4*ones(1,3), 'linewidth', 2)
set(gca, 'box', 'off')
xlim([dc(1) dc(end)])
xlabel('signal strength')
ylabel('confidence')
legend('correct', 'error', 'location', 'northwest')
legend('boxoff')
set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')

set(h, 'Name', 'Signatures of Decision Confidence', 'NumberTitle', 'off')