function plot_yates(res_yates)
% plot simulation results from 'Yate_simulation.m'

close all;
h = figure;

% debug stimulus statistics
subplot(2,3,1)
histogram(sumstm)
xlabel('actual signal strength')
ylabel('trials')

subplot(2,3,4)
plot(1:lenv, mean(stm,1))
hold on;
plot(1:lenv, median(stm,1))
xlim([1 lenv])
xlabel('time during stimulus presentation')
ylabel('signal')
legend('mean', 'median')
legend('boxoff')

hdxlab = cell(1, lenhdx);
for i = 1:lenhdx
    hdxlab{i} = num2str(hdx(i));
end
p = [2 3 5 6];
for l = 1:4
    switch l
        case 1
            idx = stmsign_p2;
            tlab = 'pref stm';
        case 2
            idx = stmsign_n2;
            tlab = 'null stm';
        case 3
            idx = ch==1;
            tlab = 'ch 1';
        case 4
            idx = ch==0;
            tlab = 'ch 0';
    end
    c = zeros(nbin, lenhdx);
    for n = 1:nbin
        for m = 1:lenhdx
            c(n,m) = sum(stmmat(idx,n)==hdx(m));
        end
    end
    subplot(2,3,p(l))
    imagesc(1:lenhdx, 1:nbin, c)
    colorbar
    xlabel('hdx')
    ylabel('pulse')
    title(tlab)
    set(gca,'XTick',1:lenhdx, 'XTickLabel',hdxlab)
    set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
end

%%
% debug psth
figure;
for l = 1:2
    switch l
        case 1
            idx1 = stmsign_p2;
            idx2 = stmsign_n2;
            tlab1 = 'pref stm';
            tlab2 = 'null stm';
        case 2
            idx1 = ch==1;
            idx2 = ch==0;
            tlab1 = 'pref ch';
            tlab2 = 'null ch';
    end
    subplot(1,2,l)
    plot(time, mean(fr1(idx1,:)), '-b', 'linewidth',1)
    hold on;
    plot(time, mean(fr1(idx2,:)), '-r', 'linewidth',1)
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
    ylabel('psth')
    legend(tlab1, tlab2)
    legend('boxoff')
    set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
end

%%
% debug PK
figure;

subplot(1,2,1)
imagesc(1:nbin, hdx, tkernel/max(tkernel(:)))
colormap(copper)
c = caxis;
xlabel('time bin')
ylabel('hdx')
title(c)

subplot(1,2,2)
nom = mean(amp);
if resample_flag==1
    fill_between(1:nbin, (amp-err)/nom, (amp+err)/nom, [0 0 0])
    hold on;
end
plot(1:nbin, amp/nom, '-', 'color', [0 0 0], 'linewidth', 2)
xlim([0.5 nbin + 0.5])
ylabel({'kernel amplitude', ['(' pka_method ')']})
set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')

% yellow and green
y = [0.9576    0.7285    0.2285];
g = [0.1059    0.4706    0.2157];
figure;
%     subplot(1,4,1)
%     imagesc(1:nbin, hdx, tkernel)
%     c0 = caxis;
%     xlabel('time bin')
%     ylabel('hdx')
%     title('overall')

subplot(1,3,1)
vmax = max([max(tkernel_h(:)), max(tkernel_l(:))]);
imagesc(1:nbin, hdx, tkernel_h/vmax)
colormap(copper)
c1 = caxis;
xlabel('time bin')
ylabel('hdx')
title('high conf.')

subplot(1,3,2)
imagesc(1:nbin, hdx, tkernel_l/vmax)
colormap(copper)
c2 = caxis;
xlabel('time bin')
ylabel('hdx')
title('low conf.')

c = [min([c1 c2]) max([c1 c2])];
subplot(1,3,1)
caxis(c)
subplot(1,3,2)
caxis(c)

subplot(1,3,3)
%     nom = mean([amph, ampl]);
nom = max(amp);
if resample_flag==1
%         fill_between(1:nbin, (amp-err)/nom, (amp+err)/nom, [1 0 0])
%         hold on;
    fill_between(1:nbin, (amph-errh)/nom, (amph+errh)/nom, y)
    hold on;
    fill_between(1:nbin, (ampl-errl)/nom, (ampl+errl)/nom, g)
    hold on;
end
%     plot(1:nbin, amp/nom, '-r', 'linewidth', 2)
%     hold on;
plot(1:nbin, amph/nom, '-', 'color', y, 'linewidth', 2)
hold on;
plot(1:nbin, ampl/nom, '-', 'color', g, 'linewidth', 2)
xlim([0.5 nbin + 0.5])
ylabel({'kernel amplitude', ['(' pka_method ')']})
title(c)
set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')

function fill_between(x,y_bottom, y_top, maincolor,transparency,varargin)
if nargin < 3
        error('x, y_bottom, y_top are required as input arguments')
elseif nargin==3
        maincolor = [0 0 0];
        transparency = [];
elseif nargin==4
        transparency = [];
end

edgecolor = maincolor + (1 - maincolor)*0.55;

h = fill([x fliplr(x)],[y_bottom fliplr(y_top)],edgecolor);
set(h,'EdgeColor','none')
if ~isnan(transparency)
        alpha(transparency)
end