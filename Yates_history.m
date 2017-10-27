function [para] = Yates_history(varargin)

close all;

% input arguments
nneuron = 10;
len_tr = 500;
tmax = 400;
tau = 40;
kernelgain_s = 0.04;
kernelgain_c = 0.04;
plot_flag = 1;
j = 1;              
while  j <= length(varargin)
    switch varargin{j}
        case 'ntr'
            len_tr = varargin{j+1};
            j = j + 2;
        case 'tmax' 
            tmax = varargin{j+1};
            j = j + 2;
        case 'tau'
            tau = varargin{j+1};
            j = j + 2;
        case 'kernelgain_s'
            kernelgain_s = varargin{j+1};            
             j = j + 2;
        case 'kernelgain_c'
            kernelgain_c = varargin{j+1};            
             j = j + 2;
        case 'nofig'
            plot_flag = 0;
            j = j + 1;
    end
end

% yellow and green
y = [0.9576    0.7285    0.2285];
g = [0.1059    0.4706    0.2157];

hdx = 0.3*[-1 -0.75 -0.5 -0.25 0 0.25 0.5 0.75 1];
% hdx = 0.2*[-1 -0.5 -0.25 -0.125 0 0.125 0.25 0.5 1];
len_frame = 1050;
lenhdx = length(hdx);
nbin = 7;

% contrast
offset = 100;
co = 1*[zeros(1,offset),ones(1,len_frame),zeros(1,offset)];

% alpha function as a stimulus kernel
t = 0:tmax;
kernel1 = exp(-t/tau).*(1 - exp(-t/tau));
hn_offset = [ones(1,round(t(end)/3))*0.025, 0.025:-0.025/(round(t(end)*2/3)):0]; % manually tweaked to approximate that in Yates
kernel1 = kernel1 - 1.5*hn_offset;
kernel1(1:4) = 0;
kernel2 = kernel1;
kernel1 = kernelgain_c*kernel1/max(kernel1) ;
kernel2 = kernelgain_s*kernel2/max(kernel2);
    
para.kernel_stm = kernel2;
para.kernel_co = kernel1;

% kernel for the history term
h = 0:200;
% kernel3 = log(1+h);
% kernel3 = normalize(kernel3, -0.01, 0);
tau3 = 20;
kernel3 = exp(-h/tau3).*(1 - exp(-h/tau3));
hn_offset = [ones(1,round(h(end)/3))*0.025, 0.025:-0.025/(round(h(end)*2/3)):0]; % manually tweaked to approximate that in Yates
kernel3 = kernel3 - 3*hn_offset;
kernel3(1:2) = 0;
kernel3 = normalize(-kernel3,0, 1);

% figure;
% plot(kernel3)

time = [1:len_frame+2*offset] - offset;
    
if plot_flag == 1
    figure(1);
    subplot(2,4,1)
    plot(t, kernel2, '-b', 'linewidth',2)
    hold on;
    plot(t, -kernel2, '-r', 'linewidth',2)
    xlim([t(1) t(end)])
    title('stimulus kernel')
    set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
    
    subplot(2,4,2)
    plot(t, kernel1, '-k', 'linewidth',2)
    xlim([t(1) t(end)])
    title('contrast kernel')
    set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
    
    subplot(2,4,3)
    plot(kernel3, '-k', 'linewidth',2)
%     xlim([t(1) t(end)])
    title('history kernel')
    set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
end

% generate dynamic stimulus sequence with the 0% signal
frameperbin = len_frame/nbin;
stm = nan(len_tr, len_frame);
for i = 1:len_tr
    stm_temp = datasample(hdx, nbin, 'Replace', true);
%     while abs(sum(stm_temp)) >= 0.1
%         stm_temp = datasample(hdx, nbin, 'Replace', true);
%     end
%  
    begin = 1;
    for n = 1:nbin
        stm(i, begin:begin+frameperbin-1) = stm_temp(n)*ones(1, frameperbin);
        begin = begin + frameperbin;
    end
end

% overall stimulus sign
sumstm = sum(stm,2);
stmsign = sign(sumstm);
med_p = median(sumstm(stmsign > 0));
med_n = median(sumstm(stmsign < 0));
% stmsign_p1 = stmsign > 0 & sumstm < med_p;
stmsign_p2 = stmsign > 0 & sumstm > med_p;
% stmsign_n1 = stmsign < 0 & sumstm > med_n;
stmsign_n2 = stmsign < 0 & sumstm < med_n;


% include offset
stm = [zeros(len_tr, offset), stm, zeros(len_tr, offset)];

% if plot_flag==1
%     subplot(2,4,2)
%     plot(time, mean(stm(:, 1:length(time)),1), '-k')
%     ylabel('mean stimulus')
% end

% sensory neurons' response
for n = 1:nneuron

    % convolution with stimulus kernels
    lenv = size(stm,2);
    para.neuron(n).spk1 = nan(len_tr, lenv);
    para.neuron(n).spk2 = para.neuron(n).spk1;
    c = conv(kernel1, co);
    s = conv(kernel2, stm(1,:));
    para.neuron(n).spk1(1,:) = random('Poisson', exp(s(1:lenv) + c(1:lenv)), 1, lenv);
    s = conv(-kernel2, stm(1,:));
    para.neuron(n).spk2(1,:) = random('Poisson', exp(s(1:lenv) + c(1:lenv)), 1, lenv);

    % for debugging
    if n == 1
        figure(2);
        subplot(2,5,1)
        plot(kernel1);
        title('contrast')

        subplot(2,5,2)
        plot(kernel2,'-b')
        hold on;
        plot(-kernel2, '-r')
        title('stimulus')

        subplot(2,5,3)
        plot(kernel3)
        title('history')

        subplot(2,5,6)
        c = conv(kernel1, co);
        plot(c(1:lenv));
        hold on;
        plot(co)

        subplot(2,5,7)
        s = conv(kernel2, stm(2,:));
        plot(s(1:lenv))
        hold on;
        plot(stm(2,:))

        subplot(2,5,8)
        h = conv(kernel3, para.neuron(n).spk1(1,:));
        plot(h(1:lenv))
        hold on;
        plot(para.neuron(n).spk1(1,:))

        subplot(2,5,9)
        plot(c(1:lenv) + s(1:lenv) + h(1:lenv))
        title('sum of all')

        subplot(2,5,10)
        plot(exp(c(1:lenv) + s(1:lenv) + h(1:lenv)))
        title('nonlinearity')
    end

    for i = 2:len_tr
        s = conv(kernel2, stm(i,:));
        h = conv(kernel3, para.neuron(n).spk1(i-1,:));
        para.neuron(n).spk1(i,:) = random('Poisson', exp(s(1:lenv) + c(1:lenv) + h(1:lenv))/10, 1, lenv);
        s = conv(-kernel2, stm(i,:));
        h = conv(kernel3, para.neuron(n).spk2(i-1,:));
        para.neuron(n).spk2(i,:) = random('Poisson', exp(s(1:lenv) + c(1:lenv) + h(1:lenv))/10, 1, lenv);
    end
end

% % noise
% cvstm1(:, offset+1:offset+len_frame) = ...
%     cvstm1(:, offset+1:offset+len_frame) + normrnd(0, 5, [len_tr, len_frame]);
% cvstm2(:, offset+1:offset+len_frame) = ...
%     cvstm2(:, offset+1:offset+len_frame) + normrnd(0, 5, [len_tr, len_frame]);

% evidence
dv1 = zeros(len_tr, len_frame);
dv2 = zeros(len_tr, len_frame);
for n = 1:nneuron
    dv1 = dv1 + cumsum(para.neuron(n).spk1(:, offset+1:offset+len_frame),2);
    dv2 = dv2 + cumsum(para.neuron(n).spk2(:, offset+1:offset+len_frame),2);
    
    % mean firing rate
    para.neuron(n).mfr1 = mean(para.neuron(n).spk1(:));
    para.neuron(n).mfr2 = mean(para.neuron(n).spk2(:));
end

ev = dv1(:,end) - dv2(:, end);

% choice
ch = sign(ev);
ch(ch==0) = datasample([-1 1], sum(ch==0), 'Replace', true);
ch(ch==-1) = 0;

disp(['far choice: ' num2str(sum(ch==1)) ...
    ', near choice: ' num2str(sum(ch==0))])

% % normalization
% cvstm1 = normalize(cvstm1, 0, 1);
% cvstm2 = normalize(cvstm2, 0, 1);

% debug psth
k = datasample(1:nneuron, 1);
begin = offset + 1;
psth_stm = nan(1, nbin);
psth_ch = nan(1, nbin);
for i = 1:nbin
    psth_stm(i) = sum(mean(para.neuron(k).spk1(stmsign_p2, begin:begin+frameperbin-1),1))...
    	- sum(mean(para.neuron(k).spk1(stmsign_n2, begin:begin+frameperbin-1),1));
    psth_ch(i) = sum(mean(para.neuron(k).spk1(ch==1, begin:begin+frameperbin-1),1))...
    	- sum(mean(para.neuron(k).spk1(ch==0, begin:begin+frameperbin-1),1));
    begin = begin + frameperbin;
end

para.psth.stm_diff = psth_stm;
para.psth.ch_diff = psth_ch;
para.psth.stm_raw0 = mean(para.neuron(k).spk1(stmsign_n2, 1:length(time)),1);
para.psth.stm_raw1 = mean(para.neuron(k).spk1(stmsign_p2, 1:length(time)),1);
para.psth.ch_raw0 = mean(para.neuron(k).spk1(ch==1, 1:length(time)),1);
para.psth.ch_raw1 = mean(para.neuron(k).spk1(ch==0, 1:length(time)),1);


if plot_flag==1    
        
    figure(1);
    subplot(2,4,[5 6])
    % raster
%     rng(19891220)
    tr = datasample(1:len_tr,1);
    ras1 = zeros(nneuron, len_frame+2*offset);
    ras2 = zeros(nneuron, len_frame+2*offset);
    spkch1 = zeros(1, len_frame + 2*offset);
    spkch2 = zeros(1, len_frame + 2*offset);
    for n = 1:nneuron
            ras1(n,:) = 1*(para.neuron(n).spk1(tr,:) > 0);
            ras2(n,:) = 1*(para.neuron(n).spk2(tr,:) > 0);
            spkch1 = spkch1 + mean(para.neuron(n).spk1(ch==1, 1:(len_frame+2*offset)), 1);
            spkch2 = spkch2 + mean(para.neuron(n).spk1(ch==0, 1:(len_frame+2*offset)), 1);
    end
    imagesc(time, 1:2*nneuron, [ras2; ras1])
    colormap(gca, flipud(bone))
    hold on;
    plot(time, nneuron*ones(1, length(time)), '-m')
    hold on;
    plot(time, nneuron+ 3*nneuron*stm(tr,:), '-g', 'linewidth',2)
    hold on;
    yy = get(gca, 'YLim');
    begin = 1;
    for i = 1:nbin+1
        hold on;
        plot(begin*[1 1],yy, '-m')
        begin = begin + frameperbin;
    end
    xlim([-offset len_frame+offset])
    ylim([0 2*nneuron+1])    
    ylabel('neurons')
    set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')

    % psth
    subplot(2,4,6)
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
    ylabel('psth (by choice)')
    set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
    
    subplot(2,4,7)
    plot(time, spkch1/nneuron, '-b', 'linewidth',1)
    hold on;
    plot(time, spkch2/nneuron, '-r', 'linewidth',1)
    hold on;
%     plot(time, mean(cvstm2(ch==0, 1:length(time)),1), '--b', 'linewidth',2)
%     hold on;
%     plot(time, mean(cvstm2(ch==1, 1:length(time)),1), '--r', 'linewidth',2)
    yy = get(gca, 'YLim');
    begin = 1;
    for i = 1:nbin+1
        hold on;
        plot(begin*[1 1],yy, ':k')
        begin = begin + frameperbin;
    end
    xlim([-offset len_frame+offset])
    ylim(yy)    
    ylabel('psth (by choice)')
    set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
end

% confidence
conf = abs(ev);
% conf = conf + normrnd(median(conf), 1*median(conf), size(conf));
med = median(conf);

% debug confidence
figure(4);
histogram(conf)
hold on;
yy = get(gca, 'YLim');
plot(med*[1 1], yy, '-r')
xlabel('confidence')

disp(['median confidence: ' num2str(med)])

% if plot_flag==1
%     subplot(2,4,5)
%     plot(1:size(dv1,2), mean(abs(dv1(ch==1,:)), 1), '-r', 'linewidth', 2)
%     hold on;
%     plot(1:size(dv1,2), mean(abs(dv1(conf > med & ch==1, :)), 1), '-', 'color', y, 'linewidth', 2)
%     hold on;
%     plot(1:size(dv1,2), mean(abs(dv1(conf < med & ch==1, :)), 1), '-', 'color', g, 'linewidth', 2)
%     hold on;
%     plot(1:size(dv2,2), mean(abs(dv2(ch==0,:)), 1), '--r', 'linewidth', 2)
%     hold on;
%     plot(1:size(dv2,2), mean(abs(dv2(conf > med & ch==0, :)), 1), '--', 'color', y, 'linewidth', 2)
%     hold on;
%     plot(1:size(dv2,2), mean(abs(dv2(conf < med & ch==0, :)), 1), '--', 'color', g, 'linewidth', 2)
%     ylabel('mean |decision variable|')
%     xlim([0 size(dv1, 2)+1])
%     xlabel('time during stimulus presentation')
% end

% % time-averaged kernel
% pk0 = getKernel(stm(:, offset+1:offset+len_frame), ch);
% pkh = getKernel(stm(conf > med, offset+1:offset+len_frame), ch(conf > med));
% pkl = getKernel(stm(conf < med, offset+1:offset+len_frame), ch(conf < med));

% time-resolved kernel
tkernel = nan(lenhdx, nbin);
tkernel_h = nan(lenhdx, nbin);
tkernel_l = nan(lenhdx, nbin);

% kernel amplitude
amp = nan(1, nbin);
amph = nan(1, nbin);
ampl = nan(1, nbin);

begin = offset+1;
for a = 1:nbin
%     disp(['bin ' num2str(a) ': ' num2str(begin - offset) ' to ' num2str(begin+frameperbin-1-offset)])
    [tkernel(:,a)] = getKernel(stm(:, begin:begin+frameperbin-1), ch);
    [tkernel_h(:,a)] = getKernel(stm(conf > med, begin:begin+frameperbin-1), ch(conf > med));
    [tkernel_l(:,a)] = getKernel(stm(conf < med, begin:begin+frameperbin-1), ch(conf < med));
    
    begin = begin + frameperbin;
end

for a = 1:nbin
    amp(a) = tkernel(:,a)'*mean(tkernel,2);
    amph(a) = tkernel_h(:,a)'*mean(tkernel_h,2);
    ampl(a) = tkernel_l(:,a)'* mean(tkernel_l,2);
end

para.amp_h = amph;
para.amp_l = ampl;
para.amp_diff = amph - ampl;

% kernel amplitude using logistic regression
hdxmat = zeros(len_tr, nbin);
begin = offset + 1;
for n = 1:nbin
    hdxmat(:, n) = median(stm(:, begin:begin+frameperbin-1), 2);
    begin = begin + frameperbin;
end
w = logregPK(hdxmat, ch);
wh = logregPK(hdxmat(conf > med, :), ch(conf > med));
wl = logregPK(hdxmat(conf < med, :), ch(conf < med));
wn = mean([wh, wl]);

if plot_flag==1
    figure(1);
    % visualize kernels
%     subplot(2,4,4)
    subplot(2,4,4)
    imagesc(1:nbin, hdx, tkernel)
    colormap(gca, 'parula')
    xlabel('time bin')
    ylabel('hdx')
    set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
    
    % visualize kernel amplitude
%     subplot(2,4,8)
    subplot(2,4,8)
    nom = mean([amph, ampl]);
    plot(1:nbin, amp/nom, '-r', 'linewidth', 2)
    hold on;
    plot(1:nbin, amph/nom, '-', 'color', y, 'linewidth', 2)
    hold on;
    plot(1:nbin, ampl/nom, '-', 'color', g, 'linewidth', 2)
    hold on;
    plot(1:nbin, w/wn, '--r', 'linewidth', 2)
    hold on;
    plot(1:nbin, wh/wn, '--', 'color', y, 'linewidth', 2)
    hold on;
    plot(1:nbin, wl/wn, '--', 'color', g, 'linewidth', 2)
    xlim([0.5 nbin + 0.5])
    ylabel('kernel amplitude')
    set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')

    % cosmetics
%     subplot(2,4,1)
%     ylabel('hdx')
%     xlabel('time (ms)')
%     xlim([-offset offset+len_frame])
%     set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
% 
%     subplot(2,4,2)
%     xlim([-offset offset+len_frame])
%     set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')

%     subplot(2,4,3)
%     subplot(2,3,[1 2])
%     set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
%     legend('location','eastoutside')

%     subplot(2,4,4)
%     set(gca, 'XTick', 1:nbin)
%     set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')

%     subplot(2,4,5)
%     % xlim([-offset offset+len_frame])
%     set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
% 
% %     subplot(2,4,6)
%     subplot(2,3,4)
%     xlim([-offset offset+len_frame])
%     set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
% 
% %     subplot(2,4,7)
%     subplot(2,3,5)
%     xlim([-offset offset+len_frame])
%     set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
% 
% %     subplot(2,4,8)
%     subplot(2,3,6)
%     set(gca, 'XTick', 1:nbin)
%     set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')

    % additional figures for paper
    figure(10);
    subplot(1,3,1)
    imagesc(1:nbin, hdx, tkernel_h)
    colormap(pink)
    c_h = caxis;
    xlabel('time bin')
    ylabel('hdx')
    set(gca, 'XTick', 1:nbin)
    set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')

    subplot(1,3,2)
    imagesc(1:nbin, hdx, tkernel_l)
    colormap(pink)
    c_l = caxis;
    set(gca, 'XTick', 1:nbin)
    set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')

    cfix = [min([c_h, c_l]), max([c_h, c_l])];
    subplot(1,3,1)
    caxis(cfix)
    subplot(1,3,2)
    caxis(cfix)

    subplot(1,3,3)
    nom = mean([amph, ampl]);
    plot(1:nbin, amph/nom, '-', 'color', y, 'linewidth', 2)
    hold on;
    plot(1:nbin, ampl/nom, '-', 'color', g, 'linewidth', 2)
    xlim([0.5 nbin + 0.5])
    ylabel('kernel amplitude')
    set(gca, 'XTick', 1:nbin)
    set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
end

function [pk0] = getKernel(hdxmat, ch)

% trial averaged stimulus distributions
disval = unique(hdxmat);
len_d = length(disval);
svmat = zeros(size(hdxmat,1), len_d);
for r = 1:size(hdxmat,1)
    for d = 1:len_d
        svmat(r,d) = sum(hdxmat(r,:)==disval(d));
    end
end

% compute PK for 0% stimulus
pk0 = mean(svmat(ch==0,:),1) - mean(svmat(ch==1,:),1);

function w = logregPK(hdxmat, ch)
w = size(hdxmat, 2);
for i = 1:size(hdxmat, 2)
    b = glmfit(hdxmat(:, i), ch, 'binomial', 'link', 'logit', 'constant', 'on');
    w(i) = b(2);
end

function [normalized_vector] = normalize(v, newmin, newmax)

a = (newmax - newmin)/(max(max(v)) - min(min(v)));
b = newmax - max(max(v))*a;

normalized_vector = a.*v + b;

