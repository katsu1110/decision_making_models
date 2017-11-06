function [para] = Yates_toy(varargin)
%% Yates toy model ("stimulus-to-MT model") without history term
% The similar one in the supplementary at Yates et al., 2017
%
% written by Katsuhisa (02.11.17)
% ++++++++++++++++++++++++++++++++++++++++++++++++++++

close all;

% input arguments
plot_flag = 1;
nneuron = 10;
len_tr = 100;
tmax = 150;
tau = 10;
kernelgain_s = 0.08;
kernelgain_c = 0.09;

hdx = 0.3*[-1 -0.75 -0.5 -0.25 0 0.25 0.5 0.75 1];
lenhdx = length(hdx);
len_frame = 1050;
nbin = 7;

% contrast
offset = 100;
co = 1*[zeros(1,offset),ones(1,len_frame),zeros(1,offset)];

%%
% alpha function as a stimulus kernel
t = 0:tmax;
kernel1 = exp(-t/tau).*(1 - exp(-t/tau));
kernel2 = kernel1;
kernel1 = kernelgain_c*kernel1/max(kernel1) ;
kernel2 = kernelgain_s*kernel2/max(kernel2);
    
para.kernel_stm = kernel2;
para.kernel_co = kernel1;

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
    
    figure(12345);
    subplot(2,3,1)
    plot(t, kernel2, '-b', 'linewidth',2)
    hold on;
    plot(t, -kernel2, '-r', 'linewidth',2)
    xlim([t(1) t(end)])
    title('stimulus kernel')
    set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
    
    subplot(2,3,2)
    plot(t, kernel1, '-k', 'linewidth',2)
    xlim([t(1) t(end)])
    title('contrast kernel')
    set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
    
end
%%
% generate dynamic stimulus sequence with the 0% signal
frameperbin = len_frame/nbin;
[stm] = Yates_stm(hdx, len_tr, nbin, frameperbin, 1220);

% overall stimulus sign
sumstm = sum(stm,2);
stmsign = sign(sumstm);
med_p = median(sumstm(stmsign > 0));
med_n = median(sumstm(stmsign < 0));
stmsign_p2 = stmsign > 0 & sumstm > med_p;
stmsign_n2 = stmsign < 0 & sumstm < med_n;

disp('stimulus generated')

% include offset
stm = [zeros(len_tr, offset), stm, zeros(len_tr, offset)];

%%
% sensory neurons' responses
lenv = len_frame + 2*offset;
time = [1:lenv] - offset;
c = conv(kernel1, co);    
for i = 1:len_tr       
    % convolution with stimulus
    s1 = conv(kernel2, stm(i,:));
    para.tr(i).spk1 = arrayfun(@(x) poissrnd(x), repmat(exp(s1(1:lenv) + c(1:lenv)),nneuron,1));
    s2 = conv(-kernel2, stm(i,:));
    para.tr(i).spk2 = arrayfun(@(x) poissrnd(x), repmat(exp(s2(1:lenv) + c(1:lenv)),nneuron,1));
    
    % debug psth
    if i==1
        figure(12345); 
        subplot(2,3,[4 5 6])
        imagesc(time, 1:2*nneuron, [para.tr(i).spk1; para.tr(i).spk2])
        colormap(gca, flipud(bone))
        hold on;
        plot(time, (nneuron+0.5)*ones(1, length(time)), '-m')
        hold on;
        plot(time, (nneuron+0.5)+ 3*nneuron*stm(i,:), '-g', 'linewidth',2)
        hold on;
        xlim([-offset len_frame+offset])
        ylim([0 2*nneuron+1])   
        xlabel('time (ms)')
        ylabel('neurons')
        set(gca, 'TickDir', 'out')
    end
end
disp('spikes generated')

% reshape struct
for n = 1:nneuron
    for i = 1:len_tr
        para.neuron(n).spk1(i,:) = para.tr(i).spk1(n,:);
        para.neuron(n).spk2(i,:) = para.tr(i).spk2(n,:);
    end
end

%%
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

%%
% choice
ch = sign(ev);
ch(ch==0) = datasample([-1 1], sum(ch==0), 'Replace', true);
ch(ch==-1) = 0;

disp(['far choice: ' num2str(sum(ch==1)) ...
    ', near choice: ' num2str(sum(ch==0))])

%%
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
para.psth.stm_null = mean(para.neuron(k).spk1(stmsign_n2, 1:length(time)),1);
para.psth.stm_pref = mean(para.neuron(k).spk1(stmsign_p2, 1:length(time)),1);
para.psth.ch_pref = mean(para.neuron(k).spk1(ch==1, 1:length(time)),1);
para.psth.ch_null = mean(para.neuron(k).spk1(ch==0, 1:length(time)),1);

%%
% visualize neural activity
if plot_flag==1    
        
    figure(1);
    subplot(2,4,5)
    % raster
%     rng(19891220)
    tr = datasample(1:len_tr,1);
    ras1 = zeros(nneuron, len_frame+2*offset);
    ras2 = zeros(nneuron, len_frame+2*offset);
    spkch1 = zeros(1, len_frame + 2*offset);
    spkch2 = zeros(1, len_frame + 2*offset);
    spkstm1 = zeros(1, len_frame + 2*offset);
    spkstm2 = zeros(1, len_frame + 2*offset);
    for n = 1:nneuron
            ras1(n,:) = 1*(para.neuron(n).spk1(tr,:) > 0);
            ras2(n,:) = 1*(para.neuron(n).spk2(tr,:) > 0);
            spkch1 = spkch1 + mean(para.neuron(n).spk1(ch==1, 1:(len_frame+2*offset)), 1);
            spkch2 = spkch2 + mean(para.neuron(n).spk1(ch==0, 1:(len_frame+2*offset)), 1);
            spkstm1 = spkstm1 + mean(para.neuron(n).spk1(stmsign_p2, 1:(len_frame+2*offset)), 1);
            spkstm2 = spkstm2 + mean(para.neuron(n).spk1(stmsign_n2, 1:(len_frame+2*offset)), 1);
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
    ylabel('psth (by stm)')
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
    
    % PTA
    figure(234);
    col = copper(nbin);
    PTA = cell(1, nbin);
    for n = 1:nneuron
        [PTA1] = PTA_easy(para.neuron(n).spk1, stm, max(hdx),nbin, offset, frameperbin);
        [PTA2] = PTA_easy(para.neuron(n).spk2, stm, min(hdx),nbin, offset, frameperbin);
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
    xlabel('time')
    ylabel('PTA')
    set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
    
end

%%
% confidence
conf = abs(ev);
% conf = conf + normrnd(median(conf), 1*median(conf), size(conf));
med = median(conf);

%%
% debug confidence
figure(4);
histogram(conf)
hold on;
yy = get(gca, 'YLim');
plot(med*[1 1], yy, '-r')
xlabel('confidence')

disp(['median confidence: ' num2str(med)])

%%
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

%%
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

%%
if plot_flag==1
    % yellow and green
    y = [0.9576    0.7285    0.2285];
    g = [0.1059    0.4706    0.2157];

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

%% 
% subfunctions
function [stm, stmmat] = Yates_stm(hdx, len_tr, nbin, frameperbin, seed)
rng(seed)
stmidx = randi(length(hdx),len_tr,nbin);
stmmat = hdx(stmidx);
stm = reshape(repmat(stmmat,frameperbin,1),len_tr,nbin*frameperbin);


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

function [PTA] = PTA_easy(spk, stm, pref, nbin, offset, frameperbin)
PTA = cell(1, nbin);
begin = 1 + offset;
for n = 1:nbin
    PTA{n} = mean(spk(stm(:,begin)==pref,1+offset:nbin*frameperbin+offset),1)...
        - mean(spk(:,1+offset:nbin*frameperbin+offset),1);
    begin = begin + frameperbin;
end


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