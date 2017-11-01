function Yates_simple_example
%% "stimulus-to-SensoryNeuron model' as presented in
%% Yates et al., 2017
%
% written by Katsuhisa (10.2017)
% +++++++++++++++++++++++++++++++++++++++++++++++

% input arguments
len_tr = 5000;
tmax = 400;
tau = 40;
kernelgain_s = 0.03;
kernelgain_c = 0.05;

hdx = 0.3*[-1 -0.75 -0.5 -0.25 0 0.25 0.5 0.75 1];
% hdx = 0.3*[-1 -0.5 -0.25 -0.125 0 0.125 0.25 0.5 1];
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
hn_offset = [ones(1,round(t(end)/3))*0.025, 0.025:-0.025/(round(t(end)*2/3)):0]; % manually tweaked to approximate that in Yates
kernel1 = kernel1 - 1.5*hn_offset;
kernel1(1:4) = 0;
kernel2 = kernel1;
kernel1 = kernelgain_c*kernel1/max(kernel1) ;
kernel2 = kernelgain_s*kernel2/max(kernel2);

time = [1:len_frame+2*offset] - offset;    

%%
% generate dynamic stimulus sequence with the 0% signal
frameperbin = len_frame/nbin;
[stm, stmmat] = Yates_stm(hdx, len_tr, nbin, frameperbin, 1220);

% overall stimulus sign
sumstm = sum(stm,2);
stmsign_p2 = sign(sumstm) > 0 & sumstm > median(sumstm(sign(sumstm) > 0));
stmsign_n2 = sign(sumstm) < 0 & sumstm < median(sumstm(sign(sumstm) < 0));

% include offset
stm = [zeros(len_tr, offset), stm, zeros(len_tr, offset)];
disp('stimulus generated')

%%
% neural activity
lenv = size(stm,2);
c = conv(kernel1, co);
fr1 = zeros(len_tr, lenv);
fr2 = zeros(len_tr, lenv);
% res1 = zeros(len_tr, lenv);
% res2 = zeros(len_tr, lenv);
for i = 1:len_tr
    s1 = conv(kernel2, stm(i,:));
    s2 = conv(-kernel2, stm(i,:));
    fr1(i,:) = exp(s1(1:lenv)+c(1:lenv));
    fr2(i,:) = exp(s2(1:lenv)+c(1:lenv));
%     for f = 1:lenv
%         res1(i,f) = poissrnd(fr1(i,f), 1, 1);
%         res2(i,f) = poissrnd(fr2(i,f), 1, 1);
%     end
end
disp('spikes generated')

%%
% evidence
dv1 = cumsum(fr1(:, offset+1:offset+len_frame),2);
dv2 = cumsum(fr2(:, offset+1:offset+len_frame),2);
ev = dv1(:,end) - dv2(:, end);
disp('decision variables computed')

%% 
% choice
ch = sign(ev);
ch(ch==0) = datasample([-1 1], sum(ch==0), 'Replace', true);
ch(ch==-1) = 0;
disp(['ch 1: ' num2str(sum(ch==1)) ', ch2: ' num2str(sum(ch==0))])

%%
% confidence
conf = abs(ev);
med = median(conf);

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
    [tkernel(:,a)] = getKernel(hdx, stm(:, begin:begin+frameperbin-1), ch);
    [tkernel_h(:,a)] = getKernel(hdx, stm(conf > med, begin:begin+frameperbin-1), ch(conf > med));
    [tkernel_l(:,a)] = getKernel(hdx, stm(conf < med, begin:begin+frameperbin-1), ch(conf < med));
    
    begin = begin + frameperbin;
end
for a = 1:nbin
    amp(a) = tkernel(:,a)'*mean(tkernel,2);
    amph(a) = tkernel_h(:,a)'*mean(tkernel_h,2);
    ampl(a) = tkernel_l(:,a)'* mean(tkernel_l,2);
end
disp('PK computed')

%% 
close all;

% debug stimulus statistics
figure;
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
% yellow and green
y = [0.9576    0.7285    0.2285];
g = [0.1059    0.4706    0.2157];
figure;
subplot(1,4,1)
imagesc(1:nbin, hdx, tkernel)
c0 = caxis;
xlabel('time bin')
ylabel('hdx')
title('overall')

subplot(1,4,2)
imagesc(1:nbin, hdx, tkernel_h)
c1 = caxis;
xlabel('time bin')
ylabel('hdx')
title('high conf.')

subplot(1,4,3)
imagesc(1:nbin, hdx, tkernel_l)
c2 = caxis;
xlabel('time bin')
ylabel('hdx')
title('low conf.')

c = [min([c0 c1 c2]) max([c0 c1 c2])];
subplot(1,4,1)
caxis(c)
subplot(1,4,2)
caxis(c)
subplot(1,4,3)
caxis(c)

subplot(1,4,4)
nom = mean([amph, ampl]);
plot(1:nbin, amp/nom, '-r', 'linewidth', 2)
hold on;
plot(1:nbin, amph/nom, '-', 'color', y, 'linewidth', 2)
hold on;
plot(1:nbin, ampl/nom, '-', 'color', g, 'linewidth', 2)
xlim([0.5 nbin + 0.5])
ylabel('kernel amplitude')
set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')

%% 
% subfunction
function [stm, stmmat] = Yates_stm(hdx, len_tr, nbin, frameperbin, seed)
rng(seed)
stmidx = randi(length(hdx),len_tr,nbin);
stmmat = hdx(stmidx);
stm = reshape(repmat(stmmat,frameperbin,1),len_tr,nbin*frameperbin);

function [pk0] = getKernel(disval, hdxmat, ch)
len_d = length(disval);
svmat = zeros(size(hdxmat,1), len_d);
for r = 1:size(hdxmat,1)
    for d = 1:len_d
        svmat(r,d) = sum(hdxmat(r,:)==disval(d));
    end
end

pk0 = mean(svmat(ch==0,:),1) - mean(svmat(ch==1,:),1);