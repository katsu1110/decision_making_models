function [para] = Yates_simple_example(varargin)
%% "stimulus-to-SensoryNeuron model' as presented in
%% Yates et al., 2017
%
% written by Katsuhisa (10.2017)
% +++++++++++++++++++++++++++++++++++++++++++++++

% input arguments
len_tr = 50000;
tmax = 300; %was 400
tau = 50; %was 30
kernelgain_s = 0.05; %was 0.05
kernelgain_c = 0.025; %was 0.05
offset_gain = 1.5; %was 1
stm_gain = 0.3;
logreg_flag = 0;
plot_flag = 1;
resample_flag = 0;
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
        case 'offset_gain'
            offset_gain = varargin{j+1};            
             j = j + 2;
        case 'stm_gain'
            stm_gain = varargin{j+1};               
             j = j + 2;
        case 'logreg'
            logreg_flag = 1;
            j = j + 1;
        case 'resample'
            resample_flag = 1;
            j = j + 1;
        case 'nofig'
            plot_flag = 0;
            j = j + 1;
    end
end
hdx = stm_gain*[-1 -0.75 -0.5 -0.25 0 0.25 0.5 0.75 1];
% hdx = 0.3*[-1 -0.5 -0.25 -0.125 0 0.125 0.25 0.5 1];
nbin = 4;
if nbin==4
    len_frame = 600;
elseif nbin==7
    len_frame = 1050;
end
lenhdx = length(hdx);
hdxlab = cell(1, lenhdx);
for i = 1:lenhdx
    hdxlab{i} = num2str(hdx(i));
end

% contrast
offset = 100;
co = 1*[zeros(1,offset),ones(1,len_frame),zeros(1,offset)];

%%
% alpha function as a stimulus kernel
t = 0:tmax;
kernel1 = exp(-t/tau).*(1 - exp(-t/tau));
hn_offset_s = [ones(1,round(t(end)/3))*0.025, 0.025:-0.025/(round(t(end)*2/3)):0]; % manually tweaked to approximate that in Yates
hn_offset_c = [ones(1,length(t))*0.025]; % manually tweaked to approximate that in Yates for contrast kernel
kernel2 = kernel1 - offset_gain*hn_offset_s; % initially offset_gain was fixed at 1.5
kernel2(1:4) = 0;
kernel1 = kernel1-offset_gain*hn_offset_c;
kernel1(1:4) = 0;
kernel1 = kernelgain_c*kernel1/max(kernel1) ; % contrast kernel
kernel2 = kernelgain_s*kernel2/max(kernel2);  % stimulus kernel
    
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
% disp('spikes generated')

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

% % add noise on confidence judgement
% conf = conf + normrnd(median(conf), 0.2*median(conf), size(conf));
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
stmbin = nan(size(stm, 1), nbin);
begin = offset+1;
for a = 1:nbin
    stmbin(:, a) = mean(stm(:, begin:begin+frameperbin-1), 2);
    [tkernel(:,a)] = getKernel(hdx, stm(:, begin:begin+frameperbin-1), ch);
    [tkernel_h(:,a)] = getKernel(hdx, stm(conf > med, begin:begin+frameperbin-1), ch(conf > med));
    [tkernel_l(:,a)] = getKernel(hdx, stm(conf < med, begin:begin+frameperbin-1), ch(conf < med));
    
    begin = begin + frameperbin;
end
tapk = mean(tkernel,2)';
if logreg_flag==0
    % image classification
    pka_method = 'image classification';
    for a = 1:nbin
        amp(a) = tapk*tkernel(:,a);
        amph(a) = tapk*tkernel_h(:,a);
        ampl(a) = tapk*tkernel_l(:,a);
    end
else
    % logistic regression
    pka_method = 'logistic regression';
    amp = glmfit(stmbin, ch, ...
        'binomial', 'link', 'logit', 'constant', 'on');
    amph = glmfit(stmbin(conf > med, :), ch(conf > med), ...
        'binomial', 'link', 'logit', 'constant', 'on');
    ampl = glmfit(stmbin(conf < med, :), ch(conf < med), ...
        'binomial', 'link', 'logit', 'constant', 'on');
    amp = amp(2:end); amph = amph(2:end); ampl = ampl(2:end);
end
if resample_flag==1
    repeat = 1000;
    [err, errl, errh] = resamplePK(hdx, stm, ch, offset, nbin, ...
        frameperbin, conf, med, repeat, logreg_flag);
end
disp('PK computed')


%% 
if plot_flag==1
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
end

%% 
% store variables
para.kernel_stm = kernel2;
para.kernel_co = kernel1;
para.psth.stm_null = mean(fr1(stmsign_n2, 1:length(time)),1);
para.psth.stm_pref = mean(fr1(stmsign_p2, 1:length(time)),1);
para.psth.stm_null_adaptation = compute_adaptation(para.psth.stm_null, offset, nbin, frameperbin);
para.psth.stm_pref_adaptation = compute_adaptation(para.psth.stm_pref, offset, nbin, frameperbin);
para.psth.ch_pref = mean(fr1(ch==1, 1:length(time)),1);
para.psth.ch_null = mean(fr1(ch==0, 1:length(time)),1);
[cp1, cp2] = compute_cp(fr1(:, offset+1:offset+len_frame), ch, nbin, frameperbin);
para.psth.cp = [cp1; cp2];
para.pka_method = pka_method;
para.amp = amp;
para.amp_h = amph;
para.amp_l = ampl;

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

function [err, err0, err1] = resamplePK(disval, hdxmat, ch, offset, nbin, frameperbin, conf, med, repeat, logreg_flag)
conf0 = find(conf < med);
conf1 = find(conf > med);
ampr = zeros(repeat, nbin);
ampr0 = zeros(repeat, nbin);
ampr1 = zeros(repeat, nbin);
hdxmatbin = nan(size(hdxmat, 1), nbin);
for r = 1:repeat
    rtr0 = datasample(conf0, length(conf0));
    rtr1 = datasample(conf1, length(conf1));
    rtr = [rtr0; rtr1];
    tkernel = nan(length(disval), nbin);
    tkernel0 = nan(length(disval), nbin);
    tkernel1 = nan(length(disval), nbin);
    begin = offset+1;   
    for a = 1:nbin
        if r==1
            hdxmatbin(:, n) = mean(hdxmat(:, begin:begin+frameperbin-1), 2);
        end
        [tkernel(:,a)] = getKernel(disval, hdxmat(rtr, begin:begin+frameperbin-1), ch(rtr));   
        [tkernel0(:,a)] = getKernel(disval, hdxmat(rtr0, begin:begin+frameperbin-1), ch(rtr0));  
        [tkernel1(:,a)] = getKernel(disval, hdxmat(rtr1, begin:begin+frameperbin-1), ch(rtr1));  
        begin = begin + frameperbin;
    end
    tapk = mean(tkernel,2)';
    if logreg_flag==0
        % image classification
        for a = 1:nbin
            ampr(r,a) = tapk*tkernel(:,a);
            ampr0(r,a) = tapk*tkernel0(:,a);
            ampr1(r,a) = tapk*tkernel1(:,a);
        end
    else
        % logistic regression
        ampr = glmfit(hdxmatbin(rtr,:), ch, ...
            'binomial', 'link', 'logit', 'constant', 'on');
        ampr0 = glmfit(hdxmatbin(rtr0, :), ch(rtr0), ...
            'binomial', 'link', 'logit', 'constant', 'on');
        ampr1 = glmfit(hdxmatbin(rtr1, :), ch(rtr1), ...
            'binomial', 'link', 'logit', 'constant', 'on');
        amp = amp(2:end); ampr0 = ampr0(2:end); ampr1 = ampr1(2:end);
    end        
end
err = std(ampr,[],1);
err0 = std(ampr0,[],1);
err1 = std(ampr1,[],1);

function adap = compute_adaptation(mpsth, offset, binsize, frameperbin)
psth = mpsth - mean(mpsth(1:offset));
adap = nan(1, binsize);
begin = 1+offset;
for b = 1:binsize
    adap(b) = mean(psth(begin:begin+frameperbin-1))/max(psth);
    begin = begin + frameperbin;
end

function [cp1, cp2] = compute_cp(psth, ch, binsize, frameperbin)
cp1 = nan(1, binsize+1);
cp2 = nan(1, binsize+1);
begin = 1;
cp1(1) = rocN(mean(psth(ch==1, begin:begin+binsize*frameperbin-1), 2), ...
        mean(psth(ch==0, begin:begin+binsize*frameperbin-1), 2));
[~,~,~,cp2(1)] = perfcurve(ch, mean(psth(:, begin:begin+binsize*frameperbin-1),2), 1);
for b = 1:binsize
    cp1(b+1) = rocN(mean(psth(ch==1, begin:begin+frameperbin-1), 2), ...
        mean(psth(ch==0, begin:begin+frameperbin-1), 2));
    [~,~,~,cp2(b+1)] = perfcurve(ch, mean(psth(:, begin:begin+frameperbin-1),2), 1);
    begin = begin + frameperbin;
end

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

function [cp] = rocN(x,y)
N = 100;
[m n] = size(x);
x = reshape(x,1,m*n);
[m n] = size(y);
y = reshape(y,1,m*n);

zlo = min([min(x(:)) min(y(:))]);
zhi = max([max(x(:)) max(y(:))]);
z = linspace(zlo,zhi,N);
fa = zeros(1,N);	% allocate the vector
hit = zeros(1,N);
for i = 1:N
  fa(N-i+1) = sum(y > z(i));
  hit(N-i+1) = sum(x > z(i));
end
[m,ny] = size(y);
fa = fa/ny;
[m,nx] = size(x);
hit = hit/nx;
fa(1) = 0;
hit(1) = 0;
fa(N) = 1;
hit(N) = 1;
cp = trapz(fa,hit);
