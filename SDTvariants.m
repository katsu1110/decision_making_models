function [para] = SDTvariants(varargin)
%% 
% simulation of Signal detection theory (SDT) manipulating specified "noise", "weights", "db" and "leak"
% INPUT:
% 'noise' ... float; add pooling noise 
% 'weights' ... vector with its length of len_frame (100 in default); multiplied with stimulus sequence
% 'db' ... float; decision boundary
% 'leak' ... float; leaky integration
% 'nbin' ... the number of time bins to compute the time-resolved
% psychophysical kernel
% 'resample' ... perform resampling (bootstrap) to get error bars for the
% kernel amplitude
% 'plot' ... plot the results
%
% OUTPUT:
% matlab structure including relevant information



%% pre-set parameters
% hdx = [0.3 0.2 0.1 0.01 -0.01 -0.1 -0.2 -0.3];
hdx = 0.3*[-1 -0.75 -0.5 -0.25 0 0.25 0.5 0.75 1];
kernel = [0.1 0.5 0.2 0.05 0 -0.05 -0.2 -0.5 -0.1];
bias = 0;
% threshold = 30;
% sensitivity = 10/threshold;
len_tr = 5000;
len_frame = 100;
stmstrength = 0; % 0 or 0.12 <= x < 0.5

%% parameters of interest
% add Gaussian noise
noise = 0;

% weights on time-course
weights = ones(1, len_frame);

% decision boundry
db = 100;

% taking into acount decision time
dt_flag = 0;
dtweight = 1;

% leak 
leak = 0;

% number of bin
nbin = 4;

% figure
plot_flag = 0;

% resampling
resample_flag = 0;

j = 1;              
while  j<= length(varargin)
    switch varargin{j}
        case 'ntr'
            len_tr = varargin{j+1};
            j = j + 2;
        case 'stmstrength'
            stmstrength = varargin{j+1};
            j = j + 2;
        case 'noise'
            noise = varargin{j+1};
            j = j + 2;
        case 'weights' 
            weights = varargin{j+1};
            j = j + 2;
        case 'db'
            db = varargin{j+1};
            j = j + 2;
        case 'leak'
            leak = varargin{j+1};
            j = j + 2;
        case 'dt'
            dt_flag = 1;
            j = j + 1;
        case 'dtw'
            dtweight = varargin{j+1};
            j = j + 2;
        case 'nbin'
            nbin = varargin{j+1};            
             j = j + 2;
        case 'resample'
            resample_flag = 1;
            j = j + 1;
        case 'plot'
            plot_flag = 1;
            j = j + 1;
    end
end

pd = ones(length(hdx),1);
% % both sides
% if stmstrength>0  
%     pd = pd*(1 - 2*stmstrength)/7;
%     pd(2) = stmstrength;
%     pd(8) = stmstrength;
% end 
% one side
if stmstrength>0  
    pd = pd*(1 - stmstrength)/8;
    pd(2) = stmstrength;
end 

disp('----------------------------------------------------------------')
disp(['noise: ' num2str(noise) ', weightes: ' num2str((weights(end)-weights(1))/len_frame) ', DB:' num2str(db) ', leak:' num2str(leak)])

leak = leak./(linspace(1,len_frame,len_frame)/len_frame);

%% generate artificial dataset
% seed
rng(19891220);

% generate dynamic stimulus sequence with the 0% signal
stm = datasample(hdx, len_tr*len_frame, 'Replace', true, 'Weights', pd);
stm = reshape(stm, [len_tr, len_frame]);

% generate DVs
% dvs = normrnd(0, sensitivity, [len_tr, len_frame]);
dvs = stm;
for s = 1:length(hdx)
    dvs(stm==hdx(s)) = kernel(s);
end

% noise level in the stimulus
noisestm = std(stm(:));
noiseidv = std(dvs(:));

% add Gaussian noise
sensoryrep = dvs;
sensorynoise = normrnd(0, noise, [len_tr, len_frame]);
dvs = sensoryrep + sensorynoise;

% choice based on decision boundry
idvs = zeros(len_tr, len_frame);
dt = len_frame*ones(len_tr,1);
dbreach = zeros(len_tr, 1);
for c = 1:len_tr
    
    % sensory weighting
    idvs(c,:) = dvs(c,:).*weights;
    
    % (leaky) integration
    for f = 2:len_frame        
        idvs(c,f) = idvs(c,f-1)*(1 - leak(f)) + dvs(c,f);
    end
    
    % bound and decision time
    idx = find(abs(idvs(c,:)) >= db);
    if ~isempty(idx)
        idvs(c, idx(1):end) = idvs(c, idx(1));
        dt(c) = idx(1);
        dbreach(c) = 1;
    end
end

nreach0 = 100*sum(dbreach==1)/len_tr;
disp(['The % trials reaching the DB: ' num2str(nreach0)])

% median split of DVs
if dt_flag==1
    conf = (2/pi)*atan(abs(idvs(:,end))./(dtweight*dt/len_frame));
else
    conf = abs(idvs(:,end));
end

% % add noise on confidence judgement
% conf = conf + normrnd(median(conf), 0.2*median(conf), size(conf));

med = median(conf);
nreach1 = 100*sum(dbreach==1 & conf < med)/sum(conf < med);
nreach2 = 100*sum(dbreach==1 & conf > med)/sum(conf > med);
disp([num2str(nreach2) ...
    '% trials reached DB in high confidence, '...
    num2str(nreach1)...
    '% trials reached DB in low confidence '])

% choice 
ch = sign(idvs(:,end) - (bias/100)*ones(len_tr,1));

idx0 = find(ch==0);
len0 = length(idx0);
if len0 > 0
    ch(idx0) = datasample([-1 1], len0, 'Replace', true);
end
ch(ch==-1) = 0;
disp([num2str(sum(ch==0)) ' near-choices, ' num2str(sum(ch==1)) ' far-choices'])
disp('----------------------------------------------------------------')

% amplitude of the psychophysical kernel as a function of time
tkernel = nan(length(hdx), nbin);
tkernel_h = nan(length(hdx), nbin);
tkernel_l = nan(length(hdx), nbin);
begin = 1;
frameperbin = floor(len_frame/nbin);
for a = 1:nbin
    % time-binned PK
    pk0 = getKernel(stm(:,begin:begin+frameperbin-1), ch);
    pkh = getKernel(stm(conf > med, begin:begin+frameperbin-1), ch(conf > med));
    pkl = getKernel(stm(conf < med, begin:begin+frameperbin-1), ch(conf < med));
    tkernel(:,a) = pk0';
    tkernel_h(:,a) = pkh';
    tkernel_l(:,a) = pkl';
    begin = begin + frameperbin;
end

% time-averaged kernel
avkernel = getKernel(stm, ch);
% avkernel_h = getKernel(stm(conf > med,:), ch(conf > med));
% avkernel_l = getKernel(stm(conf < med,:), ch(conf < med));

amp = nan(1,nbin);
amph = nan(1,nbin);
ampl = nan(1,nbin);
for a = 1:nbin
    % amplitude of the PK
%     amp(a) = dot(tkernel(:,a), mean(tkernel,2));
%     amph(a) = dot(tkernel_h(:,a), mean(tkernel_h,2));
%     ampl(a) = dot(tkernel_l(:,a), mean(tkernel_l,2));

    amp(a) = dot(tkernel(:,a), avkernel);
    amph(a) = dot(tkernel_h(:,a), avkernel);
    ampl(a) = dot(tkernel_l(:,a), avkernel);
end
if mean(isnan(amph))
    amph = 2*amp - ampl;
elseif mean(isnan(ampl))
    ampl = 2*amp - amph;
end

% output argumant
para = struct('trKernel', tkernel, 'trKernel_highconf', tkernel_h, 'trKernel_lowconf', tkernel_l, ...
    'amplitude', amp, 'amplitude_highconf', amph, 'amplitude_lowconf', ampl,...
    'nreach',nreach0,'nreach_highconf',nreach2,'nreach_lowconf',nreach1,...
    'noisestm',noisestm,'noiseidv',noiseidv);

% visualization
if plot_flag==1
    % yellow and green
    y = [0.9576    0.7285    0.2285];
    g = [0.1059    0.4706    0.2157];
    
    close all;
    figure;
    subplot(1,2,1)
    imagesc(1:nbin, hdx, tkernel/max(tkernel(:)))
    colormap(copper)
    xlabel('time bin')
    ylabel('hdx')
    set(gca, 'XTick', 1:nbin)
    set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
    c = caxis;
    set(gca, 'XTick', 1:nbin)
    set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
    title(c)

    subplot(1,2,2)
    nom = mean(amp);
    if resample_flag==1
        repeat = 1000;
        [err, errl, errh] = resamplePK(stm, ch, nbin, frameperbin, conf, med, repeat);
        fill_between(1:nbin, (amp - err)/nom, (amp + err)/nom, [0 0 0])
        hold on;       
        plot(1:nbin, amp/nom, '-', 'color', [0 0 0], 'linewidth', 2)
    else
        plot(1:nbin, amp/nom, '-', 'color', [0 0 0], 'linewidth', 2)
    end    
    xlim([0.5 nbin + 0.5])
    ylabel('kernel amplitude')
    set(gca, 'XTick', 1:nbin)
    set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
    
    figure;
    
    subplot(1,3,1)
    vmax = max([max(tkernel_h(:)), max(tkernel_l(:))]);
    imagesc(1:nbin, hdx, tkernel_h/vmax)
    colormap(copper)
    c_h = caxis;
    xlabel('time bin')
    ylabel('hdx')
    set(gca, 'XTick', 1:nbin)
    set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')

    subplot(1,3,2)
    imagesc(1:nbin, hdx, tkernel_l/vmax)
    colormap(copper)
    c_l = caxis;
    set(gca, 'XTick', 1:nbin)
    set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')

    cfix = [min([c_h, c_l]), max([c_h, c_l])];
    subplot(1,3,1)
    caxis(cfix)
    subplot(1,3,2)
    caxis(cfix)
    title(cfix)

    subplot(1,3,3)
%     nom = mean([amph, ampl]);
    nom = max(amp);
    if resample_flag==1
        fill_between(1:nbin, (ampl - errl)/nom, (ampl + errl)/nom, g)
        hold on;
        plot(1:nbin, ampl/nom, '-', 'color', g, 'linewidth', 2)
        hold on;
        fill_between(1:nbin, (amph - errh)/nom, (amph + errh)/nom, y)
        hold on;       
        plot(1:nbin, amph/nom, '-', 'color', y, 'linewidth', 2)
    else
        plot(1:nbin, ampl/nom, '-', 'color', g, 'linewidth', 2)
        hold on;
        plot(1:nbin, amph/nom, '-', 'color', y, 'linewidth', 2)
    end    
    xlim([0.5 nbin + 0.5])
    ylabel('kernel amplitude')
    set(gca, 'XTick', 1:nbin)
    set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
end
    

%% subfunctions
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
pk0 = mean(svmat(ch==1,:)) - mean(svmat(ch==0,:));

function [err, err0, err1] = resamplePK(hdxmat, ch, nbin, frameperbin, conf, med, repeat)
conf0 = find(conf < med);
conf1 = find(conf > med);
disval = unique(hdxmat);
len_d = length(disval);
ampr = nan(repeat, nbin);
ampr0 = nan(repeat, nbin);
ampr1 = nan(repeat, nbin);
for r = 1:repeat
    rtr0 = datasample(conf0, length(conf0))';
    rtr1 = datasample(conf1, length(conf1))';
    rtr = [rtr0, rtr1];
    tkernel = nan(len_d, nbin);
    tkernel0 = nan(len_d, nbin);
    tkernel1 = nan(len_d, nbin);
    begin = 1;
    for n = 1:nbin
        tkernel(:,n) = getKernel(hdxmat(rtr,begin:begin+frameperbin-1), ch(rtr));
        tkernel0(:,n) = getKernel(hdxmat(rtr0,begin:begin+frameperbin-1), ch(rtr0));
        tkernel1(:,n) = getKernel(hdxmat(rtr1,begin:begin+frameperbin-1), ch(rtr1));
        begin = begin + frameperbin;
    end
    tapk = getKernel(hdxmat(rtr,:), ch(rtr));
    for n = 1:nbin
        ampr(r,n) = dot(tkernel(:,n), tapk);
        ampr0(r,n) = dot(tkernel0(:,n), tapk);
        ampr1(r,n) = dot(tkernel1(:,n), tapk);
    end
end
err = std(ampr, [], 1);
err0 = std(ampr0, [], 1);
err1 = std(ampr1, [], 1);

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
