function [para] = SDTvariants(varargin)
%% simulation of SDT using specified "noise", "weights", and "db"
%% pre-set parameters
% hdx = [0.3 0.2 0.1 0.01 -0.01 -0.1 -0.2 -0.3];
hdx = 0.3*[-1 -0.75 -0.5 -0.25 0 0.25 0.5 0.75 1];
kernel = [0.1 0.5 0.2 0.05 0 -0.05 -0.2 -0.5 -0.1];
bias = 0;
% threshold = 30;
% sensitivity = 10/threshold;
len_tr = 5000;
len_frame = 100;

%% parameters of interest
% add Gaussian noise
noise = 0;

% weights on time-course
weights = ones(1, len_frame);

% decision boundry
db = 100;
    
% number of bin
nbin = 4;

% figure
plot_flag = 0;

% resampling
resample_flag = 0;

j = 1;              
while  j<= length(varargin)
    switch varargin{j}
        case 'noise'
            noise = varargin{j+1};
            j = j + 2;
        case 'weights' 
            weights = varargin{j+1};
            j = j + 2;
        case 'db'
            db = varargin{j+1};
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

disp('----------------------------------------------------------------')
disp(['noise: ' num2str(noise) ', weightes: ' num2str((weights(end)-weights(1))/len_frame) ', DB:' num2str(db)])

%% generate artificial dataset
% seed
rng(19891220);

% generate dynamic stimulus sequence with the 0% signal
stm = datasample(hdx, len_tr*len_frame, 'Replace',true);
stm = reshape(stm, [len_tr, len_frame]);

% generate DVs
% dvs = normrnd(0, sensitivity, [len_tr, len_frame]);
dvs = stm;
for s = 1:length(hdx)
    dvs(stm==hdx(s)) = kernel(s);
end

% add Gaussian noise
sensoryrep = dvs;
sensorynoise = normrnd(0, noise, [len_tr, len_frame]);
dvs = sensoryrep + sensorynoise;

% weights on time-course
dvs = dvs.*repmat(weights, len_tr, 1);

% choice based on decision boundry
advs = cumsum(dvs, 2);
db_advs = advs;
dt = len_frame*ones(len_tr,1);
dbreach = zeros(len_tr, 1);
for c = 1:len_tr
    idx = find(abs(advs(c,:)) >= db);
    if ~isempty(idx)
        db_advs(c, idx(1):end) = advs(c, idx(1));
        dt(c) = idx(1);
        dbreach(c) = 1;
    end
end

disp(['The % trials reaching the DB: ' num2str(100*sum(dbreach==1)/len_tr)])

% median split of DVs
if db==100
    conf = abs(db_advs(:,end));
else
    conf = (2/pi)*atan(abs(db_advs(:,end))./(dt/len_frame));
end

% add noise on confidence judgement
conf = conf + normrnd(median(conf), 0.2*(max(max(conf)) - min(min(conf))), size(conf));

med = median(conf);
disp([num2str(100*sum(dbreach==1 & conf > med)/sum(conf > med)) ...
    '% trials reached DB in high confidence, '...
    num2str(100*sum(dbreach==1 & conf < med)/sum(conf < med))...
    '% trials reached DB in low confidence '])

% choice 
ch = sign(db_advs(:,end) - (bias/100)*ones(len_tr,1));

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

% % time-averaged kernel
% avkernel = getKernel(stm, ch);
% avkernel_h = getKernel(stm(conf > med,:), ch(conf > med));
% avkernel_l = getKernel(stm(conf < med,:), ch(conf < med));

amp = nan(1,nbin);
amph = nan(1,nbin);
ampl = nan(1,nbin);
for a = 1:nbin
    % amplitude of the PK
    amp(a) = dot(tkernel(:,a), mean(tkernel,2));
    amph(a) = dot(tkernel_h(:,a), mean(tkernel_h,2));
    ampl(a) = dot(tkernel_l(:,a), mean(tkernel_l,2));
end

% output argumant
para = struct('trKernel', tkernel, 'trKernel_highconf', tkernel_h, 'trKernel_lowconf', tkernel_l, ...
    'amplitude', amp, 'amplitude_highconf', amph, 'amplitude_lowconf', ampl);

% visualization
if plot_flag==1
    % yellow and green
    y = [0.9576    0.7285    0.2285];
    g = [0.1059    0.4706    0.2157];
    
    close all;
    figure;
    
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
    cfix

    subplot(1,3,3)
    nom = mean([amph, ampl]);
    if resample_flag==1
        repeat = 500;
        errh = resamplePK(stm(conf > med,:), ch(conf > med), nbin, frameperbin, repeat);
        errl = resamplePK(stm(conf < med,:), ch(conf < med), nbin, frameperbin, repeat);
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

function err = resamplePK(hdxmat, ch, nbin, frameperbin, repeat)
disval = unique(hdxmat);
len_d = length(disval);
ampr = nan(repeat, nbin);
for r = 1:repeat
    rtr = datasample(1:length(ch), length(ch));
    tkernel = nan(len_d, nbin);
    begin = 1;
    for n = 1:nbin
        tkernel(:,n) = getKernel(hdxmat(rtr,begin:begin+frameperbin-1), ch(rtr));
        begin = begin + frameperbin;
    end
    for n = 1:nbin
        ampr(r,n) = dot(tkernel(:,n), mean(tkernel,2));
    end
end
err = std(ampr, [], 1);