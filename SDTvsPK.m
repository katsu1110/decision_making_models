function SDTvsPK
%% simulation to test any relationships between SDT-based DV and psychophysical kernel
% INPUT: noise ... noise to encode decision variable: "noise", positive float (default: 0.1)
%             weights ... weights on time-course: "weights", vector 1 x 100 (default: ones(1, 100))
%             db ... decision boundary: "db", positive float (default: 10)
% 
%
% TEST: 
% - reduce pooling noise
% - increase weights
% - increase decision bound
%
%
% written by Katsuhisa (21.08.17)
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% yellow and green
y = [0.9576    0.7285    0.2285];
g = [0.1059    0.4706    0.2157];

s = 5;
noise = [2, 1.5, 1, 0.5, 0];
db = [1,2,3,4,100];
weights = ones(s, 100);
rng(1220);
weights(1,:) = sort(linspace(0,2,100), 'descend');
weights(2,:) = sort(linspace(0.3, 1.7, 100), 'descend');
weights(3,:) = sort(linspace(0.6, 1.4, 100), 'descend');
weights(4,:) = sort(linspace(0.9, 1.1, 100), 'descend');

close all;

for t = 1:3
    h = figure;
    p = zeros(1,s);
    for i = 1:s
        switch t
            case 1
                noise_i = noise(i);
                weights_i = weights(s,:);
                db_i = db(end);
                label = 'pooling noise';
                legendtext = {'2', '1.5', '1', '0.5', '0'};
            case 2
                noise_i = noise(end);
                weights_i = weights(i,:);
                db_i = db(end);
                label = 'sensory weights';
                legendtext = {'descend1', 'descend2','descend3','descend4','flat'};
            case 3
                noise_i = noise(end);
                weights_i = weights(s,:);
                db_i = db(i);
                label = 'decision boundary';
                legendtext = {'1','2','3','4','100'};
        end
                
        [bhvSt, kernelSt, sdtsigSt] = SDTsimulator('noise', noise_i, 'weights', weights_i, 'db', db_i);        
    
        % visualization =============================
        % median split by DV
        col = colormap(lines(s));
        subplot(5,5,2 + (i-1)*5)
        plot(1:4, kernelSt.amp_l, ':', 'color', col(i,:),'linewidth',2);
        hold on;
        plot(1:4, kernelSt.amp_h, '--', 'color', col(i,:),'linewidth',2);
        xlim([0.5 4.5])
%         yy = get(gca, 'YLim');
%         ylim([0 yy(2)])
%         title(label)
        set(gca,'XTick',[1 2 3 4])
        set(gca, 'box', 'off')
        set(gca, 'TickDir', 'out')
    %     axis square

        % original
        subplot(5,5,1 + (i-1)*5)
        p(i) = plot(1:4, kernelSt.amp, '-', 'color', col(i,:),'linewidth',2);
        xlim([0.5 4.5])
%         yy = get(gca, 'YLim');
%         ylim([0 yy(2)])
    %     if t==2
    %         xlabel('time bin')
    %     end
        if t==1
            ylabel('amplitude of PK')
        end
        title(legendtext{i})
        set(gca,'XTick',[1 2 3 4])
        set(gca, 'box', 'off')
        set(gca, 'TickDir', 'out')
    %     axis square

        % accuracy vs confidence
        subplot(5,5,3+(i-1)*5)
%         b = sdtsigSt.confVSacc;
%         x = linspace(min(bhvSt.conf), max(bhvSt.conf), 100);
%         plot(x, x*b(2) + b(1), '-k')
        plot(1:sdtsigSt.confVSacc.nbin, sdtsigSt.confVSacc.bin, '-k', 'linewidth',1)
%         xlim([min(x)-1 max(x)+1])
%         yy = get(gca, 'YLim');
%         ylim([49 yy(2)])
        set(gca, 'box', 'off'); set(gca, 'TickDir', 'out');
        xlabel('confidence bin')
        ylabel('accuracy')

        % psychometric functions separated by confidence
        subplot(5,5,4+(i-1)*5)
        for k = 1:2
            switch k
                case 1
                    col = g;
                case 2
                    col = y;
            end
            plot(sdtsigSt.pmVSconf(k).sig, sdtsigSt.pmVSconf(k).hit,...
                '-', 'color', col)
            hold on;
        end
        xlim([min(sdtsigSt.pmVSconf(k).sig)-0.001, ...
            max(sdtsigSt.pmVSconf(k).sig)+0.001])
        yy = get(gca, 'YLim');
        ylim([49 yy(2)])
        set(gca, 'box', 'off'); set(gca, 'TickDir', 'out');
        xlabel('signal')
        ylabel('accuracy')

        % confidence vs signal vs accuracy
        subplot(5,5,5+(i-1)*5)
        for k = 1:2
            switch k
                case 1
                    col = 0.6*ones(1,3);
                case 2
                    col = [0 0 0];
            end
            errorbar(sdtsigSt.confVSsig(k).sig, sdtsigSt.confVSsig(k).conf_mean, ...
                sdtsigSt.confVSsig(k).conf_sem, '-', 'color', col, 'CapSize', 0)
            hold on;
        end
        xlim([min(sdtsigSt.confVSsig(k).sig)-0.001, max(sdtsigSt.confVSsig(k).sig)+0.001])
        set(gca, 'box', 'off'); set(gca, 'TickDir', 'out');
        xlabel('signal')
        ylabel('confidence')    
    end
    
    set(h, 'Name', ['manipulated ' label], 'NumberTitle', 'off')
end        


function [bhvSt, kernelSt, sdtsigSt] = SDTsimulator(varargin)
%% simulation of SDT using specified "noise", "weights", and "db"
%% pre-set parameters
hdx = [0.3 0.2 0.1 0.01 -0.01 -0.1 -0.2 -0.3];
kernel = [0.2 0.5 0.2 0.05 -0.05 -0.2 -0.5 -0.2];
signal = -0.2:0.02:0.2;
bias = 0;
% threshold = 30;
% sensitivity = 10/threshold;
len_tr = 10000;
len_frame = 100;
len_sig = length(signal);
len_hdx = length(hdx);

%% parameters of interest
% add Gaussian noise
noise = 0.1;

% weights on time-course
weights = ones(1, len_frame);

% decision boundry
db = 100;
    
j = 1;              
while  j<= length(varargin)
    switch varargin{j}
        case 'noise'
            noise = varargin{j+1};
        case 'weights' 
            weights = varargin{j+1};
        case 'db'
            db = varargin{j+1};
    end
    j = j + 2;
end

disp('----------------------------------------------------------------')
disp(['noise: ' num2str(noise) ', weightes: ' num2str((weights(end)-weights(1))/len_frame) ', DB:' num2str(db)])

%% generate artificial dataset
% % seed
% rng(19891220);

% hdx-matrix in each signal
dvmat = nan(len_tr*len_sig, len_frame);
dt = len_frame*ones(1, len_tr*len_sig);
sig = zeros(1, len_tr*len_sig);
ch = nan(1, len_tr*len_sig);
acc = zeros(1, len_tr*len_sig);
conf = zeros(1, len_tr*len_sig);
confi = zeros(1, len_tr*len_sig);
hdxmat = nan(len_tr*len_sig, len_frame);
begin = 1;
for s = 1:len_sig
    % signal
    sig(begin:begin+len_tr-1) = signal(s)*ones(1, len_tr);
    
    % generate dynamic stimulus sequence
    weightvec = (1 - abs(signal(s)))*ones(1, len_hdx)/len_hdx;
    if signal(s) > 0
        weightvec(2) = (1 - abs(signal(s)))/len_hdx + abs(signal(s));
    elseif signal(s) < 0
        weightvec(end-1) = (1 - abs(signal(s)))/len_hdx + abs(signal(s));
    end
    hdxseq = datasample(hdx, len_tr*len_frame, 'Replace',true, 'Weights', weightvec);
    hdxseq = reshape(hdxseq, [len_tr, len_frame]);
    hdxmat(begin:begin+len_tr-1,:) = hdxseq;
    
    % generate DVs
    % dvs = normrnd(0, sensitivity, [len_tr, len_frame]);
    dv_raw = hdxseq;
    for h = 1:len_hdx
        dv_raw(hdxseq==hdx(h)) = kernel(h);
    end

    % add Gaussian noise
    sensoryrep = dv_raw;
    sensorynoise = normrnd(0, noise, [len_tr, len_frame]);
    dv_raw = sensoryrep + sensorynoise;
    
    % weights on time-course
    dv_raw = dv_raw.*repmat(weights, len_tr, 1);

    % choice based on decision boundry
    dv_acc = cumsum(dv_raw, 2);
    dbreach = zeros(len_tr, 1);
    for c = 1:len_tr
        idx = find(abs(dv_acc(c,:)) >= db);
        if ~isempty(idx)
            dv_acc(c, idx(1):end) = dv_acc(c, idx(1));
            dt(begin+c-1) = idx(1);
            dbreach(c) = 1;
        end
    end
    dvmat(begin:begin+len_tr-1, :) = dv_acc;

    % median split of DVs
    confs = abs(dv_acc(:,end));
    med = median(confs);
    conf(begin:begin+len_tr-1) = confs;
    confi(begin:begin+len_tr-1) = 1*(confs > med);

    % % median split of |DV|/dt
    % conf = (2/pi)*(atan(abs(db_advs(:,end))./(dt/len_frame)));

    % choice 
    chs = sign(dv_acc(:,end) - (bias/100)*ones(len_tr,1));

    idx0 = find(chs==0);
    len0 = length(idx0);
    if len0 > 0
        chs(idx0) = datasample([-1 1], len0, 'Replace', true);
    end
    
    % accuracy
    if signal(s)==0
        accs = datasample([0 1], len_tr, 'Replace', true);
    else
        accs  = zeros(1, len_tr);
        accs(sign(chs')==sign(signal(s))) = 1;
    end
    acc(begin:begin+len_tr-1) = accs;
    
    chs(chs==-1) = 0;
    ch(begin:begin+len_tr-1) = chs;
    
    begin = begin + len_tr;
    
    % display some info        
    disp(['stimulus: ' num2str(signal(s))])
    disp(['accuracy: ' num2str(100*sum(accs==1)/len_tr) ' %'])
    disp(['The % trials reaching the DB: ' num2str(100*sum(dbreach==1)/len_tr)])

    disp([num2str(100*sum(dbreach==1 & confs > med)/sum(confs > med)) ...
        '% trials reached DB in high confidence, '...
        num2str(100*sum(dbreach==1 & confs < med)/sum(confs < med))...
        '% trials reached DB in low confidence '])

    disp([num2str(sum(chs==0)) ' near-choices, ' num2str(sum(chs==1)) ' far-choices'])
    disp('----------------------------------------------------------------')
    
    if signal(s)==0
        % amplitude of the psychophysical kernel as a function of time
        tkernel = nan(len_hdx, 4);
        tkernel_h = nan(len_hdx, 4);
        tkernel_l = nan(len_hdx, 4);
        for a = 1:4
            % time-binned PK
            pk0 = getKernel(hdxseq, chs);
            pkh = getKernel(hdxseq(confs > med,:), chs(confs > med));
            pkl = getKernel(hdxseq(confs < med,:), chs(confs < med));
            tkernel(:,a) = pk0';
            tkernel_h(:,a) = pkh';
            tkernel_l(:,a) = pkl';
        end

        amp = nan(1,4);
        amph = nan(1,4);
        ampl = nan(1,4);
        for a = 1:4
            % amplitude of the PK
            amp(a) = dot(tkernel(:,a), mean(tkernel,2)');
            amph(a) = dot(tkernel_h(:,a), mean(tkernel_h,2)');
            ampl(a) = dot(tkernel_l(:,a), mean(tkernel_l,2)');
        end
        
        % kernel data
        kernelSt = struct('tKernel', tkernel, 'tKernel_h', tkernel_h, ...
            'tKernel_l', tkernel_l, 'amp', amp, 'amp_h', amph, 'amp_l', ampl);
    end  
end

% % normalize confidence
% conf = (conf - min(conf))/(max(conf) - min(conf));

% behavioral data
bhvSt = struct('dv', dvmat, 'dt', dt, 'sig', sig, 'ch', ch, ...
    'conf', conf,'conf_i', confi);

% signatures of confidence based on SDT
% confidence vs accuracy
b = glmfit(conf', acc', 'binomial', 'link', 'logit', 'constant', 'on');
sdtsigSt.confVSacc.b = b;
nbin = 20;
sdtsigSt.confVSacc.nbin = nbin;
sdtsigSt.confVSacc.bin = nan(1, nbin);
trperbin = len_tr/nbin;
[~, idx] = sort(conf);
a = 1;
for n = 1:nbin
    sdtsigSt.confVSacc.bin(n) = 100*sum(acc(idx(a:a+trperbin-1))==1)/trperbin;
    a = a + nbin;
end

% psychometric functions by confidence levels
[hitl, sigl] = getPM(sig(confi==0), acc(confi==0));
[hith, sigh] = getPM(sig(confi==1), acc(confi==1));
sdtsigSt.pmVSconf(1).sig = sigl;
sdtsigSt.pmVSconf(2).sig = sigh;
sdtsigSt.pmVSconf(1).hit = hitl;
sdtsigSt.pmVSconf(2).hit = hith;

% confidence vs signal vs correctness
[cf_me_l, cf_sem_l, uni_sig_l] = getCFSig(sig(acc==0), conf(acc==0));
[cf_me_h, cf_sem_h, uni_sig_h] = getCFSig(sig(acc==1), conf(acc==1));
sdtsigSt.confVSsig(1).sig = uni_sig_l;
sdtsigSt.confVSsig(2).sig = uni_sig_h;
sdtsigSt.confVSsig(1).conf_mean = cf_me_l;
sdtsigSt.confVSsig(2).conf_mean = cf_me_h;
sdtsigSt.confVSsig(1).conf_sem = cf_sem_l;
sdtsigSt.confVSsig(2).conf_sem = cf_sem_h;

%% subfunctions
% psychophysical kernel
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
pk0 = mean(svmat(ch==0,:)) - mean(svmat(ch==1,:));

% psychometric function
function [hit, uni_sig] = getPM(sig, acc)
sig = abs(sig);
uni_sig = unique(sig);
len_unis = length(uni_sig);
hit = zeros(1, len_unis);
for s = 1:len_unis
    hit(s) = 100*sum(acc==1 & sig==uni_sig(s))/sum(sig==uni_sig(s));
end

% confidence vs signals
function [cf_me, cf_sem, uni_sig] = getCFSig(sig, conf)
sig = abs(sig);
uni_sig = unique(sig);
len_unis = length(uni_sig);
cf_me = zeros(1, len_unis);
cf_sem = zeros(1, len_unis);
for s = 1:len_unis
    cf_me(s) = mean(conf(sig==uni_sig(s)));
    cf_sem(s) = std(conf(sig==uni_sig(s)))/sqrt(sum(sig==uni_sig(s)));
end

    