function para = Yates_simulation(varargin)
%%
% simulations of "stimulus-to-SensoryNeuron model' presented in Yates et al., 2017
%
% INPUT:
% 'ntr' ... int; the number of trials
% 'nneuron' ... int; the number of sensory neurons per pool
% 'tmax' ... float; time length for the stm/co kernel
% 'tau' ... float; width of the stm/co kernel
% 'kernelgain_s' ... float; magnitude of the stm kernel
% 'kernelgain_c' ... float; magnitude of the co kernel
% 'offset_gain' ... float; magnitude of the negative peak of the stm/co kernel
% 'history' ... 0; without history kernel, 1; with history kernel
% 'stm_gain'... float; range of the stimulus
% 'nbin'... int; the number of time bin for PKA
% 'cfnoise' ... float; noise on confidence judgement. 0 is default. 
% 'pkmethod' ... method to compute psychophysical kernel amplitude (PKA)
%                0; Nienborg & Cumming, 2009, 1; image classification ,2; logistic regression
% 'repeat' ... the number of repeats for resampling
% 'plot'... 0; no figure, 1; figure
%  
% OUTPUT: output metrics storing simulated results
%
% EXAMPLE: Yates_simulation('plot')
%
% +++++++++++++++++++++++++++++++++++++++++++++++

% default input arguments
nneuron = 1;
len_tr = 5000;
tmax = 300; 
tau = 50; 
kernelgain_s = 0.05; 
kernelgain_c = 0.025; 
offset_gain = 1.5; 
stm_gain = 0.3;
cfnoise = 0;
history = 0;
repeat = 0;
nbin = 4;
pkmethod = 0;
plot_flag = 0;
j = 1;              
while  j <= length(varargin)
    switch varargin{j}
        case 'ntr'
            len_tr = varargin{j+1};
            j = j + 2;
        case 'nneuron'
            nneuron = varargin{j+1};
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
        case 'cfnoise'
            cfnoise = varargin{j+1};
            j = j + 2;
        case 'history'
            history = 1;               
             j = j + 1;
        case 'nbin'
            nbin = varargin{j+1};
            j = j + 2;
        case 'pkmethod'
            pkmethod = varargin{j+1};
            j = j + 2;
        case 'repeat'
            repeat = varargin{j+1};
            j = j + 2;
        case 'plot'
            plot_flag = 1;
            j = j + 1;
    end
end

% the number of frames
len_frame = nbin*150;

% discrete stimuli
hdx = stm_gain*[-1 -0.75 -0.5 -0.25 0 0.25 0.5 0.75 1];

% contrast
offset = 100; % for PSTH
co = 1*[zeros(1,offset),ones(1,len_frame),zeros(1,offset)];

% time
time = [1:len_frame+2*offset] - offset;    

%%
% alpha function as a stm/co kernel
t = 0:tmax;
kernel1 = exp(-t/tau).*(1 - exp(-t/tau));
ker_offset_s = [ones(1,round(t(end)/3))*0.025, 0.025:-0.025/(round(t(end)*2/3)):0]; 
ker_offset_c = [ones(1,length(t))*0.025]; 
kernel2 = kernel1 - offset_gain*ker_offset_s;
kernel2(1:4) = 0;
kernel1 = kernel1 - offset_gain*ker_offset_c;
kernel1(1:4) = 0;

% normalize kernels
para.kernel_co = kernelgain_c*kernel1/max(kernel1) ; % contrast kernel
para.kernel_stm = kernelgain_s*kernel2/max(kernel2);  % stimulus kernel    

%%
% kernel for the history term
if history
    ht = 0:10;
    kernel3 = log(1+ht);
    para.kernel_ht = normalize(kernel3, -0.018, 0);
end

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
c = conv(para.kernel_co, co);
for n = 1:nneuron
    for i = 1:len_tr
        % convolution with stimulus
        s1 = conv(para.kernel_stm, stm(i,:));
        s2 = conv(-para.kernel_stm, stm(i,:));

        if history
            % poisson process to generate spikes
            para.neuron(n).spk1(i,1) = poissrnd(exp(s1(1) + c(1)));
            para.neuron(n).spk2(i,1) = poissrnd(exp(s2(1) + c(1)));

            % with history term
            for f = 2:lenv
                % cumulative history terms
                h1 = fliplr(conv(para.kernel_ht, para.neuron(n).spk1(i,f-1)));
                h2 = fliplr(conv(para.kernel_ht, para.neuron(n).spk2(i,f-1)));
                if f > ht(end)+1
                    withspk1 = para.neuron(n).spk1(i,f-1-ht(end):f-1)>0;
                    withspk2 = para.neuron(n).spk2(i,f-1-ht(end):f-1)>0;
                else
                    withspk1 = zeros(1, ht(end)+1);
                    withspk2 = zeros(1, ht(end)+1);
                    withspk1(end-f+2:end) = para.neuron(n).spk1(i,1:f-1)>0;
                    withspk2(end-f+2:end) = para.neuron(n).spk2(i,1:f-1)>0;
                end
                para.neuron(n).spk1(i,f) = poissrnd(exp(s1(f) + c(f) + sum(h1(withspk1==1))),nneuron,1);
                para.neuron(n).spk2(i,f) = poissrnd(exp(s2(f) + c(f) + sum(h2(withspk2==1))),nneuron,1);
            end
        else            
            para.neuron(n).spk1(i,:) = exp(s1(1:lenv)+c(1:lenv));
            para.neuron(n).spk2(i,:) = exp(s2(1:lenv)+c(1:lenv));
        end
    end
end
disp('sensory responses generated')

%%
% evidence (decision variable)
dv1 = zeros(len_tr, len_frame);
dv2 = zeros(len_tr, len_frame);
for n = 1:nneuron
    dv1 = dv1 + cumsum(para.neuron(n).spk1(:, offset+1:offset+len_frame),2);
    dv2 = dv2 + cumsum(para.neuron(n).spk2(:, offset+1:offset+len_frame),2);
end
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

% add noise on confidence judgement
conf = conf + normrnd(0, cfnoise*std(conf), size(conf));

%%
% PSTH
k = randi(nneuron);
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

%%
% PKA
[amp, amph, ampl] = getPKA(zeros(len_tr,1), stm(:, offset+1:offset+len_frame),...
    ch, conf, nbin, pkmethod, repeat);
disp('PK computed')

%% 
% store variables
para.kernel_stm = kernel2;
para.kernel_co = kernel1;
para.stm = stm;
para.ch = ch;
para.confidence = conf;
para.dv = ev;
para.psth.stm_diff = psth_stm;
para.psth.ch_diff = psth_ch;
para.psth.stm_null = mean(para.neuron(k).spk1(stmsign_n2, 1:length(time)),1);
para.psth.stm_pref = mean(para.neuron(k).spk1(stmsign_p2, 1:length(time)),1);
para.psth.ch_pref = mean(para.neuron(k).spk1(ch==1, 1:length(time)),1);
para.psth.ch_null = mean(para.neuron(k).spk1(ch==0, 1:length(time)),1);
para.psth.stm_null_adaptation = compute_adaptation(para.psth.stm_null, offset, nbin, frameperbin);
para.psth.stm_pref_adaptation = compute_adaptation(para.psth.stm_pref, offset, nbin, frameperbin);
[cp1, cp2] = compute_cp(para.neuron(k).spk1(:, offset+1:offset+len_frame), ch, nbin, frameperbin);
para.psth.cp = [cp1; cp2];
para.pkmethod = pkmethod;
para.pka = amp(1,:);
para.pka_highconf = amph(1,:);
para.pka_lowconf = ampl(1,:);
if repeat > 0    
    para.pka_err = amp(2,:);
    para.pka_highconf_err = amph(2,:);
    para.pka_lowconf_err = ampl(2,:);
end

%%
% visualization
if plot_flag
    plot_yates(para);
end

% subfunction
function [stm, stmmat] = Yates_stm(hdx, len_tr, nbin, frameperbin, seed)
% generate stimulus 
rng(seed)
stmidx = randi(length(hdx),len_tr, nbin);
stmmat = hdx(stmidx);
stm = reshape(repmat(stmmat,frameperbin,1), len_tr, nbin*frameperbin);

function adap = compute_adaptation(mpsth, offset, binsize, frameperbin)
% compute adaptation level in each time bin 
psth = mpsth - mean(mpsth(1:offset));
adap = nan(1, binsize);
begin = 1+offset;
for b = 1:binsize
    adap(b) = mean(psth(begin:begin+frameperbin-1))/max(psth);
    begin = begin + frameperbin;
end

function [normalized_vector] = normalize(v, newmin, newmax)
% linear vectorial nomalization
a = (newmax - newmin)/(max(max(v)) - min(min(v)));
b = newmax - max(max(v))*a;
normalized_vector = a.*v + b;

function [cp1, cp2] = compute_cp(psth, ch, binsize, frameperbin)
% compute CP based on PSTHs from two pools of sensory neurons
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

function [cp] = rocN(x,y)
% HN's function to compute CP
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