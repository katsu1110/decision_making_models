function trst = dynamic_stm_generator(ntrsig, nframe, hdx, sig)
% generate dynamic stimulus sequence 

% deal with inputs
if nargin < 1; ntrsig = 100; end
if nargin < 2; nframe = 50; end
if nargin < 3; hdx = 0.3*[-1:0.25:1]; end
if nargin < 4
    sig = 0.5*[-1 -0.5 -0.25 -0.125 -0.0625 0 0.0625 0.125 0.25 0.5 1];
end

% seed
rng(19891220);
% hdx (stimulus feature)
nhdx = length(hdx);
% percent signal
nsig = length(sig);
% generate dynamic stimulus sequence
begin = 1;
pd = zeros(nsig, nhdx);
stm = zeros(ntrsig*nsig, nframe);
signal_hdx_idx = 2;
zeropos = find(hdx==0);
C = zeros(ntrsig*nsig, 1);
for n = 1:nsig
    % signal strength
    C(begin:begin+ntrsig-1) = sig(n);
    % probability distribution of hdx in each signal
    restper = 1 - abs(sig(n));
    pd(n,:) = (restper/nhdx)*ones(1,nhdx);
    sigpos = zeropos + sign(sig(n))*signal_hdx_idx;
    pd(n, sigpos) = abs(sig(n)) + pd(n, sigpos);
    % generate dynamic stimulus sequence in each signal
    stm_temp = datasample(hdx, ntrsig*nframe, 'Replace', true, 'Weights', pd(n,:));
    stm(begin:begin+ntrsig-1, :) = reshape(stm_temp, [ntrsig, nframe]);
    begin = begin + ntrsig;
end
% randomize order
rtr = randperm(size(stm,1));
C = C(rtr);
stm = stm(rtr, :);
% into the output
trst.hdx = hdx;
trst.sig = sig;
trst.hdx_prob = pd;
trst.C = C;
trst.stm = stm;