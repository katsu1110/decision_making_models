function beh = weight_on_stm(trst, pk, pka, threshold)
% assign time-weight (pka) and dv-weight (pk) on the stimulus sequence and
% determine the choice
% 
% PK (Psychophysical kernel) ... stm to decision variable (template)
% PKA (Psychophysical Kernel Amplitude) ... weight on time
% threshold ... psychophysical threshold to determine the noise level
%

% stimulus
[ntr, nframe] = size(trst.stm);
idv = trst.stm;

% assign PK weight
lenhdx = length(trst.hdx);
idv_orig = idv;
for n = 1:lenhdx
    idv(idv_orig == trst.hdx(n)) = pk(n);
end

% noise on stm to dv
idv = idv + normrnd(0, threshold, ntr, nframe);

% weights on time
idv = idv.*repmat(pka, ntr, 1);

% integrate instantaneous decision variable
dv = cumsum(idv, 2);

% choice
ch = sign(dv(:,end));
ch(ch==0) = datasample([-1 1], sum(ch==0));
ch(ch==-1) = 0;

% accuracy
acc = datasample([0 1], ntr)';
acc(sign(ch) == sign(trst.C)) = 1;

% into output
beh.idv = idv;
beh.dv = dv;
beh.ch = ch;
beh.acc = acc;