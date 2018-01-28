function SimpleBaysianSampling
% implement a simple Bayesian sampling model proposed by Deneve, 2012

%%
% preset free parameters
hdx = 0.3*[-1 -0.75 -0.5 -0.25 0 0.25 0.5 0.75 1];
kernel = [0.1 0.5 0.2 0.05 0 -0.05 -0.2 -0.5 -0.1];
len_tr = 5000;
len_frame = 100;
stmstrength = linspace(-0.5, 0.5, 9);

q = [1 5 2 0.5 0 -0.5 -2 -5 -1]; % baseline firing rate of sensory neurons
dq = 2; % increment of sensory evidence represented by sensory neurons
L0 = 0;  % prior
DThre = 4;  % decision threshold
C0 = 0;  % initial coherence estimate
neu = 1.5;  % estimate update rate
% dynamic threshold, sensory weighting, starting point of sensory
% integration

