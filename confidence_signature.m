function [cf, acc, stm] = confidence_signature(stmdist, noise, ntr)
%% simulation of the signatures of decision confidence 
% INPUT: 
% stmdist ... type of stimulus distribution: 
% 'uniform', 'Gaussian'
% noise ... internal noise: default is 22.8
% ntr ... the number of simulated trials
%
% Monte Carlo simulation based on Hangya et al., Neural Computation (2016)
% To plot the results, use 'plot_confidence_signature.m'
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

close all;

% deal with inputs
if nargin < 1; stmdist = 'uniform'; end
if nargin < 2; noise = 22.8; end
if nargin < 3; ntr = 10^7; end

% evidence strength
dc = 0:1:50;

% stimulus distribution
switch stmdist
    case 'uniform'
        lendc = length(dc);
        weights = ones(1, lendc)/lendc;
    case 'Gaussian'
        weights = gaussmf(dc,[15 25]);
        weights = weights/sum(weights);
end

% signed stimulus
ss = sort(unique([-dc dc]));
sw = [weights weights];
sw(length(weights)) = [];

% assign stimulus
disp(['simulating ' num2str(ntr) ' trials...'])
stm = datasample(ss, ntr, 'Weights', sw);

% noisy measurements (decision variable)
dv = arrayfun(@(x) normrnd(x, noise), stm);

% choice
ch = sign(dv);

% (Bayesian) confidence 
cf = 0.5 + 0.5*erf(abs(dv)/(noise*sqrt(2)));

% accuracy
acc = 1*(ch==sign(stm));
acc(stm==0) = randi(2, 1, sum(stm==0)) - 1;