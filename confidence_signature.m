function [cf, acc, stm] = confidence_signature(noise, ntr, conftype, stmdist, overlap)
%% simulation of the signatures of decision confidence 
% INPUT: 
% noise ... internal noise: default is 22.8
% ntr ... the number of simulated trials
% conftype ... 'sdt' (Hangya et al., 2016) or 'Bayes' (Adler & Ma, 2017) 
% stmdist ... type of stimulus distribution:  'uniform' or 'Gaussian'
% overlap ... overlap (0 or >1) between P(s|C=-1) and P(s|C=1). Default is 0. 
%
% OUTPUT: cf (confidence: 0.5 - 1)
%         acc (accuracy: 0, error; 1, correct)
%         stm (signed stimulus)
%
% Monte Carlo simulation after Hangya et al., 2016, Adler & Ma, 2017 etc
%
% EXAMPLE: simply run '[cf, acc, stm] = confidence_signature;
%
% To plot the results, use 'plot_confidence_signature.m'
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

close all;

% deal with inputs
if nargin < 1; noise = 22.8; end
if nargin < 2; ntr = 10^7; end
if nargin < 3; conftype = 'sdt'; end
if nargin < 4; stmdist = 'uniform'; end
if nargin < 5; overlap = 0; end

% evidence strength
dc = -overlap:1:50;

% stimulus distribution
switch lower(stmdist)
    case 'uniform' 
        lendc = length(dc);
        pc1 = ones(1, lendc);
        pc2 = ones(1, lendc);
    case 'gaussian' 
        stmMean = 25;
        stmSD = 15;
        pc1 = gaussmf(dc,[stmSD stmMean]);
        pc2 = gaussmf(-dc,[stmSD -stmMean]);
end

% true categories
disp(['simulating ' num2str(ntr) ' trials...'])
C = datasample([-1 1], ntr);

% signed stimulus
stm = C;
stm(C==1) = datasample(dc, sum(C==1), 'Weights', pc1);
stm(C==-1) = datasample(-dc, sum(C==-1), 'Weights', pc2);

% noisy measurements (decision variable)
dv = arrayfun(@(x) normrnd(x, noise), stm);

% confidence & choice
switch lower(conftype)
    case 'sdt'
        % confidence as proportion correct (Hangya et al., 2016) 
        cf = 0.5 + 0.5*erf(abs(dv)/(noise*sqrt(2)));
        ch = sign(dv);
    case 'bayes'
        % Bayesian confidence (Adler & Ma, 2017)
        switch lower(stmdist)
            case 'uniform'
                likelihood1 = normcdf(dc(end), dv, noise) - normcdf(0, dv, noise);
                likelihood2 = normcdf(0, dv, noise) - normcdf(-dc(end), dv, noise);
            case 'gaussian'
                sumSD = sqrt(stmSD^2 + noise^2);
                likelihood1 = normpdf(dv, stmMean, sumSD);
                likelihood2 = normpdf(dv, -stmMean, sumSD);
        end
        cf = likelihood1./(likelihood1 + likelihood2);
        ch = ones(1, ntr);
        ch(cf < 0.5) = -1;
        cf(cf < 0.5) = 1 - cf(cf < 0.5);
end

% accuracy
acc = 1*(ch==C);