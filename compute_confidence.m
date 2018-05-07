function [cf, ch] = compute_confidence(dv, noise, conftype, stmdist, dc, stmMean, stmSD)
% estimate confidence based on decision variable and internal noise
%
% INPUT:
% dv ... decision variables: trials x frames
% noise ... internal noise
% conftype ... 'sdt' (Hangya et al., 2016) or 'Bayes' (Adler & Ma, 2017) 
% stmdist ... type of stimulus distribution:  'uniform' or 'Gaussian'
% dc ... presented stimulus for category 1
% stmMean ... mean of the Gaussian stimulus distribution for category 1
% stmSD ... SD of the Gaussian stimulus distribution for category 1
%
% OUTPUT:
% cf ... computed trial-by-trial decision confidence
% ch ... trial-by-trial choice (-1 or 1)
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++

% deal with inputs
if nargin < 1; error('Provide decision variables!'); end
if nargin < 2; noise = 22.8; end
if nargin < 3; conftype = 'sdt'; end
if nargin < 4; stmdist = 'uniform'; end
if nargin < 5; dc = 0:1:50; end
if nargin < 6; stmMean = 25; end
if nargin < 7; stmSD = 15; end

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
        ch = ones(size(dv,1), 1);
        ch(cf < 0.5) = -1;
        cf(cf < 0.5) = 1 - cf(cf < 0.5);
end