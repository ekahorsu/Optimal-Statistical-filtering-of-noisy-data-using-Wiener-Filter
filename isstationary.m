% +------------------------------------------------------+
% |      Signal wide-sense stationarity estimation       |
% |              with MATLAB Implementation              | 
% |                                                      |
% | Author: Ph.D. Eng. Hristo Zhivomirov        08/22/22 | 
% +------------------------------------------------------+
% 
% function: [wss_stat_flag, mean_stat_flag, var_stat_flag, cov_stat_flag] = isstationary(x, gamma)
%
% Input:
% x - signal in the time domain;
% gamma - confidence level for the hipothesis that the signal under
%         consideration is stationary (between 0 and 1 e.g., 0.9; 0.95; 0.99).
% 
% Output:
% wss_flag - a Boolean flag showing whether the signal is wide-sense   
%            stationary (WSS), that is, simultaneously stationary about   
%            its mean, variance and autocovariance;
% mean_stat_flag - a Boolean flag showing whether the signal is stationary
%                  about its mean;
% var_stat_flag - a Boolean flag showing whether the signal is stationary
%                 about its variance; 
% cov_stat_flag - a Boolean flag showing whether the signal is stationary
%                 about its autocovariance. 
%
% Note: a signal (e.g., time series) is said to be weakly stationary or
% wide sense stationary (WSS) if its (time-localized) mean and variance 
% are constant over time and if its autocovariance function Cxx(t1, t2) 
% depends only on the difference t2-t1, but not on their particular values. 
% Be aware - one does not consider the underling process or the entire
% population here! The function estimates the WSS of the signal "as it is"!

function [wss_flag, mean_stat_flag, var_stat_flag, cov_stat_flag] = isstationary(x, gamma)

% input validation
validateattributes(x, {'single', 'double'}, ...
                      {'vector', 'real', 'nonnan', 'nonempty', 'finite'}, ...
                      '', 'x', 1)
validateattributes(gamma, {'single', 'double'}, ...
                          {'scalar', 'real', 'nonnan', 'nonempty', 'finite', '>=', 0, '<=', 1}, ...
                          '', 'gamma', 2)
                      
% determine the signal length
xlen = length(x);
                      
% split the signal into two equally length parts
if rem(xlen, 2), x(end) = []; end
x1 = x(1:end/2); x2 = x(end/2+1:end);

% set the significance level
% Note: the significance level (alpha) complements the confidence level 
% (gamma) to 1.
alpha = 1-gamma;

% check for the statistical significance of the first moment (i.e., the
% mean) using the inverse form factor (i.e., the MEAN to STD ratio) and
% then empirically test for first moment stationarity
m1 = mean(x1); m2 = mean(x2);
s1 = std(x1); s2 = std(x2);
if (abs(m1)/s1 < alpha) && (abs(m2)/s2 < alpha)                             
    % when both mean values are not statistically significant...
    mean_stat_flag = true;   
else
    % when at least one of the mean values is statistically significant...
    mean_stat_flag = abs((m1 - m2)/min(m1, m2)) < alpha;                    
end
          
% empirically test for second moment (i.e., variance) stationarity   
v1 = var(x1); v2 = var(x2);
var_stat_flag = abs(v1 - v2)/min(v1, v2) < alpha;                           

% empirically test for autocovariance stationarity via: (i) computation of
% the autocovariances of the two signal halves; (ii) estimation of the
% similarity between them using the sample Pearson correlation coefficient
% and (iii) additional check by comparison of the variances of the
% autocovariance sequences of the two signal halves.
% Note: the p-value obtained by the corrcof function ranges from 0 to 1,
% where values close to 0 (less than alpha) correspond to a significant
% correlation and a low probability of observing the null hypothesis that
% there is no relationship between the sequences. For more information
% about the computation of the p-value see the script of the corrcoef
% function.
c1 = xcov(x1, 'coeff'); c2 = xcov(x2, 'coeff'); 
vc1 = var(c1); vc2 = var(c2);
[~, p] = corrcoef([c1(:) c2(:)]);
cov_stat_flag = abs(vc1 - vc2)/min(vc1, vc2) < alpha && p(1, 2) < alpha;    
                  
% check the overall wide-sense stationarity
wss_flag = mean_stat_flag && var_stat_flag && cov_stat_flag;         

end