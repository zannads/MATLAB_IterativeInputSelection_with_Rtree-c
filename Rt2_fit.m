function R = Rt2_fit(yo,ys)

% This function calculates the value of the coefficient of determinarion R2
%
% Input:
% yo = observed data
% ys = simulated (or predicted) data
%
% Output
% R = coefficient of determination

narginchk(2,2);
nargoutchk(1,1);

R = 1-(cov(yo-ys)/cov(yo));