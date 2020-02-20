%% ErrorMag.m
% Given an xdata and ydata this function tries to fit it with a thresholded
% gaussian.
%
% A gaussian put through this function
%           x - phi     if x > phi
% f(x) =        0       if |x| < phi
%           |x + phi|   if x < -phi
% The sparseness is the proportion of the weight in the distribution that
% is not set to zero.
%
% Alpha is just a conversion factor. It has units, therefore converts from
% the units of the data to numbers (the implicit units we are using for
% phi)
function [E] = ErrorMag(alpha, xdata, ydata,Sy)

% Here we can see that the ratio of phi and the varinance is set by Sy, the
% sparseness
phioversigma = sqrt(2)*erfcinv(Sy); 

% Then we fit the data using this curve
Fit = @(x) sqrt(1/(2*pi))*alpha/Sy.*exp(-0.5.*(alpha.*x+sign(x).*phioversigma).^2);

% Before finding its squared error of fit
E = sum((ydata - Fit(xdata)).^2);
end