function [ y ] = GAUSS(x, y0, I0, wG, xc)
% 10/31/2017 - E R Louden 

% GAUSS: calculates a guasian
%   y   -   dependent data (gaussian profile, here intensity)
%   x   -   independent data (usually phi, san, or theta)
%   y0  -   background
%   I0  -   integrated intensity of the peak
%   wG  -   width of Guassian (e.g. the width due to experimental resolution)
%   xc  -   center of profile

% This matches the form used by GRASP

 y = y0 + (I0 ./ (wG .* sqrt(pi ./ 2) ./ sqrt(log(4)))).*exp((-2.*(x-xc).^2)/((wG.^2)/log(4)));
 
end
