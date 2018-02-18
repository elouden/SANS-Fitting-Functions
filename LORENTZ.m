function [ y ] = LORENTZ(x, y0, I0, wL, xc)
% 10/31/2017 - E R Louden 

% LORENTZ: calculates a lorentzian
%   y   -   dependent data (lorentzian profile, here intensity)
%   x   -   independent data (usually phi, san, or theta)
%   y0  -   background
%   I0  -   integrated intensity of the peak
%   wL  -   width of Lorentzian (e.g. the intrinsic width)
%   xc  -   center of profile
 
% This matches the form used by GRASP

y = y0+((2.*I0./pi)*wL)./(4.*((x-xc).^2)+(wL.^2));
  

end

