function [voi] = voigt(x, y0, I0, wL, wG, xc)
% 10/31/2017 - E R Louden 
% updated comments on 11/29/2017

% VOIGT: defines a fit type for a voigt profile (a convolution of gaussian & lorentzian) 
%   voi   -   dependent data (voigt profile, here intensity)
%   x   -   independent data (here phi)
%   y0  -   background
%   I0  -   integrated intensity of the peak
%   wL  -   width of Lorentzian (intrinsic width)
%   wG  -   width of Guassian (experimental resolution)
%   xc  -   center of profile
 
% You will either need the GAUSS and LORENTZ functions, or adjust the code at lines 52 & 59 below
% GAUSS and LORENTZ were written to match the forms specified by GRASP

% To perform fitting with this function: 
%   1. Define the fit type using this function, 
%         problem variables must be passed and are essentially "fixed"
% 
%         ft_voigt = fittype('voigt(x, y0, I0, wL, wG, xc)','problem','wG');
%         
%   2. Specify any fit options
%         e.g. what method should be used or starting/upper bounds/lower bounds for fitting variables
%         
%         fo_voigt = fitoptions('Method','NonlinearLeastSquares',...
%                             'StartPoint', [I_s wL_s, xc_s y0_s],...    
%                             'Upper', [I_u wL_u, xc_u y0_u],...
%                             'Lower',[I_l wL_l, xc_l y0_l]);
%                         
%   3. Call the Matlab fit function, passing your fit type and fit options
%         It will return a fit object and goodness of fit parameters
%         
%         [f_voigt, gof] = fit(phi, I, ft_voigt, fo_voigt,'problem',dwL(a));
%         
%   4. Get the confidence intervals
%         This will return an upper and lower bound on the fitted value
%         
%         cfi = confint(f_voigt);


%%
% Array of x-values where Voigt is calculated - these are currently passed as an input variable
% This allow us to compute results of irregular spaced x's as this most likely the case for the RC.
%x = -10:0.05:10;

% Array of x' over which integration is done.
dxp = min([wG wL])/10;
xp = -4*wG:dxp:4*wG;

% Also computation of Gaussian vector which only has to be done once.
%gau = 1/(wG*sqrt(2*pi))*exp(-xp.^2/(2*wG^2));
gau = GAUSS(xp, 0, 1, wG, 0);

% Computation of Voigt for each value of x. 
voi = x;
for i=1:length(x)
    % Computs the Lorentzian for each value of x
    %lor = (wL/pi)./((x(i)-xc-xp).^2 + wL^2);
    lor = LORENTZ(x(i), 0, 1, wL, xc+xp);
    
    %dots it with Gaussian to do integration.
    voi(i) = dxp*dot(gau,lor);
end

% add the constant background and multiply by integrated intensity
voi = I0*voi + y0;

%% Extras
% Normalization checks 
lor = (wL/pi)./(x.^2 + wL^2);
disp(num2str(dxp*sum(gau)));
disp(num2str(0.05*sum(lor)));
disp(num2str(0.05*sum(voi)));

% Comparison of Gaussian, Lorentzian and Voigt
%lor = (wL/pi)./((x-xc).^2 + wL^2);
%gau = 1/(wG*sqrt(2*pi))*exp(-(x-xc).^2/(2*wG^2));
%hold on
%plot(x,gau)
%plot(x,lor)
%plot(x,voi)

end
