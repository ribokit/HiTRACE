function prt=predict_partials_inputwid_gaussian(pixels,f,params,dp,F)
%  prt=predict_partials_inputwid_gaussian(pixels,f,params,dp,F)

numpeaks = (length(params))/2;
xpeak    = params(  1         : numpeaks);
% Function assumes you are submitting log of the amplitudes!
amppeak  = exp(params(numpeaks+1: 2*numpeaks));
global widthpeak;

[x,xpeak_grid]     = meshgrid(pixels((numpeaks+1):end),xpeak);
[x,widthpeak_grid] = meshgrid(pixels((numpeaks+1):end),widthpeak);
[x,amppeak_grid]   = meshgrid(pixels((numpeaks+1):end),amppeak);

%Partial derivatives of a Lorentzian.
%The basis functions:
gaussian = exp( -0.5 * ((x-xpeak_grid)./widthpeak_grid).^2 );
%Partial derivatives with respect to x (peak positions).
prt(1:numpeaks,             :) = [ zeros(numpeaks,numpeaks) ...
		    -amppeak_grid.*gaussian.*(xpeak_grid- ...
						  x)./widthpeak_grid.^2 ];
%Partial derivatives with respect to amplitudes -- actually log(zeta)!
prt(numpeaks+1:2*numpeaks,  :) = [ zeros(numpeaks,numpeaks) ...
		    amppeak_grid.*gaussian ];

global xsel_constraint_strength
for k = 1:numpeaks
  prt(k,k) = xsel_constraint_strength;
end


prt = prt';
