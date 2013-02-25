function [dev_profile, jacobian, f] = ...
    evaluate_deviation_inputwid_gaussian(pixels, profiles, params)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Using existing routine.
%f = predict_profile_inputwid_gaussian( pixels, params' );
%dev_profile = f - profiles';
%jacobian = predict_partials_inputwid_gaussian( pixels, [], params', ...
%						  [], []);

params = params'; %stupid.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perhaps can save time by not repeating calculations...
numpeaks = (length(params))/2;
xpeak    = params(  1         : numpeaks);
%Function assumes you are submitting log of the amplitudes!
amppeak  = exp(params(  numpeaks+1: 2*numpeaks));
global widthpeak;

[x,xpeak_grid]    = meshgrid(pixels((numpeaks+1):end), xpeak);
[x,widthpeak_grid]= meshgrid(pixels((numpeaks+1):end), widthpeak);
[x,amppeak_grid]   = meshgrid(pixels((numpeaks+1):end),amppeak);

%gaussian = 1./(1+((x-xpeak_grid)./widthpeak_grid).^2);
gaussian = exp( -0.5 * ((x-xpeak_grid)./widthpeak_grid).^2 );

fitprofile=gaussian'*amppeak;

%Include values of xpeak at end of profile ... used to regularize
%peak positions.
global xsel_constraint_strength
fitprofile = [ xsel_constraint_strength*xpeak' fitprofile']';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Partial derivatives
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


f = fitprofile;
dev_profile = f - profiles';
jacobian = prt';
