function [dev_profile, jacobian] = ...
    evaluate_deviation_varywidonly_gaussian(pixels, profiles, params)
% [dev_profile, jacobian] = ...
%    evaluate_deviation_varywidonly_gaussian(pixels, profiles, params)


%dev_profile = predict_profile_varywidonly_gaussian( pixels, params' )  -   profiles';
%jacobian = predict_partials_varywidonly_gaussian( pixels, [], params', ...
%						  [], []);
params = params';

% Try to save time by combining calculations.

global xpeak;

numpeaks = (length(params))/2;
%Function assumes you are submitting log of the amplitudes!
amppeak  = exp(params(  numpeaks+1: 2*numpeaks));
widthpeak = params(1:numpeaks);

[x,xpeak_grid]    = meshgrid(pixels( (numpeaks+1):end ) ,xpeak);
[x,widthpeak_grid]= meshgrid(pixels( (numpeaks+1):end ),widthpeak);
[x,amppeak_grid]   = meshgrid(pixels((numpeaks+1):end),amppeak);

gaussian = exp( -0.5 * ((x-xpeak_grid)./widthpeak_grid).^2 );

fitprofile=gaussian'*amppeak;

%Include values of xpeak at end of profile ... used to regularize
%peak positions.
global width_constraint_strength
fitprofile = [ width_constraint_strength*widthpeak' fitprofile']';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Partial derivatives with respect to amplitudes -- actually log(zeta)!
prt(numpeaks+1:2*numpeaks,  :) = [zeros(numpeaks,numpeaks) ... 
		    amppeak_grid.*gaussian ];

%Partial derivatives with respect to widths. 
widthprt = +amppeak_grid.*(gaussian).*((xpeak_grid-x).^2) ./widthpeak_grid.^3;
%dwid_dminwid    = ones(1,numpeaks); 
prt(1:numpeaks,:) =  [zeros(numpeaks,numpeaks) widthprt ];

global width_constraint_strength
for k = 1:numpeaks
  prt(k,k) = width_constraint_strength;
end

f = fitprofile;
dev_profile = f - profiles';
jacobian = prt';
