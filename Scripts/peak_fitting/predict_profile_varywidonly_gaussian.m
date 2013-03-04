function fitprofile = predict_profile_varywidonly_gaussian(pixels,params)
% fitprofile = predict_profile_varywidonly_gaussian(pixels,params)

global xpeak;

numpeaks = (length(params))/2;
%Function assumes you are submitting log of the amplitudes!
amppeak  = exp(params(  numpeaks+1: 2*numpeaks));
widthpeak = params(1:numpeaks);

[x,xpeak_grid]    = meshgrid(pixels( (numpeaks+1):end ) ,xpeak);
[x,widthpeak_grid]= meshgrid(pixels( (numpeaks+1):end ),widthpeak);

basisfunction = exp( -0.5 * ((x-xpeak_grid)./widthpeak_grid).^2 );

fitprofile=basisfunction'*amppeak;

%Include values of xpeak at end of profile ... used to regularize
%peak positions.
global width_constraint_strength
fitprofile = [ width_constraint_strength*widthpeak' fitprofile']';
