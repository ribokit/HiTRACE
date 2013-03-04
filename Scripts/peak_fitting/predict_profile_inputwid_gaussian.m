function fitprofile = predict_profile_inputwid_gaussian(pixels,params)
%  fitprofile = predict_profile_inputwid_gaussian(pixels,params)

numpeaks = (length(params))/2;
xpeak    = params(  1         : numpeaks);
%Function assumes you are submitting log of the amplitudes!
amppeak  = exp(params(  numpeaks+1: 2*numpeaks));
global widthpeak;

[x,xpeak_grid]    = meshgrid(pixels((numpeaks+1):end), xpeak);
[x,widthpeak_grid]= meshgrid(pixels((numpeaks+1):end), widthpeak);

%basisfunction = 1./(1+((x-xpeak_grid)./widthpeak_grid).^2);
basisfunction = exp( -0.5 * ((x-xpeak_grid)./widthpeak_grid).^2 );


fitprofile=basisfunction'*amppeak;

%Include values of xpeak at end of profile ... used to regularize
%peak positions.
global xsel_constraint_strength
fitprofile = [ xsel_constraint_strength*xpeak' fitprofile']';
