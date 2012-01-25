function [xpeak_fit,widthpeak,amppeak,areapeak,f] ...
    = fitset_inputwid_gaussian(Profile_tofit,numres,xsel,width_input,x_disallow)

figure(2);
subplot(1,1,1);
hold off;
numpeaks = length(numres);

if ~exist( 'x_disallow')
  x_disallow = [];
end

xpeak = real(xsel(numres));
amppeak = Profile_tofit( round(xpeak) )';

% Provide default Lorentzian coefficients if not passed in as arguments
if (nargin<4) const_width = 7; end;
global widthpeak;
widthpeak = width_input';

minbin = max( round(min(xpeak)) - 100, 1);
maxbin = min( round(max(xpeak)+max(widthpeak))+150, length(Profile_tofit)) ;
pixels   = [minbin:maxbin]';
pixels = setdiff( pixels, x_disallow);
profiles = Profile_tofit( pixels )';

%Add in xsel constraints...
global xsel_constraint_strength;
xsel_constraint_strength = 0.1;
profiles = [   xsel_constraint_strength*xpeak profiles ];
pixels =   [          minbin*ones(1,numpeaks) pixels'  ]';
 
% pin=[xpeak log(amppeak)];
% 
% stol=0.001;
% niter=10;
% wt=0*pixels+1;
% dp=0*pin+1;


n = 1:length(xpeak);
param(3*(n-1) + 2) = xpeak;
param(3*(n-1) + 1) = log(amppeak);
param(3*(n-1) + 3) = 3;

% options = [0*pin; widthpeak'/100, Inf*amppeak]';
% %options = [0*pin; widthpeak'/100, 5*amppeak, 0.1,0.5]';
% %options = [0*pin; widthpeak'/10, widthpeak'/10, Inf*amppeak]';
% %options = [widthpeak/100, widthpeak/100, amppeak/100; Inf*pin]';
% params = pin;

[iter p c] =levmar('gaussianfit', 'jacgaussian', param, profiles, 50, [], 1:length(profiles));

%Try to use MATLAB's optimization toolbox for faster speed.
% OPTIM_TOOLBOX = exist( 'lsqnonlin' );
% if OPTIM_TOOLBOX
%   options = optimset( 'LevenbergMarquardt','on','Jacobian','on','TolFun',stol,'MaxIter',50,'Display','off');
%   fun = @( pin )evaluate_deviation_inputwid_gaussian( pixels, profiles, pin  );
%   params = lsqnonlin( fun, pin, [], [], options )';
%   [dev_profile, jacobian, f ] = fun( params' );
% else
%   [f,params,kvg,iter,corp,covp,covr,stdresid,Z,r2]= ...
%       leasqr(pixels, profiles,pin,'predict_profile_inputwid_gaussian',stol,niter,wt,dp,'predict_partials_inputwid_gaussian',options);
% end

% Extract the peak positions, their amplitudes and the 
% two Lorentzian coefficients from params
xpeak_fit= p( 3*(n-1) + 2)';
% Now take the exponential of the fitted peak amplitudes
amppeak  = p(params(3*(n-1) + 1))';
%width_fit  = params(2*numpeaks+1);

% Get the distance between peaks
%distpeak = getdistpeak(xpeak);
% Use our Lorentzian coefficients and the distances 
% betweeen peaks to get our Lorentzian widths
%widthpeak = width_fit + 0*xpeak_fit;

% Plot all the Lorentzians
%hold on
%for k=1:numpeaks
%    predlorentz = getlorentzian(pixels,xpeak_fit(k), widthpeak(k),amppeak(k)); 
%    %plot(pixels, predlorentz,'r');
%end

%hold on
%for k=1:numpeaks; plot([xpeak(k) xpeak(k)],[0 amppeak(k)],'k'); hold on; end;
hold off

areapeak = sqrt(2*pi) *amppeak.*widthpeak'; 
%subplot(3,1,2);
%plot(numres,areapeak,'color','k');
%axis([min(numres) max(numres) 0 max(areapeak)]);
shg
