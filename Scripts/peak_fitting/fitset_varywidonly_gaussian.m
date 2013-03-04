function [ widthpeak,amppeak,areapeak ] ...
    = fitset_varywidonly_gaussian(Profile_tofit,numres,xsel, widthpeak, ...
				  x_disallow)
%
% [ widthpeak,amppeak,areapeak ] ...
%    = fitset_varywidonly_gaussian(Profile_tofit,numres,xsel, widthpeak,  x_disallow)
%
%

figure(2);
subplot(1,1,1);
hold off;
numpeaks = length(numres);

if ~exist( 'x_disallow')
  x_disallow = [];
end

global xpeak;

xpeak = real(xsel(numres));
amppeak = Profile_tofit(round(xpeak))';

% Provide default Lorentzian coefficients if not passed in as arguments
%if (nargin<4) const_width = 7; end;
%widthpeak = xpeak'*0 + const_width;
widthpeak = widthpeak';

minbin = max( round(min(xpeak)) - 100, 1);
maxbin = min( round(max(xpeak)+max(widthpeak))+150, length(Profile_tofit)) ;
pixels   = [minbin:maxbin]';
pixels = setdiff( pixels, x_disallow);
profiles = Profile_tofit( pixels )';

%Add in width constraints...
global width_constraint_strength;
width_constraint_strength = 20.0;
profiles = [   width_constraint_strength*widthpeak' profiles ];
pixels =   [               minbin*ones(1,numpeaks) pixels'  ]';


pin=[widthpeak' log(amppeak)];

stol=0.001;
niter=50;
wt=0*pixels+1;
dp=0*pin+1;

options = [0*pin; 10*widthpeak', Inf*amppeak ]';
%options = [0*pin; widthpeak'/100, 5*amppeak, 0.1,0.5]';
%options = [0*pin; widthpeak'/10, widthpeak'/10, Inf*amppeak]';
%options = [widthpeak/100, widthpeak/100, amppeak/100; Inf*pin]';
%params = pin;


%Try to use MATLAB's optimization toolbox for faster speed.
OPTIM_TOOLBOX = exist( 'lsqnonlin' );
if OPTIM_TOOLBOX
  options = optimset( 'LevenbergMarquardt','on','Jacobian','on','TolFun',stol,'MaxIter',50);
  %options = optimset( 'Algorithm','','Jacobian','on','TolFun',stol,'MaxIter',50,'Display','off');
  fun = @( pin )evaluate_deviation_varywidonly_gaussian( pixels, ...
						  profiles, pin  );
  params = lsqnonlin( fun, pin, [], [], options )';
else
  [f,params,kvg,iter,corp,covp,covr,stdresid,Z,r2]= ...
      leasqr(pixels, profiles,pin,'predict_profile_varywidonly_gaussian',...
	     stol,niter,wt,dp,'predict_partials_varywidonly_gaussian',options);
end


% Extract the peak positions, their amplitudes and the 
% two Lorentzian coefficients from params
widthpeak =params(  1         : numpeaks)';
% Now take the exponential of the fitted peak amplitudes
amppeak  = exp(params( numpeaks+1: 2*numpeaks))';

% constant during run...
xpeak_fit = xpeak;

% Get the distance between peaks
%distpeak = getdistpeak(xpeak);
% Use our Lorentzian coefficients and the distances 
% betweeen peaks to get our Lorentzian widths
%widthpeak = width_fit + 0*xpeak_fit;

% Plot all the Lorentzians
hold on
for k=1:numpeaks
    predgauss = getgaussian(pixels,xpeak_fit(k), widthpeak(k),amppeak(k)); 
    %plot(pixels, predlorentz,'r');
end

%hold on
%for k=1:numpeaks; plot([xpeak(k) xpeak(k)],[0 amppeak(k)],'k'); hold on; end;
hold off

areapeak = sqrt(2*pi)*amppeak.*widthpeak; 
%subplot(3,1,2);
%plot(numres,areapeak,'color','k');
%axis([min(numres) max(numres) 0 max(areapeak)]);
%shg

