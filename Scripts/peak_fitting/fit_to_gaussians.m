function [ area_peak, darea_peak, prof_fit, deviation, derivatives, numerical_derivatives ] = fit_to_gaussians( d_align, xsel, const_width, xsel_error_to_width_ratio, PLOT_STUFF, DERIV_CHECK )
% FIT_TO_GAUSSIANS: Fits electrophoretic traces to sums of Gaussians -- fast due to no optimization of peak positions.
%
%  [ area_peak, darea_peak, prof_fit, deviation ] = fit_to_gaussians( d_align, xsel, xsel_error_to_width_ratio, const_width, PLOT_STUFF );
%
% Required inputs:
%  d_align     = input matrix of traces
%  xsel        = band locations
%
% Optional inputs:
%  const_width = width of gaussians [options; default = (1/4) * mean band 
%                  spacing, which appears appropriate for ABI capillaries]
%  xsel_error_to_width_ratio 
%              = Ratio of uncertainties in peak locations compared to 
%                 pek widths. [default 1.0]
%  PLOT_STUFF  = parameter used by GUI to suppress graphical 
%                   output. [0 = no plotting; 1 = start/fit/deviations; 
%                   2 = fit; 3 = show individual traces].  [default 1]
%  DERIV_CHECK = get numerical derivatives for peak areas 
%                   with respect to xsel. [default 0]
%
% Outputs:
%  area_peak   = fitted peak areas.
%  darea_peak  = error estimates on fitted peak areas, based on uncertainties
%                 of peak locations.
%
% (C) R. Das 2008-2010, 2013

if nargin < 2;  help( mfilename ); return; end;

area_peak = [];
prof_fit = 0 * d_align;
deviation = [];
if isempty(xsel); return; end;

% by default, assume that error in peak location is about the same as typical peak width.
if ~exist( 'xsel_error_to_width_ratio', 'var') || isempty( xsel_error_to_width_ratio ); xsel_error_to_width_ratio = 0.5; end;
if ~exist( 'PLOT_STUFF', 'var') ; PLOT_STUFF = 1; end;
if ~exist( 'DERIV_CHECK', 'var'); DERIV_CHECK = 0; end;

global verbose;
verbose = [0 1];
if size( xsel, 1 ) == 1; xsel = xsel'; end;
if ~exist( 'const_width', 'var' ) || isempty( const_width ) ||  const_width == 0.0
  const_width = 0.25 * median( abs( diff(xsel) ) );
end
fprintf( 'Assuming constant peak width: %5.1f\n', const_width );
%area = [];

%x = 1:size( d_align, 1);
%numpeaks = size( xsel, 1 );
num_xsel_lanes = size( xsel, 2 ) ;

% Fit xsel, using a constant width approximation
xsel_start = xsel(:, 1)';
widthpeak_start = const_width + 0 * xsel_start;

if ~exist('whichpos', 'var');  whichpos = 1:size(d_align,2); end;

area_peak = zeros( length( xsel_start ), size( d_align, 2) );
prof_fit = 0 * d_align;
deviation = zeros( 1, size( d_align, 2) ) ;

% this is a vestige of the old routine where we allowed xsel to be different in differnt lanes. No longer supported!
if num_xsel_lanes == 1
  for j = 1:size(d_align,2)
    xsel(:,j) = xsel(:,1);
  end
end

% new! derivatives
for n = 1:size( d_align,2); derivatives{n} = []; end;
for n = 1:size( d_align,2); numerical_derivatives{n} = []; end;

%%%%%%%%%%%%%%%%%%
% inner loop
%%%%%%%%%%%%%%%%%%
if parallelization_exists() && (PLOT_STUFF ~= 3)
  parfor j = whichpos
    fprintf( 1, 'Fitting profile %d\n',j);
    [ area_peak(:,j), prof_fit(:,j), deviation(j), derivatives{j}, numerical_derivatives{j} ] = do_the_fit_inner_loop( d_align(:,j), xsel(:,j)', widthpeak_start, PLOT_STUFF, DERIV_CHECK );
  end;
else
  for j = whichpos
    fprintf( 1, 'Fitting profile %d\n',j);
    [ area_peak(:,j), prof_fit(:,j), deviation(j), derivatives{j}, numerical_derivatives{j} ] = do_the_fit_inner_loop( d_align(:,j), xsel(:,j)', widthpeak_start, PLOT_STUFF, DERIV_CHECK );
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% error estimation
% amount of uncertainty expected.
xsel_error = const_width * xsel_error_to_width_ratio;
darea_peak = [];
for j = whichpos;   
  darea_peak(:,j) = sqrt( sum( derivatives{j}.^2  * xsel_error^2)  );
end

%Npeaks = length( xsel_start );
%xsel_error = const_width * xsel_error_to_width_ratio * ones( 1, Npeaks );
%% if uncertainty exceeds distance to a neighboring peak, its probably an overesimate (this happens in band-compressed G's).
%xsel_error_max = xsel_error*0;
%if length( xsel_start > 1 ) 
%  xsel_error_max(1)   = abs(xsel_start(2) - xsel_start(1)); 
%  xsel_error_max(end) = abs(xsel_start(end) - xsel_start(end-1)); 
%end;
%for k = 2:(Npeaks-1);  xsel_error_max(k) = min( abs( xsel_start(k+1)-xsel_start(k)), abs( xsel_start(k)-xsel_start(k-1)) ); end;
%plot( [1:length(xsel_start)], [xsel_error; xsel_error_max/2] ); pause;
%xsel_error = max( [xsel_error; xsel_error_max/2] );
%for j = whichpos;   
%  darea_peak(:,j) = sqrt( derivatives{j}'.^2  * (xsel_error'.^2) );
%end;


if PLOT_STUFF == 1
  h = figure(3);
  set(h, 'Position', [0, 0, 800, 600]);
  set(h, 'PaperOrientation', 'landscape', 'PaperPositionMode', 'auto', 'color', 'white');
  
  subplot(1,3,1);
  scalefactor = 40/mean(mean(d_align));
  image( scalefactor * d_align );
  ylim( [ min(xsel_start)-100, max(xsel_start)+100 ] );
  title( 'data' );
  
  subplot(1,3,2);
  image( scalefactor * prof_fit );
  ylim( [ min(xsel_start)-100, max(xsel_start)+100 ] );
  title( 'fit' );
  
  subplot(1,3,3);
  image( scalefactor * ( d_align - prof_fit) );
  %easycolorplot( 100*(d_align -prof_fit) );
  ylim( [ min(xsel_start)-100, max(xsel_start)+100 ] );
  title( 'residuals' );
  
  
  colormap( 1 - gray(100) );
elseif PLOT_STUFF == 2
  figure(3)
  scalefactor = 40/mean(mean(d_align));
  image( scalefactor * prof_fit );
  ylim( [ min(xsel_start)-100, max(xsel_start)+100 ] );
  title( 'fit' );
  colormap( 1 - gray(100) );
end

% beep notice when finished
beep on; beep; beep off;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ area_peak, prof_fit, deviation, derivatives, numerical_derivatives ] = do_the_fit_inner_loop( d_align, xsel_fit, widthpeak, PLOT_STUFF, DERIV_CHECK )

tic

%x_allow = [1:size(d_align,1)];
[ x_disallow, x_allow ] = find_x_disallow( d_align );

[xdummy, xpeak_grid]    = meshgrid(x_allow, xsel_fit);
[xdummy,widthpeak_grid] = meshgrid(x_allow, widthpeak);
gaussian = exp( -0.5 * ((xdummy-xpeak_grid)./widthpeak_grid).^2 );

A = gaussian * gaussian';
B = gaussian * d_align(x_allow);

Ainv = inv( A );
%a  = lsqr( A, B )';
a = (Ainv * B)';
area_peak = sqrt(2*pi) * a .* widthpeak;

toc

% calculate derivatives analytically
N = length(xsel_fit);
dgaussian = ((xdummy-xpeak_grid)./(widthpeak_grid.^2)) .* exp( -0.5 * ((xdummy-xpeak_grid)./widthpeak_grid).^2 );
% term involving dot product of profile  and gaussian width
derivatives =  Ainv .*  repmat( (dgaussian * d_align( x_allow ) )', N, 1 );
% terms involving cross-gaussian overlaps -- tensor madness!
D = (gaussian * dgaussian')';
derivatives = derivatives - Ainv .* repmat( a*D', N, 1) - (Ainv * D') .* repmat(a, N, 1 );
% scale up to area.
derivatives = derivatives * sqrt(2*pi) .* repmat( widthpeak', 1, N);

% numerical derivatives if desired.
numerical_derivatives = derivatives * 0;
if DERIV_CHECK;  numerical_derivatives = calculate_numerical_derivatives( xsel_fit, widthpeak, widthpeak_grid, x_allow, d_align, area_peak ); end

area = area_peak;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = 1:size( d_align, 1);
profile_fit = 0 * x;

for k = 1:N
  predgaussian(:, k) = getgaussian(x, xsel_fit(k), widthpeak(k), a(k)); 
  profile_fit = profile_fit + predgaussian(:, k)';
end;

prof_fit = profile_fit;

deviation = sum(( profile_fit( x_allow )' - d_align( x_allow ) ).^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if PLOT_STUFF;
  clf;
  plot( d_align, 'k' );
  hold on; 
  for k = 1:N; plot( predgaussian(:,k),'b'); end;

  for i = 1:length( xsel_fit )
    if ( round( xsel_fit(i)) < length( profile_fit ) )
      plot( [xsel_fit(i) xsel_fit(i)], [0 profile_fit(round(xsel_fit(i)))],'k' ); hold on;
    end;
  end;

  plot(prof_fit, 'r');
  hold off;
  set(gca, 'ylim', [-0.5 5]);
  axis([ min(xsel_fit)-100 max(xsel_fit)+100 min(prof_fit) max( prof_fit )]);
  %pause;
end;

return;

  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  numerical_derivatives = calculate_numerical_derivatives( xsel_fit, widthpeak, widthpeak_grid, x_allow, d_align, area_peak );

dx = 0.01;

for i = 1:length( xsel_fit )
  fprintf( 'Calculating numerical derivatives for %d\n',i);
  xsel_test = xsel_fit;
  xsel_test(i) = xsel_test(i) + dx;
  [xdummy,xpeak_grid]    = meshgrid(x_allow, xsel_test);
  gaussian = exp( -0.5 * ((xdummy-xpeak_grid)./widthpeak_grid).^2 );
  A = gaussian * gaussian';
  B = gaussian * d_align(x_allow);
  Ainv = inv( A );
  a = (Ainv * B)';
  area_peak_test = sqrt(2*pi) * a .* widthpeak;
  numerical_derivatives(:, i) = (area_peak_test - area_peak) / dx;
end
