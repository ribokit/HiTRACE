function [ areas, prof_fit, deviation ] = do_the_fit_fast( d_align, xsel, const_width, PLOT_STUFF );
% DO_THE_FIT_FAST: Fits electrophoretic traces to sums of Gaussians -- fast due to no optimization of peak positions.
%
%  [ areas, prof_fit, deviation ] = do_the_fit_fast( d_align, xsel, const_width, PLOT_STUFF );
%
%  d_align = input matrix of traces
%  xsel = band locations
%  const_width = width of gaussians [options; default = (1/4) * mean band spacing, which appears appropriate for ABI capillaries]
%  PLOT_STUFF = parameter used by GUI to suppress graphical output. [default 1]
%
% (C) R. Das 2008-2010

areas = [];
prof_fit = 0 * d_align;
deviation = [];
if length( xsel ) == 0; return; end;

if ~exist( 'PLOT_STUFF' ); PLOT_STUFF = 1; end;

global verbose;
verbose = [0 1];
if size( xsel, 1 ) == 1; xsel = xsel'; end;
if ~exist( 'const_width' ) | const_width == 0.0
  const_width = 0.25 * mean( abs(  xsel(2:end) - xsel(1:end-1) ) );
end
fprintf( 'Assuming constant peak width: %5.1f\n', const_width );
area = [];
VARY_XSEL = 1;
VARY_WIDTH = 1;
VARY_XSEL2 = 1;


x = 1:size( d_align,1);
numpeaks = size( xsel, 1 );

num_xsel_lanes = size( xsel, 2 ) ;

% Fit xsel, using a constant width approximation
xsel_start = xsel(:,1)';
widthpeak_start = const_width + 0 *xsel_start;

count = 0;

if ~exist('whichpos')
  whichpos = [1 : size(d_align,2)];
end

areas = zeros( length( xsel_start ), size( d_align, 2) );
prof_fit = 0 * d_align;
deviation = zeros( 1, size( d_align, 2) ) ;

if num_xsel_lanes == 1
  for j = 1:size(d_align,2)
    xsel(:,j) = xsel(:,1);
  end
end



%%%%%%%%%%%%%%%%%%
% inner loop
%%%%%%%%%%%%%%%%%%
if exist( 'matlabpool' )
  if matlabpool( 'size' ) == 0 ;   res = findResource; matlabpool( res.ClusterSize ); end
  parfor j= whichpos
    fprintf( 1, 'Fitting profile %d\n',j);
    [ areas(:,j), prof_fit(:,j), deviation(j) ] = do_the_fit_inner_loop( d_align(:,j), xsel(:,j)', widthpeak_start );
end
else
  for j= whichpos
    fprintf( 1, 'Fitting profile %d\n',j);
    [ areas(:,j), prof_fit(:,j), deviation(j) ] = do_the_fit_inner_loop( d_align(:,j), xsel(:,j)', widthpeak_start );
  end
end  


if PLOT_STUFF
  figure(3)
  
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
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ area, prof_fit, deviation ] = do_the_fit_inner_loop( d_align, xsel_start, widthpeak_start );

tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xsel_fit = xsel_start;
widthpeak = widthpeak_start;

[ x_disallow, x_allow ] = find_x_disallow( d_align );

[xdummy,xpeak_grid]    = meshgrid(x_allow, xsel_fit);
[xdummy,widthpeak_grid]= meshgrid(x_allow, widthpeak);
gaussian = exp( -0.5 * ((xdummy-xpeak_grid)./widthpeak_grid).^2 );

A = gaussian * gaussian';
B = gaussian * d_align(x_allow);
a  = lsqr( A, B )';
areapeak = sqrt(2*pi) * a .* widthpeak;

toc

area = areapeak;

plot( d_align,'k' );
hold on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = 1:size( d_align,1);
profile_fit = 0 * x;

numpeaks = length( xsel_fit );

for k=1:numpeaks
  predgaussian = getgaussian(x,xsel_fit(k), widthpeak(k),a(k)); 
  profile_fit = profile_fit + predgaussian;
  plot( predgaussian,'b');
end

prof_fit = profile_fit;

deviation = sum(( profile_fit( x_allow )' - d_align( x_allow ) ).^2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:length( xsel_fit )
  if ( round( xsel_fit(i)) < length( profile_fit ) )
    plot( [xsel_fit(i) xsel_fit(i)], [0 profile_fit(round(xsel_fit(i)))],'k' ); hold on;
  end
end

plot( prof_fit,'r');
hold off;
set(gca,'ylim',[-0.5 5]);
axis([ min(xsel_fit)-100 max(xsel_fit)+100 -0.5 max( prof_fit )]);

%pause; 
  
