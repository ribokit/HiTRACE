function [ areas, prof_fit, xsel_fit_all, width_peak_all ] = do_the_fit_GUI( d_align, xsel, whichpos );

if ( size(xsel,1) == 1  )
  xsel = xsel';
end


global verbose;
verbose = [0 1];
%const_width = 7;
const_width = 3;
area = [];
VARY_XSEL = 1;
VARY_WIDTH = 1;
VARY_XSEL2 = 1;
x = 1:size( d_align,1);
numpeaks = size( xsel, 1 );

num_xsel_lanes = size( xsel, 2 ) ;

% Fit xsel, using a constant width approximation
xsel_start = xsel(:,1)';
width_peak_start = const_width + 0 *xsel_start;
xsel_peak_all = [];
width_peak_all = [];

count = 0;

if ~exist('whichpos')
  whichpos = [1 : size(d_align,2)];
end

areas = zeros( length( xsel_start ), size( d_align, 2) );

parfor j= whichpos
  
  tic

  fprintf( 1, 'Fitting profile %d\n',j);


  x_disallow = find_x_disallow( d_align(:,j) );
  
  if ( num_xsel_lanes > 1 )
    xsel_start = xsel(:,j)';
  else
    xsel_start = xsel';
  end
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  xsel_fit = xsel_start;
  width_peak = width_peak_start;
  if VARY_XSEL
   [xsel_fit, w,a,areapeak ] = ...
    fitset_inputwid_gaussian( d_align(:,j), 1:numpeaks, xsel_start, ...
			      width_peak_start, x_disallow );
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Fit widths
  if VARY_WIDTH
    [width_peak, a,areapeak] = ...
	fitset_varywidonly_gaussian( d_align(:,j), 1:numpeaks, xsel_fit, ...
			 width_peak_start, x_disallow);
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if VARY_XSEL2
    [xsel_fit, w,a,areapeak,f ] = ...
	fitset_inputwid_gaussian( d_align(:,j), 1:numpeaks, xsel_fit, ...
				  width_peak, x_disallow );
  end

  % bizarro -- why do parameters return as complex numbers?
  % anyway, imaginary parts appear to be close zero.
  xsel_fit = real( xsel_fit );
  width_peak = real( width_peak );
  a = real( a );

  toc

  xsel_fit_all(:,j) = xsel_fit;
  width_peak_all(:,j) = width_peak;
  areas(:,j) = real( areapeak );

  plot( d_align(:,j),'k' );
  hold on;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  profile_fit = 0*x;
  for k=1:numpeaks
    predgaussian = getgaussian( x,...
				xsel_fit(k), width_peak(k),a(k)); 
    %[    xsel_fit(k), width_peak(k), a(k), sum( find( ~isreal( predgaussian)))]

    profile_fit = profile_fit + predgaussian;
    plot( predgaussian,'b');
  end
  prof_fit(:,j) = profile_fit;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  for i = 1:length( xsel_fit )
    if ( round( xsel_fit(i)) < length( profile_fit ) & round(xsel_fit(i)) > 0 )
      plot( [xsel_fit(i) xsel_fit(i)], [0 profile_fit(round(xsel_fit(i)))],'k' ); hold on;
    end
  end
  plot( prof_fit(:,j),'r');
  hold off;
  set(gca,'ylim',[-0.5 5]);

  middle_portion_min = round ( min( xsel_fit) + ( max( xsel_fit ) - min( xsel_fit ) )/3 );
  middle_portion_max = round ( min( xsel_fit) + 2*( max( xsel_fit ) - min( xsel_fit ) )/3 );
  
  ymax = max( profile_fit( middle_portion_min:middle_portion_max) );
  axis([ min(xsel_fit)-100 max(xsel_fit)+100 ...
	 -0.5 3*ymax]);

  %pause; 
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %if (count == 1 && num_xsel_lanes == 1 )
  %  xsel_start = xsel_fit;
  %  width_peak_start = width_peak;
  %end
  
end