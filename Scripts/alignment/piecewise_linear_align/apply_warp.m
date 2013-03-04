function d_out = apply_transform( d_all, x_transform_all, PLOT_STUFF );
%
%  d_out = apply_transform( d_all, x_transform_all, PLOT_STUFF );
%

if nargin == 0;  help( mfilename ); return; end;

if ~exist( 'PLOT_STUFF' ) PLOT_STUFF = 1; end;

which_lanes =  1:size( d_all, 2 );
num_pixels = size( d_all, 1 );

for i = which_lanes;  
  d_out(:,i) = interp1( x_transform_all(:,i)', d_all( 1: num_pixels, i), [1:num_pixels], 'linear',0.0);
end

if PLOT_STUFF
  colormap( 1 - gray(100 ) );
  subplot(1,2,1);
  scalefactor = 40 / mean(mean(d_all));
  image( scalefactor * d_all )
  

  subplot(1,2,2);
  image( scalefactor * d_out )
end