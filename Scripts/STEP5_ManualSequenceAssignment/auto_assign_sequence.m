function [xsel, D] = auto_assign_sequence( image_x,sequence, seqpos, offset, area_pred, ideal_spacing, input_bounds, PLOT_STUFF );
% AUTO_ASSIGN_SEQUENCE: (still experimental) automatic assignment of bands, given expected locations of marks
%
%
% (C) R. Das, 2010-2011
%
if ~exist( 'PLOT_STUFF' ); PLOT_STUFF = 1; end;
if ~exist( 'ideal_spacing' ); ideal_spacing = 12; end;
if ~exist( 'input_bounds' ) input_bounds = []; end;
num_lanes = size( image_x, 2 );
nres = length( sequence );

for i = 1:num_lanes
  %peaks = localMaximum( image_x(:,i) );
  %vals = image_x( peaks, i );
  %[vals_norm, scalefactor ] = SHAPE_normalize( vals );
  scalefactor = mean( image_x(:,i) );
  image_x(:,i) = image_x(:,i)/scalefactor/2;
end

if PLOT_STUFF
  %after normalization
  subplot(1,2,1);
  image( image_x*50 ) 
end


s = area_pred;
%s = zeros( length(sequence)+1, num_lanes );
%for i = 1:num_lanes
%  gp = find( marks(:,1) == mutpos(i) );
%  for m = gp';
%    s( find( seqpos == marks(m,2) ), i ) = 1.0;
%  end
%end

s(1,:)      = 10;
s(nres,:)   =  5;
s(nres+1,:) = 10;
s( find( s == 0.0 ) ) = 0.1;
sequence_at_bands = [sequence( end:-1:1 ), 'X'];

if PLOT_STUFF
  subplot(1,2,2);
  image( s*50 );
  colormap( 1 - gray(100) );
end

[xsel_fit, D]  = solve_xsel_by_DP( image_x, s, sequence_at_bands, ideal_spacing, input_bounds );
xsel = xsel_fit(1:end-1);
