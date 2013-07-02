function [d_out, m, s] = filter_outliers( d, stdev_cutoff )
%  d_out = filter_outliers( d, stdev_cutoff )
%

if ~exist( 'stdev_cutoff' ); stdev_cutoff = 5; end;
N_ITER = 5;
d_out = d;
for n = 1:N_ITER
  s = std( d_out );
  m = mean( d_out );
  gp = find( abs(d_out - m) / s < stdev_cutoff );
  d_out = d_out( gp );
end
