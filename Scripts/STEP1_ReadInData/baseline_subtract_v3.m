function d_sub = baseline_subtract_v3( d, ymin, ymax );

%Declare boundaries
if ~exist( 'ymin')
  ymin = 1500;
end
if ~exist( 'ymax')
  ymax = 4000;
end

window_size = 20;
d_smooth = smooth(d(ymin:ymax), window_size );
%r = [-10000:1:10000];

%[dummy, i ] = max( hist( d_smooth, r) );
%d_sub = d - r(i);
d_sub = d - min( floor(d_smooth) );
