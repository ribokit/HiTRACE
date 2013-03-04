function d_sub = baseline_subtract( d, ymin, ymax );
% BASELINE_SUBTRACT: linear baseline subtraction
%  d_sub = baseline_subtract( d, ymin, ymax );
%
% (C) R. Das 2008-2011

if nargin == 0;  help( mfilename ); return; end;

d_sub = [];

%Declare boundaries
if ~exist( 'ymin')
  ymin = 1000;
end
if ~exist( 'ymax')
  ymax = 4000;
end

window_size = 20;

for m = 1:size(d,2)
  %d_smooth = smooth(d(ymin:ymax), window_size );
  d_smooth = d(ymin:ymax, m);
  r = [-10000:1:10000];

  %[dummy, i ] = max( hist( d_smooth, r) );
  %d_sub = d - r(i);
  d_sub(:,m) = d(:,m) - mode( round(d_smooth) );
end
