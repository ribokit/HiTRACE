function [d_out, norm_factor ] = quick_norm( d, bins );
% QUICK_NORM: normalize data based on mean in specified bins
%
% [d_out, norm_factor ] = quick_norm( d, bins );
%
% (C) R. Das, 2008.
%

if nargin == 0;  help( mfilename ); return; end;

do_transpose = 0;
if size( d, 1) == 1 && size( d, 2)>1
  d = d';
  do_transpose = 1;
end

if ~exist('bins')
  bins = 1:size(d, 1);
end

%for i = 1:size( d, 2 )
%  if ( min( d( bins,i ))< 0  ) 
%    d(:,i) = d(:,i) - min( d(bins,i));
%  end
%end

norm_factor = 1.0 ./ mean( d( bins,:),1);

d_out = d * diag( norm_factor );

if do_transpose
  d_out = d_out';
end
