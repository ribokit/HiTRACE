function [d_out, norm_factor, d_out_err ] = quick_norm( d, bins, d_err, FORCE_POSITIVE )
% QUICK_NORM: normalize data based on mean in specified bins
%
% [d_out, norm_factor ] = quick_norm( d, bins );
%
% (C) R. Das, 2008.
%

if nargin == 0;  help( mfilename ); return; end;

if ~exist( 'd_err','var') d_err = 0*d; end;
if ~exist( 'FORCE_POSITIVE','var') FORCE_POSITIVE = 1; end;

do_transpose = 0;
if size( d, 1) == 1 && size( d, 2)>1
  d = d';
  do_transpose = 1;
end

if ~exist('bins','var') || isempty( bins )
  bins = 1:size(d, 1);
end

%for i = 1:size( d, 2 )
%  if ( min( d( bins,i ))< 0  ) 
%    d(:,i) = d(:,i) - min( d(bins,i));
%  end
%end

for i = 1:size(d,2)
  if FORCE_POSITIVE
    norm_factor( i ) = 1.0 ./ mean( max( d( bins,i), 0) );
  else
    norm_factor( i ) = 1.0 ./ mean( d( bins,i) );
  end
  
  if isnan( norm_factor(i) ) | isinf( norm_factor(i) ); norm_factor(i) = 0.0; end;
  
  d_out(:,i) = d(:,i) * norm_factor( i );
  d_out_err(:,i) = d_err(:,i) * norm_factor( i );

  if ~isempty( find(isnan( d_out(:,i) )) )
    i
    norm_factor(i)
  end
  
end

if do_transpose
  d_out = d_out';
  d_out_err = d_out_err';
end
