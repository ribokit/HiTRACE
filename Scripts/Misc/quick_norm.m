function d_out = quick_norm( d, bins );

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
