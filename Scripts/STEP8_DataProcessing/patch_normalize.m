function d_norm = patch_normalize( d, reference_profile, bounds );

for j = 1 : size( d, 2 )
  for k = 1:length(bounds)

    i_max = bounds(k);
    i_min = 1;
    if ( k > 1 )
      i_min = bounds(k-1)+1;
    end
      
    d_norm( i_min:i_max, j ) = match_ref( d( i_min:i_max, j ), ...
					 reference_profile(i_min:i_max  ) );
  
  end
end

return;

function d_out = match_ref(d_in, ref )
% apply offset and scaling.
scaling = std( ref ) / std( d_in );

d_out = scaling * ( d_in - mean(d_in) )  + mean( ref );

return