function  [ d_unmod, correction ] = apply_overmod_correction_EXACT( d, modification_ratio );

if ( modification_ratio == 0.0 )
  d_unmod = d;
  correction = 0*d + 1;
  return;
end

norm_factor = sum( d );

y = (d / norm_factor) * ( modification_ratio / (1 + modification_ratio ) );

attenuation = 1;
p = 0 * y;
for i = 1:length(d)
  p(i) = y(i) / attenuation;
  attenuation = attenuation * ( 1 - max( p(i), 0 ) );
end

correction = p./y;
d_unmod = d .* correction;

