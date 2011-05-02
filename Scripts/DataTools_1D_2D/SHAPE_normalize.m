function [d_norm, scalefactor, cap_value ] = SHAPE_normalize( d_for_scalefactor );

d_OK = d_for_scalefactor( find( ~isnan( d_for_scalefactor) ) );

dsort = sort( d_OK );

% this attempts to recover the normalization scheme in ShapeFinder, as reported by Deigan et al., 2008.
interquartile_range = abs( dsort( round( 0.25*length(dsort) ) ) - dsort( round( 0.75*length(dsort) ) ) );
outlier_cutoff = min( find( dsort > 1.5 ) );
%[outlier_cutoff length( dsort ) ]
if ~isempty( outlier_cutoff )
  actual_cutoff = outlier_cutoff - 1;
  %if length( dsort ) < 100; actual_cutoff = max( actual_cutoff, round(0.95*length(dsort) ) - 1 ); end;
  actual_cutoff = max( actual_cutoff, round(0.95*length(dsort) ) - 1 ); 
  cap_value = dsort( actual_cutoff+1 );
  dsort = dsort(1: actual_cutoff );
end
scalefactor = mean(dsort(  round( 0.9 * length(dsort)):end ) );

d_norm = d_for_scalefactor / scalefactor;