function [ d_unmod, modification_ratios] = overmodification_correcter( d )
%
%  [ d_unmod, modification_ratios] = overmodification_correcter( d )
%
%
%

if nargin == 0;  help( mfilename ); return; end;

d_unmod = 0*d;
modification_ratios = zeros( 1, size(d,2));

% just loop through each lane.
for i = 1: size(d, 2)
  [ d_unmod( :,i ), modification_ratios(i)] = overmodification_correction_individual( d(:,i) );
end

% make a pretty image to check what went on.
make_image( d, d_unmod );

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ d_unmod, best_modification_ratio ] = overmodification_correction_individual( d )

% "scale factor" is meant to represent the fraction of the 
%  products that is cleaved. In single-hit kinetics, this number is
%  much less than 1.
[max_modification_ratio, d_cumulative_sum ] = determine_max_modification_ratio( d );

[ d_unmod, best_modification_ratio ] = find_best_modification_ratio( d, d_cumulative_sum, max_modification_ratio );

fprintf( 1, 'Modification ratio: %6.2f [best]   %6.2f [max]\n', ...
	 best_modification_ratio, max_modification_ratio )

%fprintf( 1, 'Modification ratio: %6.2f [best] \n',...
%	 best_modification_ratio )

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ max_modification_ratio, d_cumulative_sum ]  = determine_max_modification_ratio( d )
% One way to estimate the scale factor is the ratio of 
% cleave products to uncleaved products.

% Don't read the first few residues -- dominated by fluorescence
% from "excess" primers that were unextended.
eat_in = 4;
d_cumulative_sum = 0 * d;

% d_cumulative_sum integrates band intensity from smaller to larger DNAs.
d_cumulative_sum( eat_in:end ) = cumsum( d(eat_in:end) );
d_cumulative_sum = d_cumulative_sum/d_cumulative_sum( end );

cleaved_fraction = d_cumulative_sum( end-eat_in );
max_modification_ratio = cleaved_fraction / ( 1- cleaved_fraction );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [d_unmod, best_modification_ratio ] = find_best_modification_ratio( d, d_cumulative_sum, max_modification_ratio )
% Grid search: Scan through some possible values of "modification ratio",
%  simulate corrected data, and find the one that is most uniform.  

UNCLEAVED_BAND_WEIGHT = 0.0;

modification_ratio = [0.0: 0.1: 5.0];

%modification_ratio = modification_ratio(  ...
%    find( modification_ratio < max_modification_ratio ) );

% Trying to find modification_ratio that gives the leat variance
% across the profile. 
norm_index = get_norm_index( d );

% Correction for over-modification
for i = 1:length( modification_ratio );
  d_unmod(:,i) = d .* exp ( modification_ratio(i) .* d_cumulative_sum );
  d_unmod(:,i) =  d_unmod(:,i) * mean( d(norm_index) )/ ...
      mean( d_unmod(norm_index,i) );

  x = [0:10:200];
  h(:,i) = hist( d_unmod(norm_index,i), x);

  dx = d_unmod( norm_index, i );
  dx_sort = sort( dx );

  % "half-tail" standard deviation-- no outliers with very strong signals!
  %    an alternative might be to fit a gamma function or something.
  midpt = round( (length(dx)/2)+1 );
  uniformity( i ) = mean( ( -dx_sort(1:midpt ) + dx_sort(midpt))); 
  
end

favor_low_modification_ratios = 1.2 * modification_ratio;
UNCLEAVED_BAND_WEIGHT = 1.0;
return_uncleaved_band_to_unity = UNCLEAVED_BAND_WEIGHT * abs((max_modification_ratio - modification_ratio)/3).^2;
cost_function = uniformity + favor_low_modification_ratios + return_uncleaved_band_to_unity;

PLOT_STUFF = 0;
if PLOT_STUFF
  subplot(3,1,1); plot( d_unmod( norm_index,:) );
  %subplot(3,1,2); plot( modification_ratio, uniformity );
  subplot(3,1,2); plot( x, h );
  subplot(3,1,3); plot( modification_ratio, cost_function );
  %pause;
  figure()
end

[best_variance, best_index ]  = min( cost_function );

best_modification_ratio = modification_ratio( best_index );
d_unmod = d_unmod( :, best_index );


function norm_index = get_norm_index( d )

%What bins to look at? Don't look at edges
nres = size( d, 1 );
eat_in = 4;
eat_in = min( (nres/2 - eat_in), eat_in );
norm_index = [ eat_in: nres - eat_in];

