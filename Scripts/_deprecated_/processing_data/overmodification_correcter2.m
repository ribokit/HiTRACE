function [ d_unmod, modification_ratios] = overmodification_correcter2( d, early_res, late_res )
%
%  [ d_unmod, modification_ratios] = overmodification_correcter2( d, early_res, late_res )
%

if nargin == 0;  help( mfilename ); return; end;

d_unmod = 0*d;
modification_ratios = zeros( 1, size(d,2));

% just loop through each lane.
for i = 1: size(d, 2)
  [ d_unmod( :,i ), modification_ratios(i)] = ...
      overmodification_correction_individual( d(:,i), early_res, late_res );
  d_unmod(:,i) = d_unmod(:,i) * mean( d( [early_res late_res], i ) )/ ...
      mean( d_unmod( [early_res late_res], i ) );      
end
fprintf('\n');

% make a pretty image to check what went on.
make_image( d, d_unmod );

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ d_unmod, best_modification_ratio ] = ...
    overmodification_correction_individual( d, early_res, late_res )

% "scale factor" is meant to represent the fraction of the 
%  products that is cleaved. In single-hit kinetics, this number is
%  much less than 1.
[max_modification_ratio, d_cumulative_sum ] = determine_max_modification_ratio( d );

[ d_unmod, best_modification_ratio ] = find_best_modification_ratio( ...
    d, d_cumulative_sum, early_res, late_res, max_modification_ratio  );

fprintf( 1, 'Modification ratio: %6.2f [best]   %6.2f [max]\n', ...
	 best_modification_ratio, max_modification_ratio )

%fprintf( 1, 'Modification ratio: %6.2f [best] \n',...
%	 best_modification_ratio )

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ max_modification_ratio, d_cumulative_sum ]  = determine_max_modification_ratio( d );
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
function [d_unmod, best_modification_ratio ] = ...
    find_best_modification_ratio( d, d_cumulative_sum, early_res, late_res, max_modification_ratio );
% Grid search: Scan through some possible values of "modification ratio",
%  simulate corrected data, and find the one that is most uniform.  

num_iterations = 2;

modification_ratio = 1;
d_unmod = d;
overall_modification_ratio = 0.0;

for n = 1:num_iterations

  modification_ratio = log( mean( d_unmod(early_res))/mean( d_unmod( late_res ) )  );
  modification_ratio = modification_ratio * length( d_unmod )/ ( mean(late_res)-mean(early_res));

  %if ( n == 1 ) 
  overall_modification_ratio = overall_modification_ratio + modification_ratio;
  %end
  
  overall_modification_ratio = max( overall_modification_ratio, 0.0 );
  overall_modification_ratio = min( overall_modification_ratio, max_modification_ratio );

  
  d_unmod = d_unmod .* exp ( modification_ratio .* d_cumulative_sum );

end


best_modification_ratio = overall_modification_ratio;


d_unmod = d_unmod / mean( d_unmod( [early_res, late_res] ) );

