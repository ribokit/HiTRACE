function [ area_peak_unsaturated, diluted_array_scaled ] = unsaturate( saturated_array, diluted_array, sd_cutoff, seqpos )
% UNSATURATE: Corrects area_peak arrays for saturating bands.  
%
% [ area_peak_unsaturated, diluted_array_scaled ] = unsaturate( saturated_array, diluted_array, sd_cutoff )
%
%Saturated array is the set of arrays measured at a full concentration; 
%diluted array is the same samples, in the same order, run at diluted version of the final sample.
%
%The diluted samples allow for saturating band quantitation.  REQUIRES THAT
%THE TWO SETS OF ARRAYS BE IN THE SAME ORDER (e.g. column 1 of
%saturated_array corresponds to the same sample as is in column 1 of
%diluted_array.  ALSO REQUIRES THAT AREA_PEAK ARRAYS BE LISTED WITH n ROWS IN 1
%COLUMN.
%
% sd_cutoff adjusts the number of std. devs a residual has to be from the
% mean to be rejected as a saturated band.  A higher sd_cutoff requires
% more deviation from the mean residual for a band to be excluded as an
% outlier. [default: 1.5]
%
% (c) T. Mann, 2012
% (c) T. Mann, R. Das, 2013
%

if nargin == 0;  help( mfilename ); return; end;

if ~exist( 'sd_cutoff' ) sd_cutoff = 1.5; end;

if length(saturated_array) ~= length(diluted_array)
    error('different arrays','The saturated arrays do not equal the number of diluted arrays!')
end;

[num_rows, num_cols] = size(saturated_array);

NITER = 3;



is_saturated_position = zeros( num_rows, num_cols );

for i = 1:num_cols;
  
  saturated_positions = [];
  for n = 1:NITER

    ok_positions = setdiff( [1:num_rows], saturated_positions );
    scalefactor = sum( saturated_array(ok_positions,i) ) / sum( diluted_array( ok_positions, i ) );

    is_saturated_position(:,i) = 0;
    is_saturated_position( saturated_positions, i ) = 1;
    
    diluted_array_scaled(:,i) = scalefactor * diluted_array(:,i);
    
    residuals = saturated_array(:,i) - diluted_array(:,i);
    stdev = std( residuals );
    saturated_positions = find( abs(residuals) > stdev * sd_cutoff );
  end
  
end

area_peak_unsaturated = [];

%area_peak_unsaturated gets the original value for area_peak if it was deemed
%to be non-saturating; if the peak was deemed saturating, it gets the value
%measured in diluted array scaled by the constant c (calculated to make
%other peaks overlay with saturated array)
for i = 1:num_cols;
  for j = 1:num_rows;
    if is_saturated_position(j,i)
      area_peak_unsaturated(j,i) = diluted_array_scaled(j,i);
    else
      area_peak_unsaturated(j,i) = saturated_array(j,i);
    end;
  end;
end;

% some visual feedback
if ~exist( 'seqpos' ) seqpos = [0 : size( saturated_array, 1 ) - 1]; end;

scalefactor = 40 / mean( mean( saturated_array ) );

subplot(1,3,1); make_image( saturated_array, is_saturated_position, scalefactor, seqpos );
title( 'concentrated sample');

subplot(1,3,2); make_image( diluted_array_scaled, is_saturated_position, scalefactor,  seqpos );
title( 'diluted sample, scaled');

subplot(1,3,3); make_image( area_peak_unsaturated, is_saturated_position, scalefactor, seqpos );
title( 'unsaturated');

return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function make_image( image_array, is_saturated_position, scalefactor, seqpos );

image( 1:size(image_array,2) , seqpos,  scalefactor * image_array );

for j = 1:size( image_array,2)
  hold on; plot( [j j]+0.5,  [0.5 size(image_array,1)+0.5], 'k' );
end

box_saturated_positions( is_saturated_position, seqpos );

make_lines;

hold off
colormap( 1 - gray(100));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% put red box at places where information from diluted sample replaced
% information from concentrated sample -- the saturated positions
function box_saturated_positions( is_saturated_position, seqpos );

for i = 1:size( is_saturated_position, 1 );
  for j = 1:size( is_saturated_position, 2 );

    if is_saturated_position(i,j)
      h=rectangle( 'Position', [j-0.5, seqpos(i)-0.5,1,1]);
      set(h,'edgecolor','r' );
    
    end
  end
end



