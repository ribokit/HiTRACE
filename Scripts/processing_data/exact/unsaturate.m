function [ area_peak_unsaturated, area_peak_unsaturated_error, diluted_array_scaled ] = unsaturate( undiluted_array, diluted_array, undiluted_array_error, diluted_array_error, sd_cutoff, seqpos, exclude_pos )
% UNSATURATE: Corrects area_peak arrays for saturating bands.  
%
% [ area_peak_unsaturated, diluted_array_scaled ] = unsaturate( undiluted_array, diluted_array, undiluted_array_error, diluted_array_error, sd_cutoff, seqpos, exclude_pos )
%
% undiluted_array = band intensities that may have some peaks that are saturating the detector
%                    must be ordered from 5' to 3'. First position should correspond 
%                    to fully extended cDNA.
% diluted_array   = band intensities that may have lower signal-to-noise but no saturated peaks.
%                    saturated peaks in undiluted_array will be replaced by values 
%                    from this diluted_array (scaled).  Give as empty array '[]' if unknown.
% undiluted_array_error = [optional] error corresponding to undiluted_array [taken as zero if not given]
% diluted_array_error   = [optional] error corresponding to undiluted_array [taken as zero if not given]
% sd_cutoff       = [optional] standard deviation to use as a cutoff for saturation in 'unsaturate' step [default 1.5]. 
% seqpos          = [optional] sequence numbering (used for plotting) [default is 0,1,2,...]
% exclude_pos     = [optional] sequence positions to exclude from unsaturation -- just use undiluted_array values here.
%
% (c) T. Mann, 2012
% (c) T. Mann, R. Das, 2013
%

if nargin < 2;  help( mfilename ); return; end;

if ~exist( 'undiluted_array_error', 'var' ) | isempty( undiluted_array_error ); undiluted_array_error = 0 * undiluted_array; end;
if ~exist( 'diluted_array_error', 'var' ) | isempty( diluted_array_error ); diluted_array_error = 0 * diluted_array; end;
if ~exist( 'sd_cutoff' ) | isempty( sd_cutoff) | sd_cutoff == 0; sd_cutoff = 1.5; end;
if ~exist( 'seqpos' ) seqpos = [0 : size( undiluted_array, 1 ) - 1]; end;
if ~exist( 'exclude_pos' ) exclude_pos = []; end;

if ~all( size(undiluted_array) == size(diluted_array) )
  fprintf( 'The undiluted arrays do not equal the number of diluted arrays!\n');
  error('different arrays')
end;
if ~all( size(undiluted_array) == size(undiluted_array_error) )
  fprintf( 'The undiluted_array_error does not match size of undiluted_array!\n');
  error('different sizes in undiluted_array_error and undiluted_array')
end;
if ~all( size(diluted_array) == size(diluted_array_error) )
  fprintf( 'The diluted_array_error does not match size of diluted_array!\n');
  error( 'The diluted_array_error does not match size of diluted_array!\n');
end;

[num_rows, num_cols] = size(undiluted_array);

NITER = 3;

% exclude_pos may have 'conventional numbering' -- convert based on seqpos
exclude_idx = [];
for m = 1:length( exclude_pos )
  exclude_idx = [exclude_idx, find( seqpos == exclude_pos(m) )];
end
  
is_saturated_position = zeros( num_rows, num_cols );

for i = 1:num_cols;
  
  saturated_positions = [];
  for n = 1:NITER

    ok_positions = setdiff( [1:num_rows], saturated_positions );    
    scalefactor = sum( undiluted_array(ok_positions,i) ) / sum( diluted_array( ok_positions, i ) );

    is_saturated_position(:,i) = 0;
    is_saturated_position( saturated_positions, i ) = 1;
    
    diluted_array_scaled(:,i) = scalefactor * diluted_array(:,i);
    diluted_array_scaled_error(:,i) = scalefactor * diluted_array_error(:,i);
    
    residuals = undiluted_array(:,i) - diluted_array_scaled(:,i);
    stdev = std( residuals );    

    saturated_positions = find( abs(residuals) > stdev * sd_cutoff );
    %saturated_positions = find( residuals < -stdev * sd_cutoff );
    saturated_positions = setdiff( saturated_positions, exclude_idx );
  end
  
end


area_peak_unsaturated = [];
area_peak_unsaturated_error = [];

%area_peak_unsaturated gets the original value for area_peak if it was deemed
%to be non-saturating; if the peak was deemed saturating, it gets the value
%measured in diluted array scaled by the constant c (calculated to make
%other peaks overlay with undiluted array)
for i = 1:num_cols;
  for j = 1:num_rows;
    if is_saturated_position(j,i)
      area_peak_unsaturated(j,i) = diluted_array_scaled(j,i);
      area_peak_unsaturated_error(j,i) = diluted_array_scaled_error(j,i);
    else
      area_peak_unsaturated(j,i) = undiluted_array(j,i);
      area_peak_unsaturated_error(j,i) = undiluted_array_error(j,i);
    end;
  end;
end;

% some visual feedback
if ~exist( 'seqpos' ) seqpos = [0 : size( undiluted_array, 1 ) - 1]; end;

scalefactor = 40 / mean( mean( max(area_peak_unsaturated, 0 ) ) );

set(gcf, 'Name', 'Unsaturation');
set(gcf, 'Position', [0, 0, 800, 600]);
set(gcf, 'PaperOrientation', 'Landscape', 'PaperPositionMode', 'Manual', ...
    'PaperSize', [11 8.5], 'PaperPosition', [-0.65 0.15 12 8], 'Color', 'White');

subplot(1,3,1); make_image( undiluted_array, is_saturated_position, scalefactor, seqpos );
title( 'Undiluted Sample', 'FontSize', 11, 'FontWeight', 'Bold');

subplot(1,3,2); make_image( diluted_array_scaled, is_saturated_position, scalefactor,  seqpos );
title( 'Diluted Sample, Scaled', 'FontSize', 11, 'FontWeight', 'Bold');

subplot(1,3,3); make_image( area_peak_unsaturated, is_saturated_position, scalefactor, seqpos );
title( 'Unsaturated', 'FontSize', 11, 'FontWeight', 'Bold');

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



