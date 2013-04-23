function [ reactivity, seqpos_out, area_peak_corrected, attenuation_corrected, reactionProb] = get_reactivities( undiluted_array, diluted_array, bkg_col, refpos, seqpos, exclude_pos_for_unsaturation, sd_cutoff )
%  [ reactivity, seqpos_out, area_peak_corrected,attenuation_corrected,reactionProb] = get_reactivities( undiluted_array, diluted_array, bkg_col, refpos, seqpos, sd_cutoff )
%
% Inputs:
%
% undiluted_array = band intensities that may have some peaks that are saturating the detector
%                    must be ordered from 5' to 3'. First position should correspond 
%                    to fully extended cDNA.
% diluted_array   = band intensities that may have lower signal-to-noise but no saturated peaks.
%                    saturated peaks in undiluted_array will be replaced by values 
%                    from this diluted_array (scaled).  Give as empty array '[]' if unknown.
% bkg_col         = If an integer, the number of the 'no mod' trace. [If 0, no background subtraction.]
%                    If you have different 'no mod' traces, you should specify
%                    a set of integers that is the same length as the total number of traces,
%                    containing the positions of the no mod that go with each trace.
% refpos          = nucleotide(s) to which you want to normalize reactivites.
% seqpos          = actual sequence positions specified in undiluted_array.
% sd_cutoff       = [optional] standard deviation to use as a cutoff for saturation in 'unsaturate' step. 
%
%
% Outputs:
%
% reactivity = array of reactivities, ordered 5' to 3'. The first 
%
% Fully automated reactivity referencing workflow, starting from
% measurements of area_peak.  See "unsaturate" for more details on
% undiluted_array, unsaturated_array, and sd_cutoff.
%
% Normalized reactivity values are returned 5' to 3'.
%
% Thomas Mann, November 2012.
% Updated, R. Das & S. Tian, 2013

if nargin == 0; help( mfilename ); return; end;

if ~exist( 'bkg_col' ) bkg_col = 0; end;
if ~exist( 'refpos' ) refpos = []; end;
if ~exist( 'seqpos' ) seqpos = [0 : size( undiluted_array, 1 ) - 1]; end;
if ~exist( 'exclude_pos_for_unsaturation' ); exclude_pos_for_unsaturation = []; end;
area_peak_corrected = [];
attenuation_corrected = [];
reactionProb = [];
reactivity = [];
reactionProb = {};
seqpos_out = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ( ~isempty( diluted_array ) &  ~all( size( undiluted_array ) == size( diluted_array )) ); fprintf( 'mismatch in size of undiluted_array and diluted_array'); return;  end;
nres   = size( undiluted_array, 1);
ntrace = size( undiluted_array, 2 );
if ( length(seqpos) ~= nres ); fprintf( 'Length of seqpos must match number of residues in undiluted_array' ); return; end;
if ~exist( 'sd_cutoff' ) sd_cutoff = 1.5; end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
unsaturated_array = undiluted_array;
UNSATURATED = 0;
if ~isempty( diluted_array )
  unsaturated_array = unsaturate( undiluted_array, diluted_array, sd_cutoff, seqpos, exclude_pos_for_unsaturation);
  pause;
  UNSATURATED = 1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% contains all the steps to process data after peak alignment, assignment,
% and dilution scaling.
attenuation_corrected = correct_for_attenuation( unsaturated_array, seqpos );
pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% now do background subtraction.
reactionProb = attenuation_corrected;

image_scalefactor = 2000;
BACKGROUND_SUBTRACTED = 0;
if bkg_col(1) > 0;

  if length( bkg_col ) == 1; bkg_col = bkg_col * ones( ntrace, 1 ); end;
  reactionProb  = attenuation_corrected - attenuation_corrected(:,bkg_col);
    
  image( 1:ntrace, seqpos, reactionProb * image_scalefactor );
  title( 'background corrected');
  make_lines;
  pause;

  BACKGROUND_SUBTRACTED = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normalization
ref_pos_in_array = [];
for k = 1:length( refpos )
  ref_pos_in_array = [ ref_pos_in_array, find( seqpos == refpos(k) ) ];
end

%reactivity = reactionProb/mean(mean(reactionProb));
reactivity = reactionProb;
NORMALIZED = 0;
if length( ref_pos_in_array ) > 0
  reactivity = quick_norm( reactionProb, ref_pos_in_array );
  image_scalefactor = 40;
  NORMALIZED = 1;
end

image( 1:ntrace, seqpos, reactivity * image_scalefactor );
title( 'final');
make_lines;

seqpos_out = seqpos(    2:end);
reactivity = reactivity(2:end,:);

fprintf( '\n' );
if UNSATURATED  fprintf( 'Unsaturated the data based on undiluted and diluted array.\n' );
else fprintf( 'No attempt at unsaturation.\n' ); end;

if BACKGROUND_SUBTRACTED fprintf( 'Subtracted background.\n' );
else fprintf( 'No background subtraction.\n' ); end;

if NORMALIZED fprintf( 'Normalized base on refpos.\n' );
else fprintf( 'No normalization -- ''absolute'' reactivities outputted.\n' ); end;

fprintf( 'Removed  "site 0", corresponding the fully extended cDNA [Use seqpos_out, which also removes that first data point.]\n');



