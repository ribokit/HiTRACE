function [ reactivity, seqpos_out, area_peak_corrected, attenuation_corrected, reactionProb] = get_reactivities( saturated_array, diluted_array, bkg_col, refpos, seqpos, sd_cutoff )
%  [ reactivity, seqpos_out, area_peak_corrected,attenuation_corrected,reactionProb] = get_reactivities( saturated_array, diluted_array, bkg_col, refpos, seqpos, sd_cutoff )
%
% Inputs:
%
% saturated_array = band intensities that may have some peaks that are saturating the detector
%                    must be ordered from 5' to 3'. First position should correspond to fully extended cDNA.
% diluted_array   = band intensities that may have lower signal-to-noise but no saturated peaks.
%                   saturated peakes in saturated_array will be replaced by values from this diluted_array (scaled).
% bkg_col         = If an integer, the number of the 'no mod' trace. [If 0, no background subtraction.]
%                    If you have different 'no mod' traces, you should specify
%                    a set of integers that is the same length as the total number of traces,
%                    containing the positions of the no mod that go with each trace.
% refpos          = nucleotide(s) to which you want to normalize reactivites.
% seqpos          = actual sequence positions specified in saturated_array.
% sd_cutoff       = [optional] standard deviation to use as a cutoff for saturation in 'unsaturate' step. 
%
%
% Outputs:
%
% reactivity = array of reactivities, ordered 5' to 3'. The first 
%
% Fully automated reactivity referencing workflow, starting from
% measurements of area_peak.  See "unsaturate" for more details on
% saturated_array, unsaturated_array, and sd_cutoff.
%
% Normalized reactivity values are returned 5' to 3'.
%
% Thomas Mann, November 2012.

if nargin == 0; help( mfilename ); return; end;

if ~exist( 'bkg_col' ) bkg_col = 0; end;
if ~exist( 'refpos' ) refpos = []; end;
if ~exist( 'seqpos' ) seqpos = [0 : size( saturated_array, 1 ) - 1]; end;

area_peak_corrected = [];
attenuation_corrected = [];
reactionProb = [];
reactivity = [];
reactionProb = {};
seqpos_out = [];

if ( ~all( size( saturated_array ) == size( diluted_array )) ); fprintf( 'mismatch in size of saturated_array and diluted_array'); return;  end;
nres   = size( saturated_array, 1);
ntrace = size( saturated_array, 2 );
if ( length(seqpos) ~= nres ); fprintf( 'Length of seqpos must match number of residues in saturated_array' ); return; end;

if ~exist( 'sd_cutoff' ) sd_cutoff = 1.5; end;
area_peak_corrected = unsaturate( saturated_array, diluted_array, sd_cutoff, seqpos);
pause;


% contains all the steps to process data after peak alignment, assignment,
% and dilution scaling.
attenuation_corrected = correct_for_attenuation( area_peak_corrected, seqpos );
pause;

% now do background subtraction.
reactionProb = attenuation_corrected;

if bkg_col(1) > 0;

  if length( bkg_col ) == 1; bkg_col = bkg_col * ones( ntrace, 1 ); end;
  reactionProb  = attenuation_corrected - attenuation_corrected(:,bkg_col);
    
  image( 1:ntrace, seqpos, reactionProb * 2000 );
  title( 'background corrected');
  make_lines;
  pause;
  
end

ref_pos_in_array = [];
for k = 1:length( refpos )
  ref_pos_in_array = [ ref_pos_in_array, find( seqpos == refpos(k) ) ];
end

reactivity = reactionProb/mean(mean(reactionProb));
if length( ref_pos_in_array ) > 0
  reactivity = quick_norm( reactionProb, ref_pos_in_array );
end

image( 1:ntrace, seqpos, reactivity * 40 );
title( 'final');
make_lines;

seqpos_out = seqpos(2:end);
reactivity = reactivity(2:end,:);
fprintf( 'NOTE! Reactivity does not have  "site 0", the fully extended cDNA.\nUse seqpos_out, which also removes that first data point.\n');



