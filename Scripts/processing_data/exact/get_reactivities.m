function [ reactivity, reactivity_error, seqpos_out, unsaturated, attenuation_corrected, reactionProb] = get_reactivities( undiluted, diluted, undiluted_error, diluted_error, bkg_col, refpos, seqpos, exclude_pos_for_unsaturation, sd_cutoff )
% GET_REACTIVITIES: Correct data for saturating bands; subtract background; normalize; and propagate errors.
%
%  [reactivity, reactivity_error, seqpos_out, ...
%   area_peak_corrected, attenuation_corrected, reactionProb] = 
%       GET_REACTIVITIES(undiluted, diluted, undiluted_error, diluted_error,  ...
%                        bkg_col, refpos, seqpos, sd_cutoff)
%
% Required inputs:
% undiluted = band intensities that may have some peaks that are saturating the detector
%                    must be ordered from 5' to 3'. First position should correspond 
%                    to fully extended cDNA.
% diluted   = band intensities that may have lower signal-to-noise but no saturated peaks.
%                    saturated peaks in undiluted will be replaced by values 
%                    from this diluted (scaled).  Give as empty array '[]' if unknown.
%
% Recommended inputs:
% undiluted_error = errors associated with undiluted array
% diluted_error   = errors associated with diluted array
% bkg_col         = If an integer, the number of the 'no mod' trace. [If 0, no background subtraction.]
%                    If you have different 'no mod' traces, you should specify
%                    a set of integers that is the same length as the total number of traces,
%                    containing the positions of the no mod that go with each trace.
% refpos          = nucleotide(s) to which you want to normalize reactivites.
%
% Optional inputs;
% seqpos          = actual sequence positions specified in undiluted.
% sd_cutoff       = standard deviation to use as a cutoff for saturation in 'unsaturate' step. 
%
%
% Outputs:
% reactivity = array of reactivities, ordered 5' to 3'. The first 
%
% Fully automated reactivity referencing workflow, starting from
% measurements of area_peak.
%
% Normalized reactivity values are returned 5' to 3'.
%
% ######
% Example of bkg_col
%
% Dataset of
%  {'nomod-1','SHAPE-1','DMS-1','CMCT-1',...
%   'nomod-2','SHAPE-2','DMS-2','CMCT-2',...
%   'nomod-3','SHAPE-3','DMS-3','CMCT-3'}
% In which -1/2/3 means different conditions and each number group should
%  be subtracted by its own nomod.
% bkg_col for this dataset should be:
%  [1,1,1,1,5,5,5,5,9,9,9,9]
%  in which 1,5,9 are the numbers of 'nomod-1/2/3' in the dataset.
% ######
%
% Thomas Mann, November 2012.
% Updated, R. Das & S. Tian, Apr - May 2013
%

if nargin == 0; help( mfilename ); return; end;

if ~exist( 'bkg_col','var' ); bkg_col = 0; end;
if ~exist( 'refpos','var' ); refpos = []; end;
if ~exist( 'seqpos','var' ); seqpos = [0 : size( undiluted, 1 ) - 1]; end;
if ~exist( 'exclude_pos_for_unsaturation','var' ); exclude_pos_for_unsaturation = []; end;
attenuation_corrected = [];
reactionProb = [];
reactivity = [];
reactionProb = {};
seqpos_out = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ( ~isempty( diluted ) &&  ~all( size( undiluted ) == size( diluted )) ); 
    error( 'Mismatch in size of undiluted and diluted.\n'); 
end;
if ( ~isempty( diluted_error ) &&  ~all( size( diluted ) == size( diluted_error )) ); 
    error( 'Mismatch in size of diluted and diluted_error.\n'); 
end;
if ( ~isempty( undiluted_error ) &&  ~all( size( undiluted ) == size( undiluted_error )) ); 
    error( 'Mismatch in size of undiluted and undiluted_error.\n'); 
end;
if ( length( bkg_col ) > 1  &&  ~( length( bkg_col) == size( diluted, 2 )) ); 
    error( 'Mismatch in size of bkg_col and data.\n'); 
end;
nres   = size( undiluted, 1);
ntrace = size( undiluted, 2 );
if ( length(seqpos) ~= nres ); 
    fprintf( 'Length of seqpos must match number of residues in undiluted.\n' ); 
    return; 
end;
if ~exist( 'sd_cutoff','var' ) sd_cutoff = 1.5; end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
unsaturated = undiluted;
UNSATURATED = 0;
if ~isempty( diluted )
  [unsaturated, unsaturated_error] = unsaturate( undiluted, diluted, undiluted_error, diluted_error, sd_cutoff, seqpos, exclude_pos_for_unsaturation);
  fprintf('Press any key to continue ...\n');
  pause;
  UNSATURATED = 1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% contains all the steps to process data after peak alignment, assignment,
% and dilution scaling.
[attenuation_corrected, attenuation_corrected_error] = correct_for_attenuation( unsaturated, unsaturated_error, seqpos );
fprintf('Press any key to continue ...\n');
pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% now do background subtraction.
reactionProb = attenuation_corrected;
reactionProb_error = attenuation_corrected_error;

image_scalefactor = 2000;
BACKGROUND_SUBTRACTED = 0;
if bkg_col(1) > 0;
  if length( bkg_col ) == 1; bkg_col = bkg_col * ones( ntrace, 1 ); end;
  [reactionProb, reactionProb_err ] = subtract_array( attenuation_corrected, attenuation_corrected(:,bkg_col), attenuation_corrected_error, attenuation_corrected_error(:,bkg_col), seqpos, image_scalefactor );
  BACKGROUND_SUBTRACTED = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normalization
ref_pos_in = [];
for k = 1:length( refpos )
  ref_pos_in = [ ref_pos_in, find( seqpos == refpos(k) ) ];
end

%reactivity = reactionProb/mean(mean(reactionProb));
reactivity       = reactionProb;
reactivity_error = reactionProb_error;
NORMALIZED = 0;
if ~isempty( ref_pos_in )
  [reactivity, norm_scalefactors, reactivity_error ] = quick_norm( reactionProb, ref_pos_in, reactionProb_error );
  image_scalefactor = 40;
  NORMALIZED = 1;
end


image( 1:ntrace, seqpos, reactivity * image_scalefactor );
title( 'Final', 'FontSize', 11, 'FontWeight', 'Bold');
if ( ntrace < 100 ) make_lines; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
REMOVED_FULL_EXTENSION_SITE = 0;
seqpos_out = seqpos(2:end);
reactivity = reactivity(2:end,:);
reactivity_error = reactivity_error(2:end,:);
REMOVED_FULL_EXTENSION_SITE = 1;

fprintf( '\n' );
if UNSATURATED;  fprintf( 'Unsaturated the data based on undiluted and diluted array.\n' );
else fprintf( 'No attempt at unsaturation.\n' ); end;

if BACKGROUND_SUBTRACTED; fprintf( 'Subtracted background.\n' );
else fprintf( 'No background subtraction.\n' ); end;

if NORMALIZED; fprintf( 'Normalized base on refpos.\n' );
else fprintf( 'No normalization -- ''absolute'' reactivities outputted.\n' ); end;

if REMOVED_FULL_EXTENSION_SITE; fprintf( 'Removed  "site 0", corresponding the fully extended cDNA [Use seqpos_out, which also removes that first data point.]\n'); end;



