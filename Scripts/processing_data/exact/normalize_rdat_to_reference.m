function r = normalize_rdat_to_reference( r, outfilename, ref_sequence, USE_LAST );
% r = normalize_rdat_to_reference( r, outfilename, ref_sequence, USE_LAST );
%
% Inputs:
%  r            = RDAT object or name of RDAT file
%  outfilename  = name of output file
%  ref_sequence = [default: 'GAGUA'] sequence of referencing segment 
%  USE_LAST     = [default: 1] use last instance, rather than all instances, of reference segment.
%
% (C) R. Das, Stanford University, 2013

if ( nargin < 2 ) help( mfilename ); return; end;

if ~exist( 'ref_sequence', 'var' ) ref_sequence = 'GAGUA'; end;
if ~exist( 'USE_LAST', 'var' ) USE_LAST = 1; end;

if ischar(r);
  rfilename = r;
  r = read_rdat_file( rfilename );
end

for i = 1:size( r.reactivity, 2 )

  sequence = r.sequences{i};
  sequence_in_r = sequence( r.seqpos - r.offset );

  modifier = get_tag( r.data_annotations{i}, 'modifier' );
  if length( modifier ) == 0;  modifier = get_tag( r.annotations, 'modifier' ); end;

  [ r.reactivity(:,i), r.reactivity_error(:,i) ] = ...
      normalize_to_reference( r.reactivity(:,i), r.reactivity_error(:,i),...
			      sequence_in_r, ref_sequence, USE_LAST, modifier );
end

output_rdat_to_file( outfilename, r );

