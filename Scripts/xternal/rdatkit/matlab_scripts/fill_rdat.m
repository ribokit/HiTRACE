function rdat = fill_rdat( name, sequence, offset, seqpos, reactivity, structure, ...
			   annotations, data_annotations, reactivity_error, trace_in, xsel, xsel_refine, comments );
% rdat = fill_rdat( name, sequence, offset, seqpos, reactivity, structure, ...
%			   annotations, data_annotations, reactivity_error, trace, xsel, xsel_refine, comments );
%
% Copyright R. Das, P. Cordero, Stanford University, 2010,2011
%

if nargin==0; help( mfilename ); return; end;

rdat = RDATFile;
rdat.name = ''; % signal that we're not filled yet.

if ~exist( 'reactivity' );  fprintf( 'Must specify six variables: filename, name, sequence, offset, seqpos, reactivity'); end
if ~exist( 'structure' ); structure=''; end;
if ~exist( 'annotations' ); annotations = {}; end;
if ~exist( 'data_annotations' ); data_annotations = {}; end;
if ~exist( 'reactivity_error' ); reactivity_error = {}; end;
if ~exist( 'trace_in' ); trace_in = []; end;
if ~exist( 'xsel' ); xsel = []; end;
if ~exist( 'xsel_refine' ); xsel_refine = []; end;
if ~exist( 'comments' ); comments = {}; end;

if length( xsel_refine ) > 0 & ( size( reactivity, 2) ~= size( xsel_refine, 2 ) )
  xsel_refine = xsel_refine';
end


% Current version of fill_rdat script... work in progress!
rdat.version = 0.33;
rdat.comments = comments;
rdat.name = name;
rdat.sequence = sequence;
rdat.structure = structure;
rdat.offset = offset;
rdat.seqpos = seqpos;
rdat.annotations = annotations;
rdat.data_annotations = data_annotations;
rdat.reactivity = reactivity;
rdat.reactivity_error = reactivity_error;
rdat.xsel = xsel;
rdat.xsel_refine = xsel_refine;
rdat.trace = trace_in;

% reorder if seqpos is in weird order!
if length(rdat.seqpos) > 1 & rdat.seqpos(1) > rdat.seqpos(2); rdat = reverse_rdat_seqpos_order( rdat ); end

rdat = fill_sequences_and_structures( rdat );

check_rdat( rdat );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rdat = reverse_rdat_seqpos_order( rdat ); 

rdat.seqpos          = rdat.seqpos( end:-1:1 );
rdat.reactivity       = rdat.reactivity( end:-1:1, : );
rdat.reactivity_error = rdat.reactivity_error( end:-1:1, : );  
rdat.xsel            = rdat.xsel( end:-1:1 );
rdat.xsel_refine            = rdat.xsel_refine( end:-1:1, : );

