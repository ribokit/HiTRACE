function rdat = output_workspace_to_rdat_file( filename, name, sequence, offset, seqpos, reactivity, ...
					       structure, ...
					       annotations, data_annotations, ...
					       reactivity_error, ...
					       trace_in, xsel, xsel_refine, comments )
%
% output_workspace_to_rdat_file( filename, name, sequence, offset, seqpos, reactivity, ...
%					       structure, ...
%					       annotations, data_annotations, ...
%					       reactivity_error, ...
%					       trace, xsel, xsel_refine, comments );
%
% Copyright R. Das, P. Cordero, Stanford University, 2010,2011
%

if nargin==0; help( mfilename ); return; end;

if exist( 'structure', 'var' ) && ~ischar(structure) && isnumeric( structure )
  fprintf( 'WARNING!WARNING!WARNING!!!\n');
  fprintf( 'NO LONGER ACCEPTING MUTPOS!\n');
  fprintf( 'WARNING!WARNING!WARNING!!!\n');
  error( 'No longer accepting mutpos'); 
  return;
end

% set defaults.
if ~exist( 'reactivity','var' );  fprintf( 'Must specify six variables: filename, name, sequence, offset, seqpos, reactivity'); end;
if ~exist( 'annotations','var' ); annotations = {}; end;
if ~exist( 'data_annotations','var' ); data_annotations = {}; end;
if ~exist( 'reactivity_error','var' ); reactivity_error = {}; end;
if ~exist( 'xsel','var' ); xsel = []; end;
if ~exist( 'xsel_refine','var' ); xsel_refine = []; end;
if ~exist( 'comments','var' ); comments = {}; end;
if ~exist( 'trace_in','var' ); trace_in = []; end;
trace = trace_in; % this is necessary because matlab has a function called trace, of course

% (t47) correct data_annotations to cell of cell instead of cell of string
for i = 1:length(data_annotations)
    if ~iscell(data_annotations{i}); data_annotations{i} = {data_annotations{i}}; end;
end;

rdat = fill_rdat( name, sequence, offset, seqpos, reactivity, structure, ...
		  annotations, data_annotations, reactivity_error, trace, xsel, xsel_refine, comments );

if isempty( rdat.name ); 
  fprintf( 'PROBLEM filling rdat!\n' );
  return;
end;

output_rdat_to_file( filename, rdat );
