function r_out = slice_rdat( r, search_tag, filename );
% r_out = slice_rdat( r, search_tag ,filename );
%
% Slices out a subset rdat from your rdat.
%
%  Input: 
%   r          = name of RDAT file, or RDAT object read into MATLAB.
%   search_tag = numbers of data rows that you want, or 
%                 text to look for in file.
%   filename   = [optional] output filename 
%
% (c) R. Das, Stanford University, 2013.
%

if nargin < 2; help( mfilename ); return; end;

if ischar( r )
 r = read_rdat_file( r );
end

if isnumeric( search_tag );
  data_cols = search_tag;
else
  assert( ischar( search_tag ) );
  data_cols = get_data_cols( r, search_tag ); 
end

r_out = r;
r_out.sequences  = r.sequences( data_cols );
r_out.structures = r.structures( data_cols );
r_out.data_annotations = r.data_annotations( data_cols );
r_out.reactivity = r.reactivity( :, data_cols );
if length( r.reactivity_error ) > 0
  r_out.reactivity_error = r.reactivity( :, data_cols );
end

if exist( 'filename', 'var' )
  fprintf( 'Outputting to: %s\n', filename );
  output_rdat_to_file( filename, r_out );
end
