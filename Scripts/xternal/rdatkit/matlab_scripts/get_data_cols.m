function data_cols = get_data_cols( r, search_tag )
%
% data_cols = get_data_cols( r, search_tag )
%
%
fprintf( 'Searching for tag: %s\n', search_tag );

data_cols = [];

for n = 1:length( r.data_annotations );

  anot = r.data_annotations{n};

  found_it = 0;
  for m = 1:length( anot )
    if ~isempty( strfind( anot{m}, search_tag ) )
      found_it = 1;
      data_cols = [data_cols, n ];
      break;
    end
  end

end

