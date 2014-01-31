function values = find_annotation_tag( annotation, tag )
% values = find_annotation_tag( annotation, tag )
%
% for use in parsing data_annotations or annotations cells.
%

values = {};
for m = 1:length( annotation )
  anot = annotation{m};
  idx = strfind( anot, tag );
  if ~isempty( idx ) & idx(1) == 1
    values = [values, anot( length(tag)+2: end )];
  end
end