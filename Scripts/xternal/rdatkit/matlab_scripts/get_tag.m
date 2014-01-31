function values = get_tag( annotations, tag )
% GET_TAG: Figure out values of specific annotations from RDAT data_annotations.
%
% value = get_tag( annotations, tag )
% value = get_tag( annotation, tag )
% value = get_tag( rdat, tag )
%
% Can also give the two arguments in reverse order.
%
% for use in parsing data_annotations or annotations cells.
%

if isobject( annotations ) & sum( strcmp( fieldnames( annotations ), 'data_annotations' ) ) > 0; annotations = annotations.data_annotations; end;
if isobject( tag ) & sum( strcmp( fieldnames( tag ), 'data_annotations' ) ) > 0; tag = tag.data_annotations; end;

if ischar( annotations ) & iscell( tag ) % user has given in opposite order.
  tag_tmp = annotations;
  annotations = tag;
  tag = tag_tmp;
end
assert( ischar( tag ) );
assert( iscell( annotations ) );

if tag(end) == ':'; tag = tag(1:end-1); end

if length( annotations ) > 0 & ischar( annotations{1} )
  values = get_tag_one_annotation( annotations, tag );
  return 
else
  for i = 1:length( annotations )
    values{i} = get_tag_one_annotation( annotations{i}, tag );
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%`
function value = get_tag_one_annotation( annotation, tag )

vals = find_annotation_tag( annotation, tag );
if length( vals ) == 0; 
  value = ''; 
else
  value = vals{1};
end