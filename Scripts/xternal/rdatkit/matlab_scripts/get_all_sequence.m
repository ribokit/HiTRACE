function all_sequence = get_all_sequence( d );

for i = 1:length( d.data_annotations )
  anot = d.data_annotations{i};
  for j = 1:length( anot )
    if ~isempty( strfind( anot{j}, 'sequence:' ) )
      [dummy, r ] = strtok( anot{j}, 'sequence:' );
      all_sequence{i} = dummy;
    end
  end
end
