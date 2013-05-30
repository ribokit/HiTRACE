function d_separate = separate_data( data, bounds );

for i = 1:length( bounds )
  minpos = 1;
  if i >1; minpos = bounds(i-1) + 1; end;
  maxpos = bounds(i);
  
  d_separate{i} = data(minpos:maxpos,:);
end