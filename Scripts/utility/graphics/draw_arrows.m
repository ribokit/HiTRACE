function draw_arrows( point_values_SHAPE, residue_locations, offset, ...
		      colorcode, RADIUS );
%  draw_arrows( point_values_SHAPE, residue_locations, offset,  colorcode, RADIUS );

startpos = unique(point_values_SHAPE(1,:));
for k = startpos
  gp = find( k == point_values_SHAPE(1,:)) ;

  seqpos = sort( point_values_SHAPE(2,gp) ); 
  seqpos_plot = seqpos(1);
  for i = seqpos(2:end)
    if ( i > seqpos_plot(end) + 1 )
      draw_arrow( residue_locations(:, k-offset),...
		  residue_locations(:, seqpos_plot-offset),...
		  colorcode...
		  );
      seqpos_plot = [i];
    else
      seqpos_plot = [ seqpos_plot i];
    end
  end
  draw_arrow( residue_locations(:, k-offset),...
	      residue_locations(:, seqpos_plot-offset),...
	      colorcode, RADIUS...
	      );
end
