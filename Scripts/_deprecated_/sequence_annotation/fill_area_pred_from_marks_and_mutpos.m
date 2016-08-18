function       area_pred = fill_area_pred_from_marks_and_mutpos( marks, mutpos, seqpos, offset );

area_pred = zeros( length( seqpos ), length( mutpos ) );

for i = 1:length( mutpos )
  if ~isnan( mutpos( i ) )
    gp = find( mutpos(i) == marks( :,1 ) ) ;
    if length( gp ) > 0
      for m = gp'
	area_pred( find( seqpos == marks(m,2) ), i ) = 1.0;
      end
    end
  end
end