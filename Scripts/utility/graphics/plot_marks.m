function plot_marks( mutpos, seqpos, marks )
%  plot_marks( mutpos, seqpos, marks )
if nargin == 0;  help( mfilename ); return; end;

hold on
for j = 1: size( mutpos, 2 );
  gp = find( mutpos(j) == marks(:,1) );
  if ( length( gp ) > 0 )
    for m = gp'
      yp = find( seqpos == marks(m, 2) );
      if ( length( yp ) > 0 )
	for n = yp
	  plot( j, n, 'x','markersize',5 ); 
	end
      end
    end
  end
end
hold off
