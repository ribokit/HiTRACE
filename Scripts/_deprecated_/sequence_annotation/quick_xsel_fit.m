function xsel = quick_xsel_fit( profiles, xsel );

stop_sel = 0;

normbins = [2300:2700];

while ( ~stop_sel )

  clf;
  plot( 0.2 * profiles/ mean( profiles( normbins,1) ) ) ;
  hold on;
  for i = 1:length( xsel )
    plot( [xsel(i) xsel(i)], [0 10],'r' ); 
  end
  hold off;
  axis([ 1500  4500 0 10]);

  drawnow;

  [xselpick, yselpick, button ]  = ginput(1);
  
  drawnow;

  switch( button )
   case 1
    xsel = [xsel xselpick];
    xsel = sort( xsel );
   case 2
    xsel = remove_pick( xsel, xselpick );
  end
  
  switch char(button)
   case {'q','Q','q','Z'}
    stop_sel = 1;
   case {'e','E'}
    xsel = remove_pick( xsel, xselpick );
   case {'j','J'}
    xsel = xsel -1;
   case {'l','L'}
    xsel = xsel +1;
   case {'i','i'}
    xsel = 1.005*(xsel-min(xsel))+min(xsel);
   case {'k','K'}
    xsel = (1/1.005)*(xsel-min(xsel))+min(xsel);
  end
  
  drawnow;
  
end


function xsel = remove_pick( xsel, xselpick )

if length( xsel ) > 0
  [dummy, closestpick] = min( abs(xsel - xselpick) );
  xsel = xsel( [1:(closestpick-1) ...
		(closestpick+1):length(xsel)] ...
	       );
end



