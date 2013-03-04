function draw_arrow( x1, x2, colorcode, RADIUS );
%  draw_arrow( x1, x2, colorcode, RADIUS );

if ~exist('RADIUS')
  RADIUS = 15;
end

colorcode_fade = get_colorcode_fade( colorcode, size( x2, 2) );

h=rectangle( 'Position', [x1(1)-RADIUS x1(2)-RADIUS 2*RADIUS 2*RADIUS],...
	     'Curvature',[1 1]);
set(h,'edgecolor',colorcode_fade,'linewidth',1.5 );

startrectx = min( x2(1,:) );
startrecty = min( x2(2,:) );
maxrectx = max( x2(1,:) );
maxrecty = max( x2(2,:) );

delx = [maxrectx - startrectx]+2*RADIUS;
dely = [maxrecty - startrecty]+2*RADIUS;

h=rectangle( 'Position', ...
	     [startrectx-RADIUS startrecty-RADIUS delx dely ],...
	     'Curvature',[0.5 0.5]);
set(h,'edgecolor',colorcode_fade,'linewidth',1.5 );

meanx = mean( x2(1,:)); 
meany = mean( x2(2,:)); 
x2 = [meanx; meany];

x1_shift = x1 + RADIUS*(x2-x1)/norm( x2-x1);
x2_shift = x2 - RADIUS*(x2-x1)/norm( x2-x1);
h = arrow( x1_shift, x2_shift,5 );

set(h,'facecolor',colorcode_fade,'edgecolor',colorcode_fade,'linewidth',1.5);

function colorcode_fade = get_colorcode_fade( colorcode, numpts );

fade_factor = exp( - 0.25*( numpts-1) );

colorcode_fade = [1.0,1.0,1.0] - fade_factor * ([1.0, 1.0,1.0] - colorcode );