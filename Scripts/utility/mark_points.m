function [point_values, point_locations] = ...
    mark_points(imagex, ...
		mutpos, xsel, probepos,...
		point_locations );

stop_pick = 0;

if ~exist( 'point_locations' ) 
  point_locations = [];
end
count=size( point_locations,2 );

[xsize,ysize,zsize]=size(imagex);
colormap( 1- gray(100));

zoomedin = 0;
view_axes = [0 ysize 0 xsize ];

while (stop_pick<1)

  plot_stuff( imagex, xsize, ysize, ...
	      view_axes, point_locations  );
  title( 'Pick a point');
  
  [xpick,ypick,button] = ginput(1);
  switch button
   case 1
    count = count+1;
    point_locations(:,count) = [xpick;ypick];
    x = get_values( point_locations(:,count), mutpos, xsel, probepos );
    fprintf( 1, 'Selected point: %d %d\n', x(1), x(2) );
   case 2
    if size( point_locations,2) > 0 
      dist2 = ( (xpick - point_locations(1,:)).^2 + ...
		(ypick - point_locations(2,:)).^2 );
      [dummy, closestpos] = min( dist2 );

      x = get_values( point_locations(:,closestpos), mutpos, xsel, probepos );
      fprintf( 1, 'Erasing point: %d %d\n', x(1), x(2) );

      point_locations = point_locations( :, [1:(closestpos-1) ...
		    (closestpos+1):end]);
      count = count - 1;
    end
   case 3
    if (zoomedin == 0)
      view_axes = [xpick-ysize/5 xpick+ysize/5 ...
		   ypick - xsize/5 ypick+xsize/5];
      zoomedin= 1;
    else
      view_axes = [0 ysize 0 xsize ];
      zoomedin = 0;
    end       
  end
  switch char(button)
   case {'q','Q','z','Z'}
    stop_pick=1;  
  end
end
    
point_values = [];
for k = 1:count
  x = get_values( point_locations(:,k), mutpos, xsel, probepos );
  point_values(:,k) = x;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  plot_stuff( imagex, xsize, ysize, ...
		      view_axes,...
		      point_locations  );
clf;
figure(1); subplot(1,1,1);
image(imagex); hold on
axis([0 ysize 0 xsize]); 
hold on

for k = 1:size( point_locations,2)
  plot( point_locations(1,k), point_locations(2,k),'rx');
  text( point_locations(1,k), point_locations(2,k),...
	num2str(k) );
end
hold off
axis(view_axes);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  axis_locations = fix_order_axis_locations( axis_locations );
if axis_locations(2) < axis_locations(1)
  axis_locations([1 2]) = axis_locations([2 1]);
end
if axis_locations(4) < axis_locations(3)
  axis_locations([3 4]) = axis_locations([4 3]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  point_values = get_values( point_locations, mutpos, xsel, probepos );

point_values(1) =  mutpos( round( point_locations(1) ) );
[dummy, closestres] = min( abs(point_locations(2)  - xsel) );
point_values(2) = probepos( closestres );
