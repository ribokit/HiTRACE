function [point_values, axis_locations, point_locations] = mark_points_on_image(imagex, axis_values, axis_locations, point_locations);


stop_pick = 0;


if ~exist( 'point_locations' ) 
  point_locations = [];
end
count=size( point_locations,2 );

[xsize,ysize,zsize]=size(imagex);

zoomedin = 0;
view_axes = [0 ysize 0 xsize ];

if ~exist( 'axis_locations' ) | isempty(axis_locations)
  axis_locations = [];
  plot_stuff( imagex, xsize, ysize, view_axes, axis_locations, point_locations  );
  axis_locations = select_axis_locations(xsize,ysize);
end  
axis_locations = fix_order_axis_locations( axis_locations );


while (stop_pick<1)

  plot_stuff( imagex, xsize, ysize, ...
	      view_axes,...
	      axis_locations, point_locations  );
  title( 'Pick a point');
  
  [xpick,ypick,button] = ginput(1);
  switch button
   case 1
    count = count+1;
    point_locations(:,count) = [xpick;ypick];
   case 2
    if size( point_locations,2) > 0 
      dist2 = ( (xpick - point_locations(1,:)).^2 + ...
		(ypick - point_locations(2,:)).^2 );
      [dummy, closestpos] = min( dist2 );
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
   case {'x'} % reorder by x
    [ dummy, sort_index ] = sort( point_locations(1,:) );
    point_locations = point_locations( :, sort_index );
   case {'r','R'}
    % Refine axis.
    axis_dists = [ abs(xpick - axis_locations(1)) ...
		   abs(xpick - axis_locations(2)) ...
		   abs(ypick - axis_locations(3)) ...
		   abs(ypick - axis_locations(4)) ];
    [dummy,closest_axis] = min( axis_dists );

    switch closest_axis
     case {1,2}
      hold on; 
      plot( axis_locations(closest_axis)*[1 1],[0 xsize],'m');
     case {3,4}
      hold on; 
      plot( [0 ysize],axis_locations(closest_axis)*[1 1],'m');
    end

    [xpick,ypick,button] = ginput(1);
    switch closest_axis
     case {1,2}
      axis_locations(1) = xpick;
     case {3,4}
      axis_locations(4) = ypick;
    end
  
  end
end
    
point_values = [];
for k = 1:count
  point_values(1,k) =  (axis_values(2)-axis_values(1)) * ...
      (point_locations(1,k) - axis_locations(1))/ ...
      (axis_locations(2)-axis_locations(1)) + axis_values(1);

  point_values(2,k) =  (axis_values(4)-axis_values(3)) * ...
      (point_locations(2,k) - axis_locations(4))/ ...
      (axis_locations(3)-axis_locations(4)) + axis_values(3);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  axis_locations = select_axis_locations(xsize, ysize);

hold on
title( 'Select *left* y axis' )
[xpick,ypick,button] = ginput(1);
axis_locations(1) = xpick;
plot( [xpick xpick], [0, xsize],'r');

title( 'Select *right* y axis' )
[xpick,ypick,button] = ginput(1);
axis_locations(2) = xpick;
plot( [xpick xpick], [0, xsize],'r');

title( 'Select *bottom* x axis' )
[xpick,ypick,button] = ginput(1);
axis_locations(3) = ypick;
plot( [0, ysize],[ypick ypick],'r');

title( 'Select *top* x axis' )
[xpick,ypick,button] = ginput(1);
axis_locations(4) = ypick;
plot( [0, ysize],[ypick ypick],'r');

hold off

function  plot_axes( axis_locations, xsize, ysize );
hold on
plot( [axis_locations(1) axis_locations(1)], [0, xsize],'r');
plot( [axis_locations(2) axis_locations(2)], [0, xsize],'r');
plot( [0, ysize],[axis_locations(3) axis_locations(3)],'r');
plot( [0, ysize],[axis_locations(4) axis_locations(4)],'r');

hold off


function  plot_stuff( imagex, xsize, ysize, ...
		      view_axes,...
		      axis_locations, point_locations  );
clf;
figure(1); subplot(1,1,1);
image(imagex); hold on
axis([0 ysize 0 xsize]); 
hold on
if length( axis_locations) == 4
  plot_axes( axis_locations, xsize, ysize );
  set(gca,'xtick',axis_locations([1 2]),'xticklabel',{'L','R'});
  set(gca,'ytick',axis_locations([3 4]),'yticklabel',{'T','B'});
end
hold on
for k = 1:size( point_locations,2)
  plot( point_locations(1,k), point_locations(2,k),'rx');
  handle = text( point_locations(1,k), point_locations(2,k),...
		 num2str(k) );
  set(handle, 'clipping','on');
end
hold off
axis(view_axes);


function  axis_locations = fix_order_axis_locations( axis_locations );
if axis_locations(2) < axis_locations(1)
  axis_locations([1 2]) = axis_locations([2 1]);
end
if axis_locations(4) < axis_locations(3)
  axis_locations([3 4]) = axis_locations([4 3]);
end
