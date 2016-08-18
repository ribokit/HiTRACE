function [xsel_all] = refine_peak_assignment( image_x, xsel_all );

% Altered to not zoom in and out - CVL

numlanes = size(image_x,2);
contrast_factor = 1;
stop_sel = 0;

xmin = -0.5;
xmax = numlanes+0.5;
ymin = min(min( xsel_all )) - 50;
ymax = max(max( xsel_all )) + 50;


xmin_min = xmin; 
xmax_max = xmax;
ymin_min = ymin; 
ymax_max = ymax;

while ~stop_sel

  make_plot( image_x, xsel_all, xmin, xmax, ymin, ymax, contrast_factor );

  colormap( 1 - gray(100) );
  
  [ysel_pick, xsel_pick, button ]  = ginput(1);
  
  switch( button )
   case 1
    xsel_all = add_pick( xsel_all, xsel_pick, ysel_pick );
  end

  switch char(button)
   case {'q','Q','z','Z'}
    stop_sel = 1;
   case {'1'}
    contrast_factor = contrast_factor * sqrt(2);
   case {'2'}
    contrast_factor = contrast_factor / sqrt(2);
   case {'j','J'}
    %current_relative_pos =  (xselpick - ymin)/(ymax-ymin);
    current_relative_pos = 0.5;
    yscale = (ymax - ymin)*0.75;
    ymin = xsel_pick - yscale * (current_relative_pos);
    ymax = xsel_pick + yscale * ( 1- current_relative_pos);

    %xscale = (xmax - xmin)*0.75;
    %if ( xscale < ( xmax_max - xmin_min ) )
    %  xmin = ysel_pick - xscale * (current_relative_pos);
    %  xmax = ysel_pick + xscale * ( 1- current_relative_pos);
    %end
   case {'l','L'}
    %current_relative_pos =  (xselpick - ymin)/(ymax-ymin);
	current_relative_pos = 0.5;
    yscale = (ymax - ymin)/0.75;
    ymin = xsel_pick - yscale * (current_relative_pos);
    ymax = xsel_pick + yscale * ( 1- current_relative_pos);

    %xscale = (xmax - xmin)/0.75;
    %if ( xscale < ( xmax_max - xmin_min ) )
    %  xmin = ysel_pick - xscale * (current_relative_pos);
    %  xmax = ysel_pick + xscale * ( 1- current_relative_pos);
    %end
	case {'i','I'}
    yscale = (ymax - ymin);
    ymin = ymin - yscale*0.05;
    ymax = ymax - yscale*0.05;
   case {'k','K'}
    yscale = (ymax - ymin);
    ymin = ymin + yscale*0.05;
    ymax = ymax + yscale*0.05;
  end

  %if xmax > xmax_max
  %  shift = xmax - xmax_max;
  %  xmin = xmin - shift;
  %  xmax = xmax - shift;
  %end
  %if xmin < xmin_min
  %  shift = xmin_min - xmin;
  %  xmin = xmin + shift;
  %  xmax = xmax + shift;
  %end
  %if ymax > ymax_max
  %  shift = ymax - ymax_max;
  %  ymin = ymin - shift;
  %  ymax = ymax - shift;
  %end
  %if ymin < ymin_min
  %  shift = ymin_min - ymin;
  %  ymin = ymin + shift;
  %  ymax = ymax + shift;
  %end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xsel_all = add_pick( xsel_all, xsel_pick, ysel_pick )

closestlane = round( ysel_pick );
[x,closestpeak] = min( abs(xsel_all(:,closestlane) -xsel_pick ) );

hold on
plot( closestlane, xsel_all(closestpeak, closestlane),'x','color',[1 0.5 0]);
hold off

closestpeak_to_erase = closestpeak;

xsel_new = xsel_all( [1:(closestpeak_to_erase-1) ...
		    (closestpeak_to_erase+1):end], closestlane );

%%%%%%%%%%%%%%%%%%%%%%%
[ysel_pick, xsel_pick, button ]  = ginput(1);

xsel_new = sort( [xsel_new; xsel_pick] );

xsel_all( :, closestlane) = xsel_new;

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  make_plot( image_x, xsel_all, xmin, xmax, ymin, ymax, contrast_factor );

numlanes = size( image_x, 2 ) ;

clf;
image( contrast_factor * abs(image_x));

hold on
plot( xsel_all', 'r-');
%axis off
hold off;

axlim = [ xmin xmax ymin ymax];

axis( axlim );

