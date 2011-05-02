function [xsel,ysel,pickedpeak,grayscaleon, maxprof] = quickpolyxline(imx,dontplotthegel,yselpick,plottheline,grayscaleon,maxprof,ContrastScale)
% allows user to defines a polyline moving down through the gel.
%   If any number for dontplotthegel is inputed, routine doesn't
%   replot the gel.
%    Rhiju Das, August 29 2003.

pickedpeak = 0;
if (nargin<1) hold off;figure(1);image(imx);end;
if (nargin<4) grayscaleon = 1; maxprof = 1;end;

hold on;
numpixels=size(imx,1);
numxpixels=size(imx,2);
title([...
        'left                : next point of anchorline \newline',...
        'middle (opt-click)  : undo point \newline',...
        'right  (apple-click): finished']);
count=1;button=0;
if (nargin>2) xsel(1) = 0; ysel(1) = yselpick; count=2; firstsymbol=plot(1,ysel(1),'ro');end;

finishedline = 0;
while (finishedline == 0)
    [xsel(count),ysel(count),button]=ginput(1);
    if(button==2) count=count-1; count = max(count,1);
        if (count>1) set(h(count),'Visible','off');end
    end;
    
    if count==1; xsel(1)=0; plot(xsel(1),ysel(1),'ro');
    elseif (button==1)
        h(count)= plot([xsel(count-1) xsel(count)],[ysel(count-1) ysel(count)],'r'); 
        if xsel(count) == xsel(count-1) xsel(count) = xsel(count-1)+0.5;    end;
        pickedpeak=1;  count=count+1;
    end    
    if (button == 3) finishedline = 1; end;
        
    if (nargin>3)
        switch char(button)
	 case '1'
	  maxprof = maxprof/ContrastScale;
	  setcolormap(grayscaleon, maxprof);    
	 case '2'
	  maxprof = maxprof*ContrastScale;
	  setcolormap(grayscaleon, maxprof);    
	 case {'c','C'}
	  grayscaleon = 1 - grayscaleon;
	  setcolormap(grayscaleon, maxprof);    
	 case {'r','R'}
	  count = count - 1;
	  set( h( count ) , 'Visible','off');
	  range_min = max(round(ysel(count)) -100,1);
	  range_max = min(round(ysel(count)) +100, size(imx,1));
	  range_bins = [range_min:range_max];
	  [dummy, maxindex ] = max( imx( range_bins, ...
					 round( xsel(count) ) ) );
	  ysel( count ) = range_bins( maxindex );
	  h(count)= plot([xsel(count-1) xsel(count)],[ysel(count-1) ysel(count)],'r'); 
	  count = count + 1;
	end
    end
end

count=count-1; count = max(count,1);
xsel(count)=numxpixels;
lastsymbol = plot(xsel(count),ysel(count),'ro');
if (count>1) 
    set(h(count),'visible','off');
    h(count)= plot([xsel(count-1) xsel(count)],[ysel(count-1) ysel(count)],'r'); 
        if (plottheline == 0) 
        for i=1:count
            set(h(i),'visible','off');
            set(firstsymbol,'visible','off');
            set(lastsymbol,'visible','off');
        end
    end
end;
if (count<2) pickedpeak=0; set(firstsymbol,'visible','off');set(lastsymbol,'visible','off');end;

xsel=xsel(1:count);    ysel=ysel(1:count);
hold off
