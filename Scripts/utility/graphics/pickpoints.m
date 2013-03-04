function residue_locations = pickpoints(imagex,offset,residue_locations,square_width,GRIDSIZE)
% residue_locations = pickpoints(imagex,offset,residue_locations,square_width,GRIDSIZE)
%
% imagex            = RGB image [M1 x M2 x 3 matrix] read in from, say a tif file with 
%                      the 'imread' command.
% offset            = integer to add to 1, 2, ... N to get the actual positions on your image
% residue_locations = 2 x N matrix with the (x,y) positions of each 'sequence' position 
%                      on the image. Initially you can set this to the empty set [].
% squarewidth       = [default 24] how big to make squares.
% GRIDSIZE          = [default 1] 'snap' squares to grid with this pixel resolution. 
%
%

if nargin == 0;  help( mfilename ); return; end;

figure(1); subplot(1,1,1); hold off; image(imagex); hold on
[xsize,ysize,zsize]=size(imagex);
axis([0 ysize 0 xsize]); zoomedin = 0;

if ~exist('square_width')
  square_width = 40;
end

stop_pick = 0;
count=1;
if ~exist( 'GRIDSIZE'); GRIDSIZE=1; end;

if (nargin>2)
    count = length(residue_locations)+1;
    for k=1:count-1
        xpick = residue_locations(1,k);
        ypick = residue_locations(2,k);
        h(k) = rectangle('Position',...
            [xpick - square_width/2, ypick-square_width/2,...
                square_width,square_width]);
        set(h(k),'edgecolor','b');
    end
end  

set(gcf,'color','w');
axis off

while (stop_pick<1)
    title(['left-click, select;  middle-click, replace any box; \newline i,j,k,l, adjust last box; q, done; next: ', num2str(count+offset)])
    [xpick,ypick,button] = ginput(1);
    switch button
        case 1
            xpick = xpick - mod(xpick,GRIDSIZE);
            ypick = ypick - mod(ypick,GRIDSIZE);
            residue_locations(:,count) = [xpick;ypick];
            h(count) = rectangle('Position',...
                [xpick - square_width/2, ypick-square_width/2,...
                square_width,square_width]);
            set(h(count),'edgecolor','r');
            count= count+1;
            title(['next: ', num2str(count+offset)])
        case 2
	 %count= count-1;
	 %set(h(count),'visible','off');
	    %case {'e','E'}
            [dummy, erasesquare] = min( (residue_locations(1,:) - xpick).^2 + (residue_locations(2,:) - ypick).^2 );
            if (erasesquare>0)
	      %get(h(erasesquare))
                set(h(erasesquare),'visible','off');
                title(['Replace Pick ', num2str(erasesquare+offset), ' Now'])
                [xpick,ypick,button] = ginput(1);
                xpick = xpick - mod(xpick,GRIDSIZE);
                ypick = ypick - mod(ypick,GRIDSIZE);
                residue_locations(:,erasesquare) = [xpick;ypick];
                h(erasesquare) = rectangle('Position',...
                    [xpick - square_width/2, ypick-square_width/2,...
                    square_width,square_width]);
                set(h(erasesquare),'edgecolor','r');
            end
        case 3
            if (zoomedin ==0)
                axis([xpick-ysize/5 xpick+ysize/5 ...
                    ypick - xsize/5 ypick+xsize/5])
                zoomedin= 1;
            else
                axis([0 ysize 0 xsize]); zoomedin = 0;

            end
    end

    switch char(button)
        case {'q','z','Q','Z'}
            stop_pick=1;
        case {'i','w','I','W'}
            residue_locations(2,count-1) = residue_locations(2,count-1) - GRIDSIZE;
        case {'k','s','K','S'}
            residue_locations(2,count-1) = residue_locations(2,count-1) + GRIDSIZE;
        case {'l','d','L','D'}
            residue_locations(1,count-1) = residue_locations(1,count-1) + GRIDSIZE;
        case {'j','a','J','A'}
            residue_locations(1,count-1) = residue_locations(1,count-1) - GRIDSIZE;
    end
    switch char(button)
        case {'i','j','k','l','I','J','K','L','a','s','d','w','A','S','D','W'}
            xpick = residue_locations(1,count-1);
            ypick = residue_locations(2,count-1);
            set(h(count-1),'Position',  [xpick - square_width/2, ypick-square_width/2,...
                square_width,square_width]);
    end


end
