function residue_locations = pickpoints(imagex,offset,residue_locations,square_width)

figure(1); subplot(1,1,1); hold off; image(imagex); hold on
[xsize,ysize,zsize]=size(imagex);
axis([0 ysize 0 xsize]); zoomedin = 0;

if ~exist('square_width')
  square_width = 40;
end

stop_pick = 0;
count=1;
GRIDSIZE=1;

if (nargin>2)
    count = length(residue_locations)+1;
    for k=1:count-1
        xpick = residue_locations(1,k);
        ypick = residue_locations(2,k);
        h(count) = rectangle('Position',...
            [xpick - square_width/2, ypick-square_width/2,...
                square_width,square_width]);
        set(h(count),'edgecolor','b');
    end
end  

while (stop_pick<1)
    title(['next: ', num2str(count+offset)])
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
	      size( h )
	      erasesquare
	      get(h(erasesquare))
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
            residue_locations(2,count-1) = residue_locations(2,count-1) - 1;
        case {'k','s','K','S'}
            residue_locations(2,count-1) = residue_locations(2,count-1) + 1;
        case {'l','d','L','D'}
            residue_locations(1,count-1) = residue_locations(1,count-1) + 1;
        case {'j','a','J','A'}
            residue_locations(1,count-1) = residue_locations(1,count-1) - 1;
    end
    switch char(button)
        case {'i','j','k','l','I','J','K','L','a','s','d','w','A','S','D','W'}
            xpick = residue_locations(1,count-1);
            ypick = residue_locations(2,count-1);
            set(h(count-1),'Position',  [xpick - square_width/2, ypick-square_width/2,...
                square_width,square_width]);
    end


end
