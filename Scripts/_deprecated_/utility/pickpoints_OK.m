function residue_locations = pickpoints_OK(imagex,offset,residue_locations);
%
% residue_locations = pickpoints(imagex,offset,residue_locations);
%
if nargin == 0;  help( mfilename ); return; end;

figure(1); subplot(1,1,1); hold off; image(imagex); hold on
[xsize,ysize,zsize]=size(imagex);
axis([0 ysize 0 xsize]); zoomedin = 0;

square_width = 24;
stop_pick = 0;
count=1;

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
            xpick = xpick - mod(xpick,4);
            ypick = ypick - mod(xpick,4);
            residue_locations(:,count) = [xpick;ypick];
            h(count) = rectangle('Position',...
                [xpick - square_width/2, ypick-square_width/2,...
                    square_width,square_width]);
            set(h(count),'edgecolor','r');
            count= count+1;
            title(['next: ', num2str(count+offset)])    
        case 2
            count= count-1;
            set(h(count),'visible','off');
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
            case 'q'
                stop_pick=1;
        end
    end
    
