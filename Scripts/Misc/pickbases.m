function base_locations = pickbases(imagex,offset,residue_locations, base_locations);
figure(1); subplot(1,1,1); hold off; image(imagex); hold on
[xsize,ysize,zsize]=size(imagex);
axis([0 ysize 0 xsize]); zoomedin = 0;

square_width = 40;
stop_pick = 0;

count = length(residue_locations)+1;
for k=1:count-1
    xpick = residue_locations(1,k);
    ypick = residue_locations(2,k);
    hres = rectangle('Position',...
        [xpick - square_width/2, ypick-square_width/2,...
            square_width,square_width]);
    set(hres,'edgecolor','b');
    
    %    hnum = text(xpick,ypick,num2str(k+offset));
    
end

count=1;

if (nargin>3)
    count = length(base_locations)+1;
    for k=1:count-1
        xpick = base_locations(1,k);
        ypick = base_locations(2,k);
        h(k) = rectangle('Position',...
            [xpick - square_width/2, ypick-square_width/2,...
                square_width,square_width]);
        set(hres,'edgecolor','r');
        
        %    hnum = text(xpick,ypick,num2str(k+offset));
        
    end
end


while (stop_pick<1)
    title(['next: ', num2str(count+offset)])    
    [xpick,ypick,button] = ginput(1);
    switch button
        case 1
            distx = xpick-residue_locations(1,count) ;
            disty = ypick-residue_locations(2,count) ;
            dist = sqrt(distx^2+disty^2); 
            desired_dist = square_width;
            distx = distx*desired_dist/dist;
            disty = disty*desired_dist/dist;
            coarseness = desired_dist/2;
            distx = coarseness*round(distx / coarseness);    
            disty = coarseness*round(disty / coarseness);    
            xpick = residue_locations(1,count)  + distx;
            ypick = residue_locations(2,count)  + disty;
            base_locations(:,count) = [xpick;ypick];
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
    
