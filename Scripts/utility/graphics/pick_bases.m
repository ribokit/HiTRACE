function base_locations = pickbases(imagex,offset,residue_locations, base_locations,square_width);
% base_locations = pickbases(imagex,offset,residue_locations, base_locations,square_width);
%
% imagex            = RGB image [M1 x M2 x 3 matrix] read in from, say a tif file with 
%                      the 'imread' command.
% offset            = integer to add to 1, 2, ... N to get the actual positions on your image
% residue_locations = 2 x N matrix with the (x,y) positions of each 'sequence' position 
%                      on the image. Use 'pickpoints' to get this.
% base_locations    = 2 x N matrix with the (x,y) positions of each 'base stub' position 
%                      on the image. The script will try to align these vertically or
%                      horizontally with residue_locations. Initially you can set this to the empty set [].
% squarewidth       = [default 24] how big to make squares.

if nargin == 0;  help( mfilename ); return; end;

figure(1); subplot(1,1,1); hold off; image(imagex); hold on
[xsize,ysize,zsize]=size(imagex);
axis([0 ysize 0 xsize]); zoomedin = 0;
figure_full_screen(); axis equal;

if ~exist( 'square_width','var' ); square_width = 40; end;
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
    count = size(base_locations,2)+1;
    for k=1:count-1
        xpick = base_locations(1,k);
        ypick = base_locations(2,k);
        h(k) = rectangle('Position',...
            [xpick - square_width/2, ypick-square_width/2,...
                square_width,square_width]);
        set(h(k),'edgecolor','r');        
        %    hnum = text(xpick,ypick,num2str(k+offset));        
    end
end

set(gcf,'color','w');
axis off

while (stop_pick<1)
  title(['next: ', num2str(count+offset)])    
  [xpick,ypick,button] = ginput(1);
  switch button
   case 1
    if ( count <= size( residue_locations, 2 ) ) 
        [xpick,ypick] = grid_fix( xpick, ypick, residue_locations, count, square_width );
        base_locations(:,count) = [xpick;ypick];
        h(count) = rectangle('Position',...
            [xpick - square_width/2, ypick-square_width/2,...
            square_width,square_width]);
        set(h(count),'edgecolor','r');
        count= count+1;
        title(['next: ', num2str(count+offset)])
    else
        fprintf( 'Cannot pick another point ... number of base locations cannot exceed number of residue locations\n' );
    end
   case 2
    [dummy, erasesquare] = min( (base_locations(1,:) - xpick).^2 + (base_locations(2,:) - ypick).^2 );
    if ( erasesquare > 0 )
      set(h(erasesquare),'visible','off');
      title(['Replace Pick ', num2str(erasesquare+offset), ' Now'])
      [xpick,ypick,button] = ginput(1);
      [xpick,ypick] = grid_fix( xpick, ypick, residue_locations, erasesquare, square_width );
      base_locations(:,erasesquare) = [xpick;ypick];
      set( h(erasesquare), 'Position', ...
			[xpick - square_width/2, ypick-square_width/2,...
			 square_width,square_width]);      
      set(h(erasesquare),'visible','on');
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
   case 'q'
    stop_pick=1;
  end
end
    

    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [xpick,ypick] = grid_fix( xpick, ypick, residue_locations, count, square_width );

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
