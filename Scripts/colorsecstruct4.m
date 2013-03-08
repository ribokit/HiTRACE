function imagex_color = colorsecstruct4(imagex,offset,residue_locations, whichres,whattoplot,maxplot,minplot,colorscheme,makelegend,boxsize)
%
% imagex_color = colorsecstruct3(imagex,offset,residue_locations, whichres,whattoplot,maxplot,minplot,colorscheme,makelegend,boxsize)
%
% imagex            = RGB image [M1 x M2 x 3 matrix] read in from, say a tif file with 
%                      the 'imread' command.
% offset            = integer to add to 1, 2, ... N to get the actual positions on your image
% residue_locations = 2 x N matrix with the (x,y) positions of each letter 
%                      on the image. Can be selected with 'pickpoints'
%                      function.
% whichres          = sequence positions of the input data.
% whattoplot        = the input data (you may want to add or subtract a 
%                      constant so that the 'baseline' value is 0]. 
% maxplot           = [default 1.0] value at which colors should saturate in the positive
%                      direction
% minplot           = [default 1.0] value at which color should saturate in the negative direction 
% colorscheme       = [default 1] integer reflecting coloring -- type 'help getcolor'
%                      for list
% makelegend        = [default 1] 0 or 1 -- make a legend of the colors!
% boxsize           = [default 24] how big to make squares.
% 
% A useful command:
% hold off; image(imagex); hold on; axis equal; 
% 
% (C) R. Das, 2004-2012, Thomas Mann 2012.
[xsize,ysize,zsize]=size(imagex);
axis([0 ysize 0 xsize]); zoomedin = 0;

numres = length(residue_locations);

if (nargin<6) maxplot = max(whattoplot);end;
%if (nargin<7) maxplot2 = maxplot;end;
if (nargin<7) minplot = min(whattoplot);end;
%if (nargin<8) centerplot = (maxplot + minplot) / 2;end;
if ~exist('colorscheme') colorscheme = 1;end;
if ~exist( 'boxsize') boxsize = 24; end;
if size( whichres, 1) > 1; whichres = whichres'; end;

square_width = boxsize;
imagex_color = double(imagex);
count = 1;
for k = whichres
    x=residue_locations(1,k - offset);
    y=residue_locations(2,k - offset);
    %   h=rectangle('Position', ...
    %       [x - square_width/2, y - square_width/2, square_width,square_width]);
    colorplot = getcolor(whattoplot(count),maxplot,minplot,colorscheme); %converts the reactivity values in whattoplot into RGB values
    xbins = max(floor(y-square_width/2),1) :  min(floor(y+square_width/2),xsize); 
    ybins = max(floor(x-square_width/2),1) :  min(floor(x+square_width/2),ysize); 
    for n=1:3
      imagex_color(xbins,ybins,n) = double(imagex(xbins,ybins,n))*colorplot(n);
    end
    count=count+1;
end

%for squares where there is no data 
% for k=1:size(residue_locations,2)
%     if (isempty(find(k+offset == whichres)))
%         x=residue_locations(1,k);
%         y=residue_locations(2,k);
%         xbins = floor(y-square_width/2) :  min(floor(y+square_width/2),xsize); 
%         ybins = floor(x-square_width/2) :  min(floor(x+square_width/2),ysize); 
%         for n=1:3
%             imagex_color(xbins,ybins,n) = double(imagex(xbins,ybins,n))*0.7;
%         end
%     end
% end

%Make a "legend"
if ~exist('makelegend') makelegend=1;end;

if (makelegend)
    numlines = 10;
    x_offset = size( imagex,2) - 2*square_width;
    %x_offset = 5*square_width;
    sizebar = 4;
    bin_move_count = 1 % keeps track of where to put the individual color bars in the legend; there are always 21.
    for k=maxplot:-(0.05*(maxplot-minplot)):minplot %legend now goes from maxplot to minplot rather than 1 to -1; take 20 steps from maxplot to minplot
        %bin_move_count = -1;
        ybins = square_width:2*square_width;
        xbins = x_offset + ...
            [round(sizebar*-1*bin_move_count*square_width) : round(sizebar*-1*bin_move_count*square_width)+sizebar*square_width/numlines];
        colorplot = getcolor(k, maxplot,minplot,colorscheme);
        bin_move_count = bin_move_count - 0.1;
        for n=1:3
            imagex_color(xbins,ybins,n) = double(imagex(xbins,ybins,n))*colorplot(n);
        end
    end
end

hold off; image(imagex_color/256); hold on; axis equal

if (makelegend)
    k=-1;
    text(2*square_width,x_offset+round(sizebar*-1*k*square_width),['-',num2str(abs(minplot))]);
    k=+1;
    text(2*square_width,x_offset+round(sizebar*-1*k*square_width),['+',num2str(maxplot)]);
end

axis off
