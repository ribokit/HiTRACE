function imagex_color = colorbases(imagex,offset,base_locations, residue_locations, whichres,whattoplot,maxplot,maxplot2,colorscheme,makelegend,square_width,maxlabel,minlabel)
%
% imagex_color = colorbases(imagex,offset,base_locations, residue_locations, whichres,whattoplot,maxplot,maxplot2,colorscheme,makelegend,square_width,maxlabel,minlabel)
%
% imagex            = RGB image [M1 x M2 x 3 matrix] read in from, say a tif file with 
%                      the 'imread' command.
% offset            = integer to add to 1, 2, ... N to get the actual positions on your image
% residue_locations = 2 x N matrix with the (x,y) positions of each 'base stub' position 
%                      on the image. Can be selected with 'pickbases'
% residue_locations = 2 x N matrix with the (x,y) positions of each letter 
%                      on the image. Can be selected with 'pickpoints'
% whichres          = sequence positions of the input data.
% whattoplot        = the input data (you may want to add or subtract a 
%                      constant so that the 'baseline' value is 0]. 
% maxplot           = [default 1.0] value at which colors should saturate in the positive
%                      direction
% maxplot2          = [default 1.0] value at which color should saturate in the negative direction 
% colorscheme       = [default 1] integer reflecting coloring -- type 'help getcolor'
%                      for list
% makelegend        = [default 1] 0 or 1 -- make a legend of the colors!
% squarewidth       = [default 24] how big to make squares.
% minlabel          = text or number that corresponds to maxplot2
% maxlabel          = text or number that corresponds to maxplot
%
% A useful command:
% hold off; image(imagex); hold on; axis equal; 
% 
% (C) R. Das, 2004-2012

[xsize,ysize,zsize]=size(imagex);
axis([0 ysize 0 xsize]); zoomedin = 0;

numres = length(base_locations);

if (nargin<7) maxplot = max(abs(whattoplot));end;
if (nargin<8) maxplot2 = maxplot;end;
if ~exist('colorscheme') colorscheme = 1;end;
    
if ~exist( 'square_width'); 
        square_width = 12;
end;
deviation = 5;
if length( square_width ) == 1
    square_extent = square_width/2;
    square_width = square_width/4;
else
    if length( square_width ) > 2
        deviation = square_width(3);
    end
    square_extent = square_width(1); 
    square_width = square_width(2);
end



imagex_color = double(imagex);
count = 1;
for k=whichres
    x=base_locations(1,k - offset);
    y=base_locations(2,k - offset);
    
    distx = base_locations(1,k - offset)-residue_locations(1,k-offset);
    disty = base_locations(2,k - offset)-residue_locations(2,k-offset);
    %   h=rectangle('Position', ...
    %       [x - square_width/2, y - square_width/2, square_width,square_width]);
    colorplot = getcolor(whattoplot(count),maxplot,maxplot2,colorscheme);
        
    if abs(distx)>abs(disty)
        xthickness = square_extent;
        ythickness = square_width;
        xdeviation = deviation*sign(distx);         ydeviation = 0;       
    else
        xthickness = square_width;
        ythickness = square_extent;
         xdeviation = 0; ydeviation = deviation*sign(disty);   
    end    
    
    bottom = max(floor(y+ydeviation-ythickness),2); top=min(floor(y+ydeviation+ythickness),xsize-1);
    left = max(floor(x+xdeviation-xthickness),2) ;right= min(floor(x+xdeviation+xthickness),ysize-1); 
    xbins =  bottom:top; 
    ybins = left:right;
    for n=1:3
      imagex_color(xbins,ybins,n) = double(imagex(xbins,ybins,n))*colorplot(n);
        %imagex_color(xbins,ybins,n) = 255*colorplot(n);
    end
    
    MAKE_BOX_OUTLINE = 0;
    if MAKE_BOX_OUTLINE;    
      xbins =  bottom:top; 
      ybins = left-1;
      for n=1:3
        imagex_color(xbins,ybins,n) = 0;
      end
      
      xbins =  bottom:top; 
      ybins = right+1;
      for n=1:3
        imagex_color(xbins,ybins,n) = 0;
      end
      
      xbins =  bottom-1; 
      ybins = left:right;
      for n=1:3
        imagex_color(xbins,ybins,n) = 0;
      end
      
      xbins =  top+1; 
      ybins = left:right;
      for n=1:3
        imagex_color(xbins,ybins,n) = 0;
      end
    end
    
    count=count+1;
end

%for squares where there is no data 
%for k=1:size(base_locations,2)
%    if (isempty(find(k+offset == whichres)))
%    x=base_locations(1,k);
%    y=base_locations(2,k);
%    xbins = floor(y-square_width/2) :  min(floor(y+square_width/2),xsize); 
%    ybins = floor(x-square_width/2) :  min(floor(x+square_width/2),ysize); 
%        for n=1:3
%            imagex_color(xbins,ybins,n) = double(imagex(xbins,ybins,n))*0.7;
%        end
%    end
%end

%Make a "legend"
if ~exist('makelegend') 
  makelegend=1; 
end

if (makelegend)
numlines = 10;  
x_offset = size( imagex,2) - 2*square_width;
y_offset = 5*square_width;
sizebar = 4;
for k=1:-(1/numlines):-1
    ybins = x_offset -square_width + [1:square_width];
    xbins = y_offset + ...
        [round(sizebar*-1*k*square_width) : round(sizebar*-1*k*square_width)+sizebar*square_width/numlines];
    if k > 0
      colorplot = getcolor( maxplot*k, maxplot,maxplot2,colorscheme);
    else
      colorplot = getcolor( maxplot2*k, maxplot,maxplot2,colorscheme);
    end
    for n=1:3
        imagex_color(xbins,ybins,n) = double(imagex(xbins,ybins,n))*colorplot(n);
    end
end
end
if ~exist( 'maxlabel' ) maxlabel = ['+',num2str(abs(maxplot))]; end;

if ~exist( 'minlabel' ) maxlabel = ['-',num2str(abs(maxplot2))]; end;

hold off; image(imagex_color/256); hold on; axis equal
if (makelegend)
  if colorscheme == 8
    k=0;
    text(x_offset,y_offset+round(sizebar*-1*k*square_width),'0');
  else
    k=-1;
    text(x_offset,y_offset+round(sizebar*-1*k*square_width),minlabel);
  end
  k=+1;
  text(x_offset,y_offset+round(sizebar*-1*k*square_width),maxlabel);
end

axis off


