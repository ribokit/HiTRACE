function imagex_color = colorsecstruct(imagex,offset,residue_locations, whichres,whattoplot,maxplot,COLORCODE)
%
% imagex_color = colorsecstruct(imagex,offset,residue_locations, whichres,whattoplot,maxplot,COLORCODE)
%

if nargin == 0;  help( mfilename ); return; end;

if ~exist( 'COLORCODE' ) COLORCODE = 1; end;
%subplot(1,1,1); hold off; 
image(imagex); hold on; axis equal; 
[xsize,ysize,zsize]=size(imagex);
axis([0 ysize 0 xsize]); zoomedin = 0;

numres = length(residue_locations);

if (nargin<6) maxplot = max(abs(whattoplot));end;

square_width = 24;
imagex_color = double(imagex);
for k=whichres
    x=residue_locations(1,k - offset);
    y=residue_locations(2,k - offset);
    %   h=rectangle('Position', ...
    %       [x - square_width/2, y - square_width/2, square_width,square_width]);
    
    colorplot = getcolor(whattoplot(k),maxplot,COLORCODE);
    xbins = floor(y-square_width/2) :  min(floor(y+square_width/2),xsize); 
    ybins = floor(x-square_width/2) :  min(floor(x+square_width/2),ysize); 
    for n=1:3
        imagex_color(xbins,ybins,n) = double(imagex(xbins,ybins,n))*colorplot(n);
    end
    
end

%Make a "legend"
numlines = 10;  
sizebar = 4;
%x_offset =  size(imagex,1) - sizebar*2*square_width;
x_offset =  sizebar*2*square_width;

bounds = [1 -1];
if COLORCODE == 2; bounds  = [1 0]; end;

for k=bounds(1):-(1/numlines):bounds(2)
  ybins = square_width:2*square_width;
  xbins = x_offset + ...
	  [round(sizebar*-1*k*square_width) : round(sizebar*-1*k*square_width)+sizebar*square_width/numlines];
  colorplot = getcolor(k, 1.0,COLORCODE);
  for n=1:3
    imagex_color(xbins,ybins,n) = double(imagex(xbins,ybins,n))*colorplot(n);
  end
end

maxplotstring = ['+',num2str(maxplot)];
minplotstring = ['-',num2str(maxplot)];
if COLORCODE == 2; 
  maxplotstring = num2str(maxplot);
  minplotstring = '0'; 
end;

%subplot(1,1,1); 
hold off; image(imagex_color/256); hold on; axis equal
k=bounds(2);
text(2.2*square_width,x_offset+round(sizebar*-1*k*square_width),minplotstring);
k=bounds(1);
text(2.2*square_width,x_offset+round(sizebar*-1*k*square_width),maxplotstring);


rectangle( 'Position', [square_width, x_offset+round(sizebar*-1*bounds(1)*square_width),  ...
		    square_width, sizebar*abs(bounds(2)-bounds(1))*square_width], 'edgecolor','k' );

axis off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  colorplot = getcolor(colorvalue, maxplot,COLORCODE);
if isnan( colorvalue ); colorplot = [0.7 0.7 0.7]; return;end;

if COLORCODE == 1
  if (colorvalue>0)
    colorplot = [1, max(1-colorvalue/maxplot,0), max(1-colorvalue/maxplot,0)] ;
  else  
    colorplot = [max(1+colorvalue/maxplot,0), 1, max(1+colorvalue/maxplot,0) ] ;
  end
else
  %cmap = interp1( [1 10 80 100], [1 1 1; 1 1 1; 1 0.5333 0; 1 0 0 ], [1:100] );
  scalefactor = max(min(colorvalue/maxplot,1.0), 0.0 ) * 100;  
  colorplot = interp1( [0 10 80 100], [1 1 1; 1 1 1; 1 0.5333 0; 1 0 0 ], scalefactor );
end

