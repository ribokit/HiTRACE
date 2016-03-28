function  colorplot = getcolor(colorvalue, maxplot,maxplot2,colorscheme)
%
% colorplot = getcolor(colorvalue, maxplot,maxplot2,colorscheme);
%
% colorvalue = some numerical value that you want to convert to an (R,G,B)
%                color
% maxplot    = maximum value in positive direction.
% maxplot2   = maxmium value in negative direction.
% colorscheme = one of the following integers:
%
%            ( 1) red-white-blue
%            ( 2) cyan-green-yellow
%            ( 3) blue-gray-yellow
%            (-3) blue-gray-yellow [faded]
%            ( 4) Magenta-lightblue-cyan
%            ( 5) Red-gray-cyan
%            ( 6) Red-gray-cyan
%            ( 7) Red-green-blue
%            (17) Red-green-blue [FADED]
%            ( 9) Red-green [faded]
%            ( 8) Red-orange-white
%            (10) Green-white-red
%            (11) blue-white-red
%            (12) red-yellow-white [t47]
%            (13) red-white-blue [t47]
%
% (C) R. Das, 2008-2012
%

if nargin == 0;  help( mfilename ); return; end;

if ~exist('colorscheme','var'); colorscheme = 1; end;

if isnan( colorvalue); colorplot = [0.7 0.7 0.7]; return;end;

switch colorscheme
    case 1 %Red-white-blue
        if (colorvalue>0)
            colorplot = [1, max(1-colorvalue/maxplot,0), max(1-colorvalue/maxplot,0)] ;
        else
            colorplot = [max(1+colorvalue/abs(maxplot2),0),  max(1+colorvalue/abs(maxplot2),0),1 ] ;
        end
        return
    case 2 % Cyan-green-yellow
        colorval1 = [0 1 1];
        colorval2 = [1 1 0];
        maxplot2 = -abs(maxplot2);
        frac_color = (colorvalue - maxplot2)/(maxplot-maxplot2);
        if (frac_color<0); frac_color = 0; end;
        if (frac_color>1); frac_color = 1; end;
        colorplot = frac_color*colorval2 + (1-frac_color)*colorval1;
    case 3 % Blue-gray-yellow
        colorval1 = [0 0 1];
        colorval2 = [1 1 0];
        maxplot2 = -abs(maxplot2);
        frac_color = (colorvalue - maxplot2)/(maxplot-maxplot2);
        if (frac_color<0); frac_color = 0; end;
        if (frac_color>1); frac_color = 1; end;
        colorplot = frac_color*colorval2 + (1-frac_color)*colorval1;
        
    case -3  % Blue-gray-yellow, faded
        colorval1 = [0.5 0.5 1];
        colorval2 = [1 1 0];
        maxplot2 = -abs(maxplot2);
        frac_color = (colorvalue - maxplot2)/(maxplot-maxplot2);
        if (frac_color<0); frac_color = 0; end;
        if (frac_color>1); frac_color = 1; end;
        colorplot = frac_color*colorval2 + (1-frac_color)*colorval1;
        
    case 4 % Magenta-lightblue-cyan
        colorval1 = [1 0 1];
        colorval2 = [0 1 1];
        maxplot2 = -abs(maxplot2);
        frac_color = (colorvalue - maxplot2)/(maxplot-maxplot2);
        if (frac_color<0); frac_color = 0; end;
        if (frac_color>1); frac_color = 1; end;
        colorplot = frac_color*colorval2 + (1-frac_color)*colorval1;
    case 5 % Red-gray-cyan
        colorval1 = [1 0 0];
        colorval2 = [0 1 1];
        maxplot2 = -abs(maxplot2);
        frac_color = (colorvalue - maxplot2)/(maxplot-maxplot2);
        if (frac_color<0); frac_color = 0; end;
        if (frac_color>1); frac_color = 1; end;
        colorplot = frac_color*colorval2 + (1-frac_color)*colorval1;
    case 6 % Red-gray-cyan
        colorval1 = [1 0 0];
        colorval2 = [0 0 1];
        maxplot2 = -abs(maxplot2);
        frac_color = (colorvalue - maxplot2)/(maxplot-maxplot2);
        if (frac_color<0); frac_color = 0; end;
        if (frac_color>1); frac_color = 1; end;
        colorplot = frac_color*colorval2 + (1-frac_color)*colorval1;
    case 7 % Red-green-blue
        if (colorvalue>0)
            frac_color = min(colorvalue/maxplot,1);
            if (frac_color < 0.5)
                colorplot = [2*frac_color,1,0] ;
            else
                colorplot = [1,2*(1-frac_color),0] ;
            end;
        else
            frac_color = min(abs(colorvalue/maxplot2),1);
            if (frac_color < 0.5)
                colorplot = [0,1,2*frac_color] ;
            else
                colorplot = [0,2*(1-frac_color),1] ;
            end;
        end
        GREENSCALE=0.8;
        colorplot(2) = colorplot(2)*GREENSCALE;
    case 17 % Red-green-blue [FADED]
        if (colorvalue>0)
            frac_color = min(colorvalue/maxplot,1);
            if (frac_color < 0.5)
                colorplot = [2*frac_color,1,0] ;
            else
                colorplot = [1,2*(1-frac_color),0] ;
            end;
        else
            frac_color = min(abs(colorvalue/maxplot2),1);
            if (frac_color < 0.5)
                colorplot = [0,1,2*frac_color] ;
            else
                colorplot = [0,2*(1-frac_color),1] ;
            end;
        end
        GREENSCALE=0.8;
        colorplot(2) = colorplot(2)*GREENSCALE;
        colorplot = 1.0 - 0.7*(1.0 - colorplot);
    case 9 % Red-green [faded]
        colorval2 = [1 0.3 0.3];
        colorval1 = [0.3 0.9 0.5];
        maxplot2 = -abs(maxplot2);
        frac_color = (colorvalue - maxplot2)/(maxplot-maxplot2);
        if (frac_color<0); frac_color = 0; end;
        if (frac_color>1); frac_color = 1; end;
        colorplot = frac_color*colorval2 + (1-frac_color)*colorval1;
    case 8 %Red-orange-white
        a = (colorvalue - maxplot2)/(maxplot-maxplot2);
        
        if (a >= 0.5);
            colorplot = [1, min(max(1-a,0),1) ,0];
        else
            colorplot = [1, min(max(1-a,0),1) , min(max((0.5-a)*2,0),1)];
        end;
        return
    case 12 %Red-yellow-white [t47]
        a = (colorvalue - maxplot2)/(maxplot-maxplot2);
        a = min(max(a,0),1);
        
        if (a >= 0.33);
            colorplot = [1, 0.76-0.5*2*(a-0.33),0.21*2*(a-0.33)];
        else
            colorplot = [1, 1-0.24*3.33*a, 1-a*3.33];
        end;
        colorplot(colorplot < 0) = 0;
        colorplot(colorplot > 1) = 1;
        return        
    case 13 %Red-white-blue [t47]
        a = (colorvalue - maxplot2)/(maxplot-maxplot2);
        a = min(max(a,-1),1);
        
        if (a >= 0.66);
            colorplot = [1, 0.76-0.5*3.33*(a-0.66),0.21*3.33*(a-0.66)];
        elseif (a >= 0.5);
            colorplot = [1, 1-0.24*6*(a-0.5), 1-6*(a-0.5)];
        elseif (a >= 0.16);
            colorplot = [(a-0.16)*3, 0.66+(a-0.16)*3*0.34, 1];
        else
            colorplot = [0.33*(1-a*6), 0.58+a*6*0.08, 0.84+a*6*0.16];
        end;
        colorplot(colorplot < 0) = 0;
        colorplot(colorplot > 1) = 1;
        return        

    
    case 10 %Green-white-red
        colorplot = [1,1,1];
        a = colorvalue/maxplot;
        if (colorvalue>0)
            colorplot = [1, max(1-colorvalue/maxplot,0), max(1-colorvalue/maxplot,0)] ;
        else
            colorplot = [max(1+colorvalue/abs(maxplot2),0),  1.0, max(1+colorvalue/abs(maxplot2),0) ] ;
        end
        return
    case 11 %blue-white-red
        colorplot = [1,1,1];
        a = colorvalue/maxplot;
        if (colorvalue>0)
            colorplot = [1, max(1-colorvalue/maxplot,0), max(1-colorvalue/maxplot,0)] ;
        else
            colorplot = [max(1+colorvalue/abs(maxplot2),0),  max(1+colorvalue/abs(maxplot2),0), 1 ] ;
        end
%         colorplot = 1 - 0.5*( 1-colorplot);
        return
end

