function  colorplot = getcolor(colorvalue, maxplot,maxplot2,colorscheme);

if ~exist('colorscheme') colorscheme = 1; end;

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
        if (frac_color<0) frac_color = 0; end;
        if (frac_color>1) frac_color = 1; end;
        colorplot = frac_color*colorval2 + (1-frac_color)*colorval1;
    case 3 % Blue-gray-yellow
        colorval1 = [0 0 1];
        colorval2 = [1 1 0];
        maxplot2 = -abs(maxplot2);
        frac_color = (colorvalue - maxplot2)/(maxplot-maxplot2);
        if (frac_color<0) frac_color = 0; end;
        if (frac_color>1) frac_color = 1; end;
        colorplot = frac_color*colorval2 + (1-frac_color)*colorval1;

    case 4 % Magenta-lightblue-cyan
        colorval1 = [1 0 1];
        colorval2 = [0 1 1];
        maxplot2 = -abs(maxplot2);
        frac_color = (colorvalue - maxplot2)/(maxplot-maxplot2);
        if (frac_color<0) frac_color = 0; end;
        if (frac_color>1) frac_color = 1; end;
        colorplot = frac_color*colorval2 + (1-frac_color)*colorval1;
    case 5 % Red-gray-cyan
        colorval1 = [1 0 0];
        colorval2 = [0 1 1];
        maxplot2 = -abs(maxplot2);
        frac_color = (colorvalue - maxplot2)/(maxplot-maxplot2);
        if (frac_color<0) frac_color = 0; end;
        if (frac_color>1) frac_color = 1; end;
        colorplot = frac_color*colorval2 + (1-frac_color)*colorval1;
    case 6 % Red-gray-cyan
        colorval1 = [1 0 0];
        colorval2 = [0 0 1];
        maxplot2 = -abs(maxplot2);
        frac_color = (colorvalue - maxplot2)/(maxplot-maxplot2);
        if (frac_color<0) frac_color = 0; end;
        if (frac_color>1) frac_color = 1; end;
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
end

