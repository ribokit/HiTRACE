function setcolormap(grayscaleon,maxprof)

if (grayscaleon)
    colorscale = 1-gray(maxprof);
else
    colorscale = jet(ceil(maxprof));
end

colormap(colorscale);
