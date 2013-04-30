function imagex_color = colorbases(imagex, offset, base_locations, residue_locations, which_res, what2plot, plot_limit, color_scheme, if_legend, pos_legend, square_width, labels, font_size, if_output)
%
% imagex_color = colorbases(imagex, offset, base_locations, residue_locations, which_res, what2plot, ...
%                           plot_limit, color_scheme, if_legend, pos_legend, square_width, ...
%                           labels, font_size, if_output)
%
% Colors base stubs on secondary structure diagram. Stubs orientation will be automatically 
%  determined by comparing base_locations and residue_locations. Reactivities from what2plot will be
%  plotted on positions specified by base_locations in an order of which_res with box size specified 
%  by square_width. Coloring scheme is specified by colorscheme and color ranges by plot_limit. Legend
%  will be added at position specified by pos_legend, with orientation specified by if_legend. Legend
%  labels text is specified by labels and size is specified by font_size. Colored image will be 
%  exported to hi-res non-compressed tiff file if asked.
% 
%
% Input
% =====
% imagex            = RGB image [M1 x M2 x 3 matrix] read in from, say a tif file with the 'imread' 
%                      command.
% offset            = integer to add to 1, 2, ... N to get the actual positions on your image.
% base_locations    = 2 x N matrix with the (x,y) positions of each 'base stub' position on the image.
%                      Can be selected with 'pickbases' or move_base_locations.
% residue_locations = 2 x N matrix with the (x,y) positions of each letter on the image. Can be 
%                      selected with 'pickpoints'.
% which_res         = sequence positions of the input data.
% what2plot         = the input data (you may want to add or subtract a constant so that the
%                      'baseline' value is 0]. 
% plot_limit        = [maxplot maxplot2]
%     maxplot       = [default 1.0] value at which colors should saturate in the positive direction.
%     maxplot2      = [default 1.0] value at which color should saturate in the negative direction. 
% color_scheme      = [default 1] integer reflecting coloring -- type 'help getcolor' for list.
% if_legend         = [makelegend orient]
%     makelegend    = [default 1] 0 or 1, for whether make a color legend.
%     orient        = [default 1] 0 or 1, for horizontal or vertical.
% pos_legend        = [corner x_offset y_offset size]
%     corner        = [default 0] 0 or 1 or 2 or 3, for top-left, top-right, bottom-left, and 
%                      bottom-right corners.
%     x_offset      = [default 2] value for distance of legend from vertical border, in unit of
%                      square_width.
%     y_offset      = [default 2] value for distance of legend from horizontal border, in unit of
%                      square_width.
%     size          = [default 4] value for length of legend, in unit of square_width.
% square_width      = [default 24] how big to make squares.
% labels            = {minlabel maxlabel}   
%     min_label     = text or number that corresponds to maxplot2.
%     max_label     = text or number that corresponds to maxplot.
% font_size         = [default 0.8] value for legend label font size, in unit of 
%                      square_width.
% if_output         = [default 1] 0 or 1, for whether output to hi-res tiff images.
%
%
% A useful command:
% hold off; image(imagex); hold on; axis equal; 
% 
% (C) R. Das, 2004-2012
% modified by T47, Apr 2013.
%

if nargin == 0; help(mfilename); return; end;

[x_image_size, y_image_size, z_image_size] = size(imagex);
if z_image_size ~= 3; fprintf('WARNING: invalid image input.\n'); end;
axis([0 y_image_size 0 x_image_size]); 
tiff_file_name = 'color_base_output.tiff';

%zoomedin = 0;
%numres = length(base_locations);

if size(which_res, 1) > 1; which_res = which_res'; end;
if ~exist('plot_limit','var') || length(plot_limit) < 2; 
    maxplot = max(abs(what2plot));
    maxplot2 = maxplot;
else
    maxplot = plot_limit(1);
    maxplot2 = plot_limit(2);
end;
if ~exist('color_scheme','var'); color_scheme = 1; end;

if ~exist('if_legend', 'var') || length(if_legend) < 2; if_legend = [1 1]; end;
makelegend = if_legend(1); legend_orient = if_legend(2);

if ~exist('pos_legend','var') || length(pos_legend) < 4; pos_legend = [0 2 2 4]; end;
legend_corner = pos_legend(1); legend_x_offset = pos_legend(2);
legend_y_offset = pos_legend(3); legend_size = pos_legend(4);

if ~exist('square_width','var'); square_width = 12; end;
if ~exist('labels', 'var') || length(labels) < 2; 
    max_label = ['+', num2str(abs(maxplot))];
    min_label = ['-', num2str(abs(maxplot2))];
else
    max_label = labels{2};
    min_label = labels{1};
end;
if ~exist('font_size','var'); font_size = 0.8; end;
if ~exist('if_output','var'); if_output = 1; end;

deviation = 5;
if length(square_width) == 1;
    square_extent = square_width / 2;
    square_width = square_width / 4;
else
    if length(square_width) > 2
        deviation = square_width(3);
    end;
    square_extent = square_width(1); 
    square_width = square_width(2);
end;


% output parameters to screen
fprintf('\n');
fprintf(['input image = ', num2str(x_image_size), ' x ', num2str(y_image_size), ' x ' ...
    num2str(z_image_size), ' matrix.\n']);
fprintf(['offset = ', num2str(offset), ', base_locations = ', num2str(size(base_locations, 1)), ' x ', num2str(size(base_locations, 2)), ...
    ' matrix, residue_locations = ', num2str(size(residue_locations, 1)), ' x ', num2str(size(residue_locations, 2)), ' matrix.\n']);
fprintf(['which_res = 1 x ', num2str(length(which_res)),' matrix, what2plot = 1 x ', num2str(length(what2plot)), ' matrix.\n']);
fprintf('\n');
fprintf(['color_scheme = ', num2str(color_scheme), ' , color saturation at ', num2str(maxplot), ' (labled ', max_label, ...
    ') and ', num2str(maxplot2), ' (labeled ', min_label, ').\n']);
fprintf(['square_width = ', num2str(square_width * 4), '.\n']);
label_orient = {'vertical', 'horizontal'};
fprintf(['make legend (', num2yn(makelegend), '), with ', label_orient{legend_orient + 1}, ' orientation.\n']);
label_pos = {'top-left', 'top-right', 'bottom-left', 'bottom-right'};
fprintf(['legend is located ', num2str(legend_x_offset), ' x ', num2str(legend_y_offset), ' square_width away from ', ...
    label_pos{legend_corner + 1}, ' corner,\n  with size 1 x ', num2str(legend_size), ' square_width, and label font size ', ...
    num2str(font_size), 'x square_width.\n']);
fprintf(['output to tiff file "', tiff_file_name, '" (', num2yn(if_output), ').\n']);
fprintf('\n');
fprintf('Press any key to continue...\n');
pause;

imagex_color = double(imagex);
count = 1;
for k = which_res
    x_pos = base_locations(1, k - offset);
    y_pos = base_locations(2, k - offset);
    
    x_dist = base_locations(1, k - offset) - residue_locations(1, k - offset);
    y_dist = base_locations(2, k - offset) - residue_locations(2, k - offset);
    %   h=rectangle('Position', ...
    %       [x - square_width/2, y - square_width/2, square_width,square_width]);
    colorplot = getcolor(what2plot(count), maxplot, maxplot2, color_scheme);
    
    % (t47)
    % determine tick orientation from the relative position of
    % base_locations and residue_locations
    if abs(x_dist) > abs(y_dist)
        x_thick = square_extent;
        y_thick = square_width;
        x_dev = deviation * sign(x_dist);         
        y_dev = 0;
    else
        x_thick = square_width;
        y_thick = square_extent;
        x_dev = 0;
        y_dev = deviation * sign(y_dist);
    end;    
    
    
    % (t47)
    % imagex is flipped... switch xbins and ybins
    bottom = max(floor(y_pos + y_dev - y_thick), 2); 
    top = min(floor(y_pos + y_dev + y_thick), x_image_size - 1);
    left = max(floor(x_pos + x_dev - x_thick), 2);
    right= min(floor(x_pos + x_dev + x_thick), y_image_size - 1); 
    
    xbins = bottom:top; 
    ybins = left:right;
    
    for n = 1:3
      imagex_color(xbins, ybins, n) = double(imagex(xbins,ybins,n)) * colorplot(n);
    end
    
%     MAKE_BOX_OUTLINE = 0;
%     if MAKE_BOX_OUTLINE;    
%       xbins = bottom:top; 
%       ybins = left-1;
%       for n = 1:3
%         imagex_color(xbins,ybins,n) = 0;
%       end;
%       
%       xbins = bottom:top; 
%       ybins = right+1;
%       for n = 1:3
%         imagex_color(xbins,ybins,n) = 0;
%       end;
%       
%       xbins =  bottom-1; 
%       ybins = left:right;
%       for n = 1:3
%         imagex_color(xbins,ybins,n) = 0;
%       end;
%       
%       xbins =  top+1; 
%       ybins = left:right;
%       for n = 1:3
%         imagex_color(xbins,ybins,n) = 0;
%       end;
%     end;
    
    count = count + 1;
end;

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


% make a "legend"
if makelegend
    numlines = legend_size * 5;
    square_width = square_width * 5;
    
    % (t47)
    % legend position determined by pos_legend
    switch legend_corner;
        case 0      % north-west
            x_offset = legend_x_offset * square_width;
            y_offset = legend_y_offset * square_width;
        case 1      % north-east
            x_offset = y_image_size - (legend_x_offset + legend_orient) * square_width - (1 - legend_orient) * legend_size * square_width;
            y_offset = legend_y_offset * square_width;
        case 2      % south-west
            x_offset = legend_x_offset * square_width;
            y_offset = x_image_size - (legend_y_offset + 1 - legend_orient) * square_width - legend_orient * legend_size * square_width;
        case 3      % south-east
            x_offset = y_image_size - (legend_x_offset + legend_orient) * square_width - (1 - legend_orient) * legend_size * square_width;
            y_offset = x_image_size - (legend_y_offset + 1 - legend_orient) * square_width - legend_orient * legend_size * square_width;
        otherwise
            fprintf('WARNING: invalid pos_legend corner parameter.\n');
            x_offset = legend_x_offset * square_width;
            y_offset = legend_y_offset * square_width;
    end;
    
    % color the legend 
    for k = 1:-(1/numlines):-1
        ybins_range = [1:square_width];
        xbins_range = legend_size * (1 - k) * square_width / 2 + [0:legend_size*square_width/numlines];
        
        % (t47)
        % bins area deteremined by legend_orient
        if legend_orient
            ybins = x_offset + ybins_range;
            xbins = y_offset + xbins_range;
        else
            ybins = x_offset + xbins_range;
            xbins = y_offset + ybins_range;
        end;
        
        % (t47)
        % found this necessary to avoid non-integer coordinates
        xbins = round(xbins); ybins = round(ybins);
        
        if k > 0;
            colorplot = getcolor( maxplot * k, maxplot, maxplot2, color_scheme);
        else
            colorplot = getcolor( maxplot2 * k, maxplot, maxplot2, color_scheme);
        end;
        for n = 1:3
            imagex_color(xbins, ybins, n) = double(imagex(xbins, ybins, n)) * colorplot(n);
        end;
    end;
end;


% (t47)
% embed text as part of image matrix
if legend_orient
    label_pos_1 = [x_offset + square_width, y_offset - square_width / 2];
    label_pos_2 = [x_offset + square_width, y_offset + (legend_size / 2 - 0.5) * square_width];
    label_pos_3 = [x_offset + square_width, y_offset + (legend_size - 0.5) * square_width];
else
    label_pos_1 = [x_offset - square_width, y_offset + square_width];
    label_pos_2 = [x_offset + (legend_size / 2 - 0.5) * square_width, y_offset + square_width];
    label_pos_3 = [x_offset + (legend_size - 0.5) * square_width, y_offset + square_width];    
end;

% (t47)
% add legend label
% text is rasterized into image
if makelegend
    ft_sz = round(square_width * font_size);
    if color_scheme == 8
        k = 0;
        H = vision.TextInserter('0');
        H.Color = [0 0 0]; H.Location = label_pos_2; H.FontSize = ft_sz;
        imagex_color = step(H, imagex_color);
    else
        H = vision.TextInserter(min_label);
        H.Color = [0 0 0]; H.Location = label_pos_3; H.FontSize = ft_sz;
        imagex_color = step(H, imagex_color);
    end;
    H = vision.TextInserter(max_label);
    H.Color = [0 0 0]; H.Location = label_pos_1; H.FontSize = ft_sz;
    imagex_color = step(H, imagex_color);
end;

hold off; image(imagex_color / 256); hold on; axis equal; axis off;

% (t47)
% output to hi-res non-compressed tiff
if if_output
    d_tiff = im2uint8(imagex_color / 256);
    imwrite(d_tiff, tiff_file_name, 'TIFF', ...
        'Resolution', 300, 'Compression', 'lzw');
end;



