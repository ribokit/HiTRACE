function imagex_color = color_legend(imagex, color_profile, square_width, orientation, position, labels, font_size)

%
% imagex_color = COLOR_LEGEND(imagex, color_profile, square_width...
%                             orientation, position, labels, font_size)
%
% Colors legend bar at position specified by position, with orientation specified by 
%  orientation. Color scheme and range are according to color_profile. Legend labels text is 
%  specified by labels and size is specified by font_size, all in units of square_width.
% 
%
% Input
% =====
%   imagex             Required        Provides the [N x N x 3] RGB image read from file with
%                                       imread command.
%   colof_profile      Optional        Provides the coloring scheme and ranges. Default is 
%                                       rainbow, offset 0, saturation values are average +/-
%                                       3 standard deviation.
%   square_width       Optional        Provides the size of boxes. Default is 24. Serves as
%                                       unit for position and font_size.
%   orientation        Optional        Provides the orientation of legend bar. Default is 
%   [orient direct]                     [1 1].
%                                      'ORIENT' specifies the orientation (0 for horizontal;
%                                        1 for vertical).
%                                      'DIRECT' specifies color gradient direction (0 for 
%                                        light to dark; 1 for dark to light).
%   position           Optional        Provides the position of legend bar on the image.
%   [corner, x_offset,                  Default is [0 2 2 4].
%    y_offset, size]                   'CORNER' specifies on which corner will the legend be
%                                        added (0 for top-left; 1 for top-right; 2 for 
%                                        bottom-left; 3 for bottom-right).
%                                      'X_OFFSET' specifies the distance of legend from the
%                                        vertical border, in units of square_width.
%                                      'Y_OFFSET' specifies the distance of legend from the
%                                        horizontal border, in units of square_width.
%                                      'SIZE' specifies the length of legend bar, in units of
%                                        square_width.
%   labels             Optional        Provides labels for color saturation value. Default is
%   {min_label,                         {'color_min', 'color_max'} from color_profile.
%    max_label}                        'MIN_LABEL' specifies the low-end color saturation value.
%                                      'MAX_LABEL' specifies the high-end color saturation 
%                                        value.
%   font_size          Optional        Provides the legend label font size, in units of 
%                                       square_width. Default is 0.8.
%
% Output
% ======
%   imagex_color                       Gives the [N x N x 3] RGB image with legend.
%
%
% (C) R. Das, 2004-2012
% modified by T47, Apr 2013 - May 2013.
%

if nargin == 0; help(mfilename); return; end;

[x_image_size, y_image_size, z_image_size] = size(imagex);
if z_image_size ~= 3; fprintf('WARNING: invalid image input.\n'); end;
axis([0 y_image_size 0 x_image_size]); 

[color_scheme, color_center, max_color, min_color] = parse_color_profile(color_profile);
min_color = abs(min_color);
if ~exist('square_width','var'); square_width = 12; end;

orientation_org = [1 1];
if ~exist('orientation','var') || isempty(orientation) ; 
    orientation = orientation_org;
else
    if length(orientation) < length(orientation_org);
        orientation_org(1:length(orientation)) = orientation;
        orientation = orientation_org;
    end;
end;
legend_orient = orientation(1); legend_direct = orientation(2);
if legend_direct == 0; legend_direct = -1; end;

position_org = [0 2 2 4];
if ~exist('position','var') || isempty(position) ; 
    position = position_org;
else
    if length(position) < length(position_org);
        position_org(1:length(position)) = position;
        position = position_org;
    end;
end;
legend_corner = position(1); legend_x_offset = position(2);
legend_y_offset = position(3); legend_size = position(4);

if ~exist('labels', 'var') || length(labels) < 2; 
    max_label = ['+', num2str(max_color + color_center)];
    min_label = ['-', num2str(- min_color + color_center)];
else
    max_label = labels{2};
    min_label = labels{1};
end;
if ~exist('font_size','var'); font_size = 0.8; end;

label_orient = {'vertical', 'horizontal'};
color_orient = {'dark (hot)', '', 'light (cool)'};
fprintf(['Legend with ', label_orient{legend_orient + 1}, ' orientation, color from ', ...
    color_orient{2-legend_direct}, ' to ', color_orient{2+legend_direct}, '.\n']);
label_pos = {'top-left', 'top-right', 'bottom-left', 'bottom-right'};
fprintf(['Legend is located ', num2str(legend_x_offset), ' x ', num2str(legend_y_offset), ' square_width away from ', ...
    label_pos{legend_corner + 1}, ' corner,\n  with size 1 x ', num2str(legend_size), ' square_width, and label font size ', ...
    num2str(font_size), 'x square_width.\n']);


imagex_color = double(imagex);
numlines = legend_size * 5;

% legend position determined by position
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
        fprintf('WARNING: invalid position corner parameter.\n');
        x_offset = legend_x_offset * square_width;
        y_offset = legend_y_offset * square_width;
end;

% color the legend
for k = 1:-(1/numlines):-1
    ybins_range = 1:square_width;
    xbins_range = legend_size * (1 - k) * square_width / 2 + [0:legend_size*square_width/numlines];
    
    % bins area deteremined by legend_orient
    if legend_orient
        ybins = x_offset + ybins_range;
        xbins = y_offset + xbins_range;
    else
        ybins = x_offset + xbins_range;
        xbins = y_offset + ybins_range;
    end;
    
    % found this necessary to avoid non-integer coordinates
    xbins = round(xbins); ybins = round(ybins);
    
    % color direction option
    if k * legend_direct > 0;
        colorplot = getcolor( max_color * k * legend_direct, max_color, min_color, color_scheme);
    else
        colorplot = getcolor( min_color * k * legend_direct, max_color, min_color, color_scheme);
    end;
    for n = 1:3
        imagex_color(xbins, ybins, n) = double(imagex(xbins, ybins, n)) * colorplot(n);
    end;
end;

% label position adjusted to bar orientation
if legend_orient
    label_pos_1 = [x_offset + square_width, y_offset - square_width / 2];
    label_pos_2 = [x_offset + square_width, y_offset + (legend_size / 2 - 0.5) * square_width];
    label_pos_3 = [x_offset + square_width, y_offset + (legend_size - 0.5) * square_width];
else
    label_pos_1 = [x_offset - square_width, y_offset + square_width];
    label_pos_2 = [x_offset + (legend_size / 2 - 0.5) * square_width, y_offset + square_width];
    label_pos_3 = [x_offset + (legend_size - 0.5) * square_width, y_offset + square_width];    
end;

% add legend label
% text is rasterized into image
ft_sz = round(square_width * font_size);
if color_scheme == 8
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

hold off; image(imagex_color / 256); hold on; axis equal; axis off;

