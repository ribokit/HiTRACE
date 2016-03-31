function imagex_color = color_bases(imagex, base_locations, residue_locations, which_res, what2plot, color_profile, square_width)
%
% imagex_color = COLOR_BASES(imagex, base_locations, residue_locations, ...
%                            which_res, what2plot, color_profile, square_width)
%
% Colors base stubs on secondary structure diagram. Stubs orientation will be automatically 
%  determined by comparing base_locations and residue_locations. Reactivities from what2plot 
%  will be plotted on positions specified by base_locations in an order of which_res with 
%  box size specified by square_width. Coloring scheme and range are according to 
%  color_profile. 
% 
%
% Input
% =====
%   imagex             Required        Provides the [N x N x 3] RGB image read from file with
%                                       imread command.
%   base_locations     Required        Provides the [2 x N] array of (X, Y) coordinates of  
%                                       base stub positions on imagex. Can be selected with
%                                       pick_bases or move_base_locations.
%   residue_locations  Required        Provides the [2 x N] array of (X, Y) coordinates of 
%                                       each letter positions on imagex. Can be selected
%                                       with pick_points.
%   which_res          Required        Provides the sequence positions of the input data.
%   what2plot          Required        Provides the input data. 
%   colof_profile      Optional        Provides the coloring scheme and ranges. Default is 
%                                       rainbow, offset 0, saturation values are average +/-
%                                       3 standard deviation.
%   square_width       Optional        Provides the size of boxes. Default is 24.
%
% Output
% ======
%   imagex_color                       Gives the [N x N x 3] RGB image with coloring.
%
%
% (C) R. Das, 2004-2012
% modified by T47, Apr 2013 - May 2013.
%

if nargin == 0; help(mfilename); return; end;

[x_image_size, y_image_size, z_image_size] = size(imagex);
if z_image_size ~= 3; fprintf('WARNING: invalid image input.\n'); end;
axis([0 y_image_size 0 x_image_size]); 

if size(which_res, 1) > 1; which_res = which_res'; end;

if ~exist('color_profile','var'); color_profile = []; end;
color_profile_org = [17, 0, mean(what2plot) + 3 * std(what2plot), mean(what2plot) - 3 * std(what2plot)];
color_profile = is_valid_flag(color_profile, color_profile_org);
[color_scheme, d_offset, max_color, min_color] = parse_color_profile(color_profile);
if ~exist('square_width','var'); square_width = 12; end;

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
fprintf(['base_locations = ', num2str(size(base_locations, 1)), ' x ', num2str(size(base_locations, 2)), ...
    ' matrix, residue_locations = ', num2str(size(residue_locations, 1)), ' x ', num2str(size(residue_locations, 2)), ' matrix.\n']);
fprintf(['which_res = 1 x ', num2str(length(which_res)),' matrix, what2plot = 1 x ', num2str(length(what2plot)), ' matrix.\n']);
fprintf('\n');
fprintf(['color_scheme = ', num2str(color_scheme), '.\n']);
fprintf(['square_width = ', num2str(square_width), '.\n']);
fprintf('\n');

imagex_color = double(imagex);
count = 1;
for k = which_res
    x_pos = base_locations(1, k);
    y_pos = base_locations(2, k);
    
    x_dist = base_locations(1, k) - residue_locations(1, k);
    y_dist = base_locations(2, k) - residue_locations(2, k);
    colorplot = getcolor(what2plot(count) - d_offset, max_color, min_color, color_scheme);
    if what2plot(count) == -999; colorplot = [0.4 0.4 0.4]; end;
    
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

    count = count + 1;
end;

image_diagram(imagex_color);

