function color_circles(imagex, residue_locations, which_res, what2plot, color_profile, square_width, file_name)
% COLOR_CIRCLES(imagex, residue_locations, which_res, what2plot, color_profile, square_width, file_name)
%
% Makes colored circles on secondary structure diagram. Reactivities from what2plot will be plotted
%  on positions specified by residue_locations in an order of which_res with box size
%  specified by square_width. Coloring scheme and range are according to color_profile.
% Saves to EPS file of given file name.
%
%
% Input
% =====
%   imagex             Required        Provides the [N x N x 3] RGB image read from file with
%                                       imread command.
%   residue_locations  Required        Provides the [2 x N] array of (X, Y) coordinates of
%                                       each letter positions on imagex. Can be selected
%                                       with pick_points.
%   which_res          Required        Provides the sequence positions of the input data.
%   what2plot          Required        Provides the input data.
%   colof_profile      Optional        Provides the coloring scheme and ranges. Default is
%                                       rainbow, offset 0, saturation values are average +/-
%                                       3 standard deviation.
%   square_width       Optional        Provides the size of boxes. Default is 24.
%   file_name          Optional        Provides file name for EPS file. Default is 'circles'.

if ~exist('square_width', 'var') || isempty(square_width); square_width = 24; end;
if ~exist('file_name','var') || isempty(file_name); file_name = 'circles'; end;

[color_scheme, d_offset, max_color, min_color] = parse_color_profile(color_profile);
[xsize, ysize, zsize] = size(imagex);

figure(); hold on; axis equal; axis off;
axis([0 ysize 0 xsize]);
set(gca, 'ydir', 'reverse')

count = 1; h = [];
for k = which_res
    x = residue_locations(1, k);
    y = residue_locations(2, k);
    colorplot = getcolor(what2plot(count) - d_offset, max_color, min_color, color_scheme);
    
    h(count) = rectangle('Position', [x-square_width/2, y-square_width/2, square_width, square_width], 'Curvature', [1 1]);
    set(h(count), 'edgecolor', colorplot, 'facecolor', colorplot)
    count = count + 1;
end;

set_print_page(gcf, 1);
print_save_figure(gcf, file_name);
