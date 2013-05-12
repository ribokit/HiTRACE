function imagex_color = color_residues(imagex, residue_locations, which_res, what2plot, color_profile, square_width)
%
% imagex_color = COLOR_RESIDUES(imagex, residue_locations, ...
%                               which_res, what2plot, color_profile, square_width)
%
% Colors letters on secondary structure diagram. Reactivities from what2plot will be plotted 
%  on positions specified by residue_locations in an order of which_res with box size 
%  specified by square_width. Coloring scheme and range are according to color_profile.
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
%
% Output
% ======
%   imagex_color                       Gives the [N x N x 3] RGB image with coloring.
%
%
% (C) R. Das, 2004-2012, Thomas Mann, 2012.
% modified by T47, Apr 2013 - May 2013.
%

if nargin == 0;  help( mfilename ); return; end;

[x_image_size, y_image_size, z_image_size] = size(imagex);
if z_image_size ~= 3; fprintf('WARNING: invalid image input.\n'); end;
axis([0 y_image_size 0 x_image_size]); 

if size( which_res, 1) > 1; which_res = which_res'; end;

color_profile_org = [17, 0, mean(what2plot) + 3 * std(what2plot), mean(what2plot) - 3 * std(what2plot)];
if ~exist('color_profie','var') || length(color_profile) < length(color_profile_org);
    color_profile_org(1:length(color_profile)) = color_profile;
    color_profile = color_profile_org;
end;
[color_scheme, d_offset, max_color, min_color] = parse_color_profile(color_profile);
if ~exist('square_width','var'); square_width = 12; end;


imagex_color = double(imagex);
count = 1;
for k = which_res
    x = residue_locations(1, k);
    y = residue_locations(2, k);
    
    colorplot = getcolor( what2plot(count) - d_offset, max_color, min_color, color_scheme);
    
    xbins = max(floor(y-square_width/2),1) : min(floor(y+square_width/2), x_image_size); 
    ybins = max(floor(x-square_width/2),1) : min(floor(x+square_width/2), y_image_size); 
    
    for n = 1:3
      imagex_color(xbins, ybins, n) = double(imagex(xbins, ybins, n)) * colorplot(n);
    end;
    count=count+1;
end;

hold off; image(imagex_color/256); hold on; axis equal; axis off
