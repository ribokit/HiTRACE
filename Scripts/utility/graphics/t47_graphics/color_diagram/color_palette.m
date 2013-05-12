function color_profile = color_palette(d_plot, color_max, color_min, color_scheme)

%
% color_profile = COLOR_PALETTE(d_plot, color_max, color_min, color_scheme)
%
% Plots reactivities with color scheme (color_scheme) and color range of (color_max and
%  color_min) as a preview to optimize color saturation values before actual coloring 
%  the boxes on the image.
%
%
% Input
% =====
%   d_plot          Required        Provides the reactivity data.
%   color_max       Optional        Provides the color saturation value in high-end.
%                                    Default is average plus 3 standard deviation.
%   color_min       Optional        Provides the color saturation value in low-end.
%                                    Default is average minus 3 standard deviation.
%   color_scheme    Optional        Provides the color scheme. See getcolor for more
%                                    details. Default is 17 (rainbow).
%
% Output
% ======
%   color_profile                   Gives the integrated color_profile for following 
%   [color_scheme, d_offset,         coloring steps.
%    max_color, min_color]          'COLOR_SCHEME' specifies the color scheme, same
%                                     as input.
%                                   'D_OFFSET" specifies the vertical offset value all 
%                                     reactivities subtract by to make it centered at 0
%                                     when coloring.
%                                   'MAX_COLOR' specifies the color saturation value at 
%                                     high-end based on d_offset directly used by 
%                                     getcolor when coloring.
%                                   'MIN_COLOR' specifies the color saturation value at 
%                                     low-end based on d_offset directly used by 
%                                     getcolor when coloring.
%                                 
%
% by T47, May 2013.
%

if nargin == 0; help(mfilename); return; end;

if ~exist('color_max','var') || isempty(color_max) || isnan(color_max);
    color_max = mean(d_plot) + 3 * std(d_plot);
end;
if ~exist('color_min','var') || isempty(color_min) || isnan(color_min);
    color_min = mean(d_plot) - 3 * std(d_plot);
end;
if ~exist('color_scheme','var') || isempty(color_scheme); color_scheme = 17; end;


color_center = mean([color_max color_min]);
maxplot = color_max - color_center;
maxplot2 = color_min - color_center;
plot_interval = - (maxplot - maxplot2) / 100;

y_max = max(max(d_plot), color_max);
y_min = min(min(d_plot), color_min);
y_margin = (y_max - y_min) / 4;


figure(1); clf;
set_print_page(gcf, 0, [0 0 800 600], '');
axis([0 length(d_plot) (y_min - y_margin) (y_max + y_margin)]);

for i = y_max:plot_interval:y_min
    color_value = i - color_center;
    make_lines_horizontal(i - 0.5, ...
        getcolor(color_value, maxplot, maxplot2, color_scheme), ...
        1);
    hold on;
end;

plot(d_plot, 'k', 'Linewidth', 2); hold on;
make_lines_horizontal(color_max - 0.5); hold on;
make_lines_horizontal(color_min - 0.5); hold on;
make_lines_horizontal(color_center - 0.5); hold on;
plot(d_plot, 'kd', 'Linewidth', 2); hold off;


figure(2); clf;
set_print_page(gcf, 0, [100 100 800 600], '');
axis([0 length(d_plot) (y_min - y_margin) (y_max + y_margin)]);



color_profile = [color_scheme, color_center, maxplot, maxplot2];
