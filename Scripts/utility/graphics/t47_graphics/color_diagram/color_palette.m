function color_profile = color_palette(d_plot, color_max, color_min, color_scheme, seqpos, sequence)

%
% color_profile = COLOR_PALETTE(d_plot, color_max, color_min, color_scheme, seqpos, sequence)
%
% Plots reactivities with color scheme (color_scheme) and color range of (color_max and
%  color_min) as a preview to optimize color saturation values before actual coloring
%  the boxes on the image.
% Two figures will be shown for preview: a background-colored plot and a foreground-
%  colored bargraph. Black lines denote chosen color_max and color_min as well center
%  of color range.
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
%   seqpos          Optional        Provides the seqpos for x-axis labeling. Default 
%                                    is [], means natural numbering will be used.
%   sequence        Optional        Provides the sequence for x-axis labeling. Default 
%                                    is '', means natural numbering will be used.
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

if ~exist('sequence','var') || isempty(sequence) || ...
        ~exist('seqpos','var') || isempty(seqpos);
    is_x_label = 0;
else
    is_x_label = 1;
end;

color_center = mean([color_max color_min]);
maxplot = color_max - color_center;
maxplot2 = color_min - color_center;
plot_interval = - (maxplot - maxplot2) / 100;

y_max = max(max(d_plot), color_max);
y_min = min(min(d_plot), color_min);
y_margin = (y_max - y_min) / 4;



figure; clf;
set_print_page(gcf, 0, [0 0 800 600], 'Color Palette Plot');

% figure of background-colored plot-graph
subplot(2,1,1);
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

% figure of foreground-colored bar-graph
subplot(2,1,2);

h = bar(d_plot);
h_child = get(h, 'Children');

% custom colormap from getcolor
my_color = zeros(length(d_plot), 3);
for i = 1:length(d_plot)
    my_color(i, :) = getcolor(d_plot(i) - color_center, maxplot, maxplot2, color_scheme);
end;
colormap(my_color);
set(h_child, 'CData', 1:length(d_plot));
set(h_child, 'EdgeColor', 'none');
set(h_child, 'CDataMapping', 'direct');

axis([0 length(d_plot) (y_min - y_margin) (y_max + y_margin)]);
make_lines_horizontal(color_max - 0.5); hold on;
make_lines_horizontal(color_min - 0.5); hold on;
make_lines_horizontal(color_center - 0.5); hold off;

% X-axis-tick, seqpos
if is_x_label;
    name = cell(1, length(seqpos));
    for i = seqpos
        name{i} = [num2str(i) sequence(i)];
    end;
    set(gca, 'FontSize', 8);
    set(gca, 'XTick', 1:length(seqpos), 'XTickLabel', char(name));
    xticklabel_rotate();
end;

color_profile = [color_scheme, color_center, maxplot, maxplot2];
