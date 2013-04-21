function plot_1d_color_ref (seqpos, sat_array, dil_array, ref_array, axs, mdfr, shadow_lines, ref_lines, caption, if_print, filename, ft_sz, line_width, clr)

% PLOT_1D_COLOR_REF(seqpos, sat_array, dil_array, ref_array, [axis], [modifier], ...
%                    [shadow_lines], [ref_lines], [title], [if_print], [filename], ...
%                    [font_size], [line_width], [color])
%
% Plots normalized reactivity of saturated lane, diluted lane and referencing corrected
%  lane together for comparison. Flanking sequences and referencing positions are 
%  marked in different shadow color.
%
% Arguments
% =========
%   seqpos               Required        Provides the sequence position for x-axis plotting.
%                                            Format in double array.
%   sat_array            Required        Provides the normalized reactivity calculated from 
%                                            saturated sample. Format in double array.
%   dil_array            Required        Provides the normalized reactivity calculated from 
%                                            diluted sample. Format in double array.
%   ref_array            Required        Provides the normalized reactivity calculated by 
%                                            referencing correction. Format in double array.
%   [axis]               Optional        Provides the axis dimensions for plot. Format in 
%                                            double array. Default is based on data range.
%   [modifier]           Optional        Provides modifier name. Format in string. Default 
%                                            is blank.
%   [shadow_lines]       Optional        Provides the region for shadowing flanking sequences.
%                                            Use numbering in reference with seqpos provided 
%                                            including offset. Format in double array. Default 
%                                            is blank.
%   [ref_lines]          Optional        Provides the position of referencing peaks chosen. 
%                                            Use numbering in reference with seqpos provided 
%                                            including offset. Format in double array. Default 
%                                            is blank.
%   [title]              Optional        Provides title for plot. Format in string. Default 
%                                            is blank.
%   [if_print]           Optional        Provides whether print to file. Default is 1 (YES).                                 
%   [filename]           Optional        Provides the filename for output images. Format in 
%                                            string. Default is 'plot.png'. Modifier name and 
%                                            '.eps' extension will automatically append.
%   [font_size]          Optional        Provides font size values for plot. Format in double 
%   [T L]                                    array Default is [25 20].
%                                            'T' (size_title) for title font size.
%                                            'L' (size_legend) for legend font size.
%   [line_width]         Optional        Provides line width values for plot. Default is
%   [S D R H F]                              [2 2 2 1 2].
%                                            'S' (size_saturated) for saturated reactivity.
%                                            'D' (size_diluted) for diluted reactivity.
%                                            'R' (size_reference) for corrected reactivity.
%                                            'H' (size_shadow) for shadowing region.
%                                            'F' (size_ref_pos) for referencing peaks.
%   [color]              Optional        Provides color code for plot. Format in string cell. 
%   [S D R H P T]                            Default is {'b', 'r', 'k', 'k', 'g', 'k'}.
%                                            'S' (color_saturated) for saturated reactivity.
%                                            'D' (color_diluted) for diluted reactivity.
%                                            'R' (color_reference) for corrected reactivity.
%                                            'H' (color_shadow) for shadowing region.
%                                            'F' (color_ref_pos) for referencing peaks.
%                                            'T' (color_title) for title.
%
% Notes
% =====
% All current opened figures will be lost. Save before run.
%
% by T47, Apr 2013
%

if ~exist('axs','var'); axs = []; end;
if ~exist('mdfr','var'); mdfr = ''; end;
if ~exist('shadow_lines','var'); shadow_lines = []; end;
if ~exist('ref_lines','var'); ref_lines = []; end;
if ~exist('caption','var'); caption = ''; end;
if ~exist('if_print','var'); if_print = 1; end;
if ~exist('filename','var') || isempty(filename); filename = 'plot'; end;
if ~exist('ft_sz','var') || isempty(ft_sz) || length(ft_sz) < 2; ft_sz = [25 20]; end;
if ~exist('line_width','var') || isempty(line_width) || length(line_width) < 5; line_width = [2 2 2 1 2]; end;
if ~exist('clr','var') || isempty(clr) || length(clr) < 6; clr = {'b', 'r', 'k', 'k', 'g', 'k'}; end;

close all;
h = figure; 
set(h,'Position',[0 0 800 600]);
set(gcf, 'PaperOrientation', 'landscape', 'PaperPositionMode', 'auto', 'color', 'white');

% plot three traces
plot(seqpos, sat_array, clr{1}, 'linewidth', line_width(1));hold on;
plot(seqpos, dil_array, clr{2}, 'linewidth', line_width(2));hold on;
plot(seqpos, ref_array, clr{3}, 'linewidth', line_width(3));hold on;

% put title and legend
if ~isempty(axs); axis(axs); end;
title(caption,'FontWeight','Bold','FontSize',ft_sz(1),'Color', clr{6});
set(gca, 'fontsize', ft_sz(2), 'fontweight', 'bold');
legend([mdfr, '-saturated'], [mdfr, '-diluted'], [mdfr, '-referenced']);

% make shadow and ref lines
plot(seqpos,zeros(length(sat_array)), 'k');hold on;
if ~isempty(shadow_lines);
    make_lines(shadow_lines - 0.5, clr{4}, line_width(4), 1, 0); hold on;
end;
if ~isempty(ref_lines)
    make_lines(ref_lines - 0.5, clr{5}, line_width(5), 1, 0); hold on;
end;

% print to file if asked
if if_print; print(h,'-depsc2',[filename, '_', mdfr, '.eps']); end;

