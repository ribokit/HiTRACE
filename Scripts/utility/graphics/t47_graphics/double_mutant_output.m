function double_mutant_output(d_align, scale_factor, labels, caption, flg, ft_sz, clr, filename)

% DOUBLE_MUTANT_OUTPUT (d_align, scale_factor, labels, [caption], ...
%                       [boolean_flag], [font_size], [color], [filename])
%
% Generates comparative output of electrophoregrams, with nomod lanes on
%   the left and SHAPE lanes on the right.
%
% =Input=
%   d_align                 Data matrix after align_by_DP or align_by_DP_fine,
%                               trimmed to best show. nomod data on the left and 
%                               SHAPE data on  the right. Format in double array.
%   scale_factor            Scales for matrix image. If 0, optimal scale is 
%                               calculated by auto_scale. Format in double.
%   labels                  X-tick labels for annotation. Only SHAPE lanes will
%                                be labeled, no repeat or space needed for input.
%                                Format in string cell, e.g. {'WT', 'A10C'}.
%   [caption]               Title of the figure. 'Double Mutants' will automatically 
%                                append. Format in string. Optional, default is ''.
%   [boolean_flag]          Flags for xlabel_rotate, make_lines, and SHAPE reagent. 
%   [R L S K P]                 Format in double array, optional, default is
%                               [1, 1, 1, 1, 0].
%                           'R' (if_rotate) for whether rorate xlabel;
%                           'L' (if_lines) for whether make lines in the center;
%                           'S' (what_modifier) for reagent selection, 0 for NMIA, 
%                               1 for 1M7, 2 for DMS, 3 for CMCT;
%                           'K' (if_xlabel) for whether display 'nomod' 'modifier'
%                               label.
%                           'P' (if_print) forwhether print to .eps file.  
%                           0 for FALSE, 1 for TRUE.
%   [font_size]             Font sizes for all labels in the figure, and offset 
%   [T XT MD SP OFS]            for label positions. Format in double array, 
%                               optional, default is [20 12 15 54 .25].
%                           'T' (size_title) for caption font size.
%                           'XT' (size_xtick_label) for lane name font size.
%                           'MD' (size_modifier) for modifier display on x-axis.
%                           'SP' (space_modifier) for spacing between x-label
%                               'nomod' and 'modifier'.
%                           'OFS' (offset_modifier) for distance between x-label 
%                               and x-axis.
%   [color]                 Color codes for labels and lines. Format in string 
%   [T XT MD L1 L2]             cell, optional, default is {'k', 'k', 'r', 'b', 'k'}.
%                           'T' (color_title) for caption font color.
%                           'XT' (color_xtick_label) for lane name font color.
%                           'MD' (color_modifier) for modifier font color.
%                           'L1' (color_line_center) for center separating
%                               line color.
%                           'L2' (color_line_each) for the rest vertical lines.
%   [filename]              Filename for print pages. '.eps' extension will    
%                               automatically append. Format in string. Optional, only
%                               takes effect when if_print is 1; default is 'dbl_mnt'.
%
% by T47, Mar 2013.
%

if ~exist('caption','var') || isempty('caption'); caption = ''; end;
if ~exist('flg','var') || isempty('flg') || length(flg) < 5; flg = [1, 1, 1, 1, 0]; end;
if ~exist('ft_sz','var') || isempty('ft_sz') || length(ft_sz) < 5; ft_sz = [20 12 15 54 .25]; end;
if ~exist('clr','var') || isempty('clr') || length(clr) < 5; clr = {'k','k','r','b','k'}; end;
if ~exist('filename','var') || isempty('filename'); filename = 'dbl_mnt'; end;

% check if label size match d_align
if size(labels) ~= size(d_align,2)/2; fprintf('Label number mismatch!\n'); end;

% calculates auto scale and inform user
if scale_factor == 0; scale_factor = 40/mean(mean(d_align)); end;
fprintf('\n');
fprintf(['scale_factor_used = ', num2str(scale_factor),'\n']);

% append caption and generate blanks for xtick
caption = [caption ' Double Mutants'];
labels = [repmat({''},1,size(d_align,2)/2), labels];

h = figure(1);
set(h, 'Position', [0, 0, 600, 800]);   
set(gcf, 'PaperOrientation', 'portrait', 'PaperPositionMode', 'auto', 'color', 'white');
image(d_align * scale_factor); colormap(1 - gray());

% make lines if asked
make_lines(0:1:size(d_align,2), clr{5}, 1);
if flg(2); make_lines(size(d_align,2)/2, clr{4}, 2); end;

% add title and xtick
title(caption,'FontWeight','Bold','FontSize',ft_sz(1),'Color', clr{1});
set(gca, 'xtick', 1:length(labels), 'xticklabel', labels,'fontsize',ft_sz(2),'color',clr{2});
if flg(1); xticklabel_rotate; end;

% add xlabel of modifier and set position if asked
if flg(4);
    if flg(3) == 0; mdfr = 'NMIA'; end;
    if flg(3) == 1; mdfr = '1M7'; end;
    if flg(3) == 2; mdfr = 'DMS'; end;
    if flg(3) == 3; mdfr = 'CMCT'; end;
    
    xlabel(['nomod', blanks(ft_sz(4)), mdfr]);
    xlabh = get(gca, 'XLabel'); set(xlabh, 'FontSize', ft_sz(3), 'Color', clr{3});
    set(xlabh, 'Position', get(xlabh, 'Position') + [0 ft_sz(5) 0]);
end;

% print to file if asked
if flg(5); print(h,'-depsc2',[filename,'.eps']); end;


