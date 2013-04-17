function sc_fc = double_mutant_output(d_align, scale_factor, labels, caption, flg)

% sc_fc = DOUBLE_MUTANT_OUTPUT (d_align, scale_factor, labels, [caption], [r l s])
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
%   [caption]               Title of the figure. Format in string. Optional,
%                                default is 'Double Mutants'.
%   [r l s]                 Flags for xlabel_rotate, make_lines, and SHAPE
%                               reagent. Format in double array: r for whether 
%                               rorate xlabel; l for whether make lines in the 
%                               center; s for SHAPE reagent selection, 0 for NMIA,
%                               1 for 1M7. Optional, default [1, 1, 1].
%
% =Output=
%   [sc_fc]                   scale_factor of figure. Format in double.
%
%
% by T47, Mar 2013.
%

if ~exist('flg');  flg = [1, 1, 1]; end;
if size(labels) ~= size(d_align,2)/2; fprintf(1, 'Label number mismatch!'); end;
if scale_factor == 0; scale_factor = 40/mean(mean(d_align)); end;
if ~exist('caption'); caption = 'Double Mutants';  end;
labels=[repmat({''},1,size(d_align,2)/2), labels];

h = figure(1);
set(h, 'Position', [0, 0, 800, 600]);   

image(d_align * scale_factor); colormap(1 - gray());
make_lines(0:1:size(d_align,2), 'k', 1);
if flg(2); make_lines(size(d_align,2)/2, 'b', 2); end;

title(caption,'FontWeight','Bold','FontSize',20);
set(gca, 'xtick', 1:length(labels), 'xticklabel', labels,'fontsize',12);
if flg(1); xticklabel_rotate; end;

xlb = 'nomod';
if flg(3)
    xla = '1M7';
else
    xla = 'NMIA';
end;
xlabel(strcat([xlb, blanks(54), xla]));
xlabh = get(gca, 'XLabel'); set(xlabh, 'FontSize', 15, 'Color', 'red');
set(xlabh, 'Position', get(xlabh, 'Position') + [0 .25 0]);

set (gcf, 'PaperPosition', [-.25, -.2, 9, 11]);

sc_fc = scale_factor;