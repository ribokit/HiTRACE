function print_xsel_split(d_align, xsel, seqpos, sequence, offset, area_pred, labels, num_flg, bol_flg, str_flg, ft_sz, clr_cd)

%
% PRINT_XSEL_SPLIT(d_align, xsel, seqpos, sequence, offset, area_pred, labels...
%               [number_flag], [boolean_flag], [string_flag], [font_size], [color_code])
%
%
% Prints the band annotation by lines and circles on electrophoregram in H x W splitted
%  figures, with X-axis labelled with lane names, Y-axis labelled with nucleotide
%  positions.
%
%
% Arguments
% =========
%   d_align         Required    Provides the electrophoregram matrix. Format in double
%                                   array.
%   xsel            Required    Provides annotation Y position. Format in double array.
%                                   Should be strict monotonicly increasing, decreasing
%                                   array will be flipped automatically.
%   seqpos          Required    Provides bands annotation array. Format in double arrat.
%                                   Should be strictly monotonic decreasing, increasing
%                                   array will be flipped automatically.
%   sequence        Required    Provides sequence. Format in string.
%   offset          Required    Provides offset. Format in double.
%   area_pred       Optional    Provides the peak position prediction  matrix to put red
%                                   circles on the plot. Format in double array. Default
%                                   is [] (no circles).
%   labels          Optional    Provides the lane lables for x-axis. Format in string
%                                   cell. Default is {} (no label).
%   [number_flag]   Optional    Provides the layout format that will split into H x W
%   [H  W  SC                       pages with vertical blank edge SP, using SC for
%    SP UO LO                       intensity scaling. D_ALIGN will be auto-trimmed by
%    L1S L1I L1E                    UO and LO if asked, two sets of lines will be added
%    L2S L2I L2E]                   at L1S:L1I:L1E and L2S:L2I:L2E. Format in double
%                                   array, default is [0 0 0 50 75 100 0 0 0 0 0 0].
%                               'H' (height) denotes the number of pages vertically,
%                                   0 means use auto_size if AS in BOOLEAN_FLAG;
%                               'W' (width) denotes the number of pages horizontally,
%                                   0 means use auto_size if AS in BOOLEAN_FLAG;
%                               'SC' (scale_factor) denotes the scaling factor of image,
%                                   default is calculated to display optimally according
%                                   to sc = 27.5 / mean(mean(d_align));
%                               'SP' (spacer) denotes the spacer width on both top and 
%                                   bottom of each sub-figure;
%                               'UO' (upper_bound) denotes the excess D_ALIGN image
%                                   shown above the first SEQPOS label;
%                               'LO' (lower_bound) denotes the excess D_ALIGN image
%                                   shown under the last SEQPOS label;
%                               'L1S:L1I:L1E' (line_1 start:interval:end) denotes the
%                                   positions for separating line set 1. All 0 means
%                                   not drawing;
%                               'L2S:L2I:L2E' (line_2 start:interval:end) denotes the
%                                   positions for separating line set 2. All 0 means
%                                   not drawing.
%   [boolean_flag]  Optional    Provides the layout format that whether one-page only, 
%   [OP AT AS                       whether auto-trim vertical boundaries, whether auto 
%    AL AC PR                       decide overall size, whether auto decide text length 
%    TL VR PN]                      for titles, whether auto adjust contrast, whether  
%                                   print to file, whether add title, version label and
%                                   page number. Format in  double  array, default is 
%                                   [0 1 1 1 0 1 1 1 1].
%                               'OP' (is_one_page) denotes whether use compact one-page 
%                                   mode. Instead of multi-pages, it will arrange by 
%                                   multi-subplots. y-axis labels will only be shown on 
%                                   right-side, W will be forced to 1;
%                               'AT' (is_auto_trim) denotes whether to auto trim top and
%                                   bottom of D_ALIGN for optimal display. Excess image
%                                   size is denoted by UO and LO in NUMBER_FLAG;
%                               'AS' (is_auto_size) denotes whether to auto decide number
%                                   of pages based on D_ALIGN and SEQUENCE. This will
%                                   override H and W in NUMBER_FLAG;
%                               'AL' (is_auto_length) denotes whether auto-decide font
%                                   size of titles to fit in page margins. This will
%                                   override T1, T2, and T3 in FONT_SIZE;
%                               'AC' (is_auto_contrast) denotes whether auto-increase
%                                   contrast between vertical panels, an increment of 
%                                   1.25^H_num will be applied for visual attenuation
%                                   alleviation, i.e. [1, 1.25, 1.5625, 1.9531]x;
%                               'PR' (if_print) denotes whether to print to .eps files,
%                                   prints will be saved in folder DN of STRING_FLAG;
%                               'TL' (is_title) denotes whether title is added to figure;
%                               'VR' (is_version) denotes whether version label is
%                                   added to bottom right corner;
%                               'PN' (is_page_number) denotes whether page number is
%                                   added to the corners of each figure;
%                               1 equals TRUE; 0 equals FALSE.
%   [string_flag]   Optional    Provides the title name, string input for output .eps file
%   [TI FN DN                       name, folder name, authorship and date onp rintout.
%    AU DT VR                       Also provides font weight options. Format in string 
%    T1 T2 VE                       cell, default is {'', '', 'print_xsel_split_output', ...
%    YT XT]                         '', 'mmm yyyy', '', 'Bold', 'Bold', 'Normal', ...
%                                   'Normal', 'Normal'}.
%                               'TI' (title_name) denotes the title to display;
%                               'FN' (file_name) denotes file name for print files;
%                                   Numbers, underscore, and '.eps' extension will
%                                   automatically append;
%                               'DN' (dir_name) denotes folder name for all files.
%                               'AU' (author_name) denotes authorship series. 'by' and
%                                   '@' will automatically append;
%                               'DT' (date_string) denotes date string appearing on top
%                                   right corner. Default will use current date;
%                               'VR' (hitrace_ver) denotes current HiTRACE subversion;
%                               'T1' font weight of title;
%                               'T2' font weight of date and author label;
%                               'VE' font weight of version label;
%                               'YT' font weight of y-axis tick label;
%                               'XT' font weight of x-axis tick label;
%   [size_flag]     Optional    Provides font size and line width values for figures. Format 
%   [T1 T2 VE                       in double array, default is [25 15 9 8 11 1 2 2].
%    YT XT                      'T1' font size of title;
%    LB L1 L2]                  'T2' font size of date and author label;
%                               'VE' font size of version label;
%                               'YT' font size of y-axis tick label;
%                               'XT' font size of x-axis tick label;
%                               'LB' line width of borders on each page;
%                               'L1' line width of separating line set 1;
%                               'L2' line width of separating line set 2;
%   [color_flag]    Optional    Provides color codes for figures. Format in string, default 
%   [T1 T2 VE                       is 'kbkkkymb'.
%    YT XT                      'T1' font color of title;
%    LB L1 L2]                  'T2' font color of date and author label;
%                               'VE' font color of version label;
%                               'YT' font color of y-axis tick label;
%                               'XT' font color of x-axis tick label;
%                               'LB' line color of borders on each page;
%                               'L1' line color of separating line set 1;
%                               'L2' line color of separating line set 2;
%
%
% e.g. PRINT_XSEL_SPLIT(d_align, xsel, seqpos, sequence, offset, area_pred);
%      PRINT_XSEL_SPLIT(d_align, xsel, seqpos, sequence, offset, area_pred, labels, ...
%                       [3], [1 0], {'', '', 'T47', 'Apr 2013'}, [], {'r', 'b', 'k'});
%
%
% Notes
% =====
% Open all .eps files in single Preview window and print using a color printer.
%   Deselect auto-rotate and auto-scale in the print dialog, make sure scale
%   is 100%. Spaces on each border are included for easy splicing. Cut at yellow
%   lines on each page, then put together the entire image.
%
% ALL current opened figures will be LOST! Save before run.
%
%
% by T47, May 2013.
%

if nargin == 0; help( mfilename ); return; end;
Script_VER = '1.3';

%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparation before splitting
fprintf(['print_xsel_split version ', print_str(Script_VER), '.\n\n']);

fprintf(['d_align (',num2str(size(d_align,1)),' x ',num2str(size(d_align,2)),') provided by user.\n']);
fprintf(['xsel (1 x ',num2str(length(xsel)),') provided by user.\n']);
fprintf(['seqpos (1 x ',num2str(length(seqpos)),') provided by user.\n']);
fprintf(['sequence (1 x ',num2str(length(sequence)),') provided by user.\n']);

if ~exist('offset','var') || isempty(offset);
    offset = 0;
    fprintf('offset (0) generated by default.\n');
else
    fprintf(['offset (',num2str(offset),') provided by user.\n']);
end;

if ~exist('area_pred','var') || isempty(area_pred);
    area_pred = zeros(length(seqpos), size(d_align, 2));
    fprintf(['area_pred (',num2str(size(area_pred,1)),' x ',num2str(size(area_pred,2)),') generated by default.\n']);
else
    fprintf(['area_pred (',num2str(size(area_pred,1)),' x ',num2str(size(area_pred,2)),') provided by user.\n']);
end;

if ~exist('labels','var') || isempty(labels);
    labels = repmat({''}, 1, size(d_align, 2));
    fprintf(['labels (1 x ',num2str(length(labels)),') generated by default.\n']);
elseif length(labels) < size(d_align, 2);
    labels = [labels, repmat({''}, 1, size(d_align, 2) - length(labels))];
    fprintf(['labels (1 x ',num2str(length(labels)),') filled up by default.\n']);
else
    fprintf(['labels (1 x ',num2str(length(labels)),') provided by user.\n']);
end;
for i = 1:length( labels ); if iscell( labels{i} ); labels{i} = char(labels{i}); end; end;
offset = round(offset);
fprintf('\n'); fprintf(['sequence  ', sequence,'\n\n']);

if ~exist('num_flg','var'); num_flg = []; end;
num_flg = is_valid_flag(num_flg, [0 0 0 50 75 100 0 0 0 0 0 0]);
num_flg([1:2, 4:7]) = round(num_flg([1:2, 4:7]));
num_flg(1) = max([num_flg(1), 0]);
num_flg(2) = max([num_flg(2), 0]);
num_flg(3) = max([num_flg(3), 0]);
if num_flg(4) < 1; num_flg(4) = 50; end;
page_num_H = num_flg(1); page_num_W = num_flg(2); scale_fc = num_flg(3);
num_sp = num_flg(4); num_up_offset = num_flg(5); num_low_offset = num_flg(6);
num_l1_s = num_flg(7); num_l1_i = num_flg(8); num_l1_e = num_flg(9);
num_l2_s = num_flg(10); num_l2_i = num_flg(11); num_l2_e = num_flg(12);

if ~exist('bol_flg','var'); bol_flg = []; end;
bol_flg = is_valid_flag(bol_flg, [0 1 1 1 0 1 1 1 1]);
for i = length(bol_flg)
    bol_flg(i) = is_valid_boolean(bol_flg(i));
end;
is_one_page = bol_flg(1);
is_auto_trim = bol_flg(2); is_auto_size = bol_flg(3); is_auto_length = bol_flg(4); is_auto_contrast = bol_flg(5);
is_print = bol_flg(6); is_title = bol_flg(7); is_version = bol_flg(8); is_page_no = bol_flg(9);

if ~exist('str_flg','var'); str_flg = {}; end;
str_flg = is_valid_flag(str_flg, {'', '', 'print_xsel_split_output', '', datestr(date, 'mmm yyyy'), '', ...
    'Bold', 'Bold', 'Normal', 'Normal', 'Normal'});
title_name = str_flg{1}; file_name = str_flg{2}; dir_name = str_flg{3};
author_str = str_flg{4}; date_str = [' @ ' str_flg{5}]; ver_hitrace = str_flg{6};
ft_w_title_1 = str_flg{7}; ft_w_title_2 = str_flg{8};
ft_w_ver = str_flg{9};
ft_w_y_tick = str_flg{10}; ft_w_x_tick = str_flg{11};
if ~isempty(file_name); file_name = [file_name '_']; end;
if isempty(dir_name); dir_name = 'print_xsel_split_output'; end;
if ~isempty(author_str); author_str = [' by ' author_str]; end;
if strcmp(date_str, ' @ '); date_str = datestr(date, 'mmm yyyy'); end;
if isempty(ver_hitrace); ver_hitrace = 'N/A'; end;

if ~exist('ft_sz','var'); ft_sz = []; end;
ft_sz = is_valid_flag(ft_sz, [25 15 9 8 11 1 2 2]);
ft_sz_title_1 = ft_sz(1); ft_sz_title_2 = ft_sz(2);
ft_sz_ver = ft_sz(3); ft_sz_y_tick = ft_sz(4); ft_sz_x_tick = ft_sz(5);
ln_wt_border = ft_sz(6); ln_wt_1 = ft_sz(7); ln_wt_2 = ft_sz(8);

if ~exist('clr_cd','var'); clr_cd = ''; end;
clr = parse_color_string(clr_cd);
clr = is_valid_flag(clr, {'k', 'b', 'k', 'k', 'k', 'y', 'm', 'b'});
for i = 1:length(clr)
    if clr{i} == 'g'; clr{i} = [0 0.5 0]; end;
end;
color_title_1 = clr{1}; color_title_2 = clr{2};
color_ver = clr{3}; color_y_tick = clr{4}; color_x_tick = clr{5};
color_line_border = clr{6}; color_line_1 = clr{7}; color_line_2 = clr{8};


% auto_size if asked
is_auto_size_force = (page_num_H == 0) && (page_num_W == 0);
[page_num_H, page_num_W] = auto_size(page_num_H, seqpos, page_num_W, size(d_align, 2), is_auto_size);

% print out summary
fprintf(['Divide into ', print_str(page_num_H)]);
if is_one_page; 
    fprintf(' panels');
else
    fprintf([' x ', print_str(page_num_W), ' pages']);
end;
fprintf([', one page (', print_yn(is_one_page), '), auto size (', ...
    print_yn(is_auto_size || is_auto_size_force), '), auto contrast (', print_yn(is_auto_contrast), ').\n']);
fprintf(['Spacer ', print_str(num_sp), ', make lines (']);
if all([num_l1_s, num_l1_i, num_l1_e, num_l2_s, num_l2_i, num_l2_e] == 0);
    fprintf([print_yn(0), ').\n']);
else
    fprintf([print_str(num_l1_s), ':', print_str(num_l1_i), ':', print_str(num_l1_e), ' and ', ...
        print_str(num_l2_s), ':', print_str(num_l2_i), ':', print_str(num_l2_e), ').\n']);
end;
fprintf(['Show title (', print_yn(is_title), '), show page number (', print_yn(is_page_no), ') show version (', ...
    print_yn(is_version), '), print to file (', print_yn(is_print), ').\n']);

% auto_trim if asked
fprintf(['Auto trim (', print_yn(is_auto_trim), ')']);
[d_align, xsel] = auto_trim(d_align, xsel, num_up_offset, num_low_offset, is_auto_trim);
fprintf('.\n');
fprintf(['Auto title size (', print_yn(is_auto_length), ').\n\n']);

% calculate auto_scale
auto_scale = 27.5 / mean(mean(d_align));
if scale_fc == 0; scale_fc = auto_scale; end;
fprintf( 'auto_scale_factor = %f\n', auto_scale);
fprintf( 'scale_used = %f\n', scale_fc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% split d_align into h*w sub matrices
% expand d_align to divisible size by h and w
if is_one_page; page_num_W = 1; end;
[d_align, h_length, w_length] = d_expand_divisible(d_align, page_num_H, page_num_W);

% add spacers on each border of each figure
d = d_split_add_spacer(d_align, page_num_H, h_length, num_sp, page_num_W, w_length, 1);

% split mutants names array into w*w_length array
name = cell(page_num_W, (w_length + 2));
for i = 1:page_num_W
    name(i, 2:(w_length + 1)) = labels((w_length * (i - 1) + 1):(w_length * i));
end;


% flip seqpos and xsel to correct order
% xsel be increasing, seqpos be decreasing, mutpos be increasing
fprintf('\n');
xsel = auto_flip_monotone(xsel, 1, 'xsel');
seqpos = auto_flip_monotone(seqpos, -1, 'seqpos');
fprintf('\n');

% read in band names (Y-axis)
% split band position array into h*h_length array
[bandpos, band] = xsel_label_split(xsel, seqpos, sequence, offset, page_num_H, h_length, num_sp);

% prepare area_pred dimensions and flip to right orientation
h_l = get_figure_h_array(band, page_num_H);
area_pred = flipud(area_pred);

% pause point
fprintf(['In each page, there are ',print_str(w_length),' lanes (X-axis), and ',...
    print_str(h_length),' of traces (Y-axis):\n']);
fprintf([' <strong>', strrep(num2str(h_l), '  ', ', ') ,'</strong> sequence positions on each page row.\n']);
fprintf('\n'); fprintf(2,'Press any key to continue...\n');
pause;



%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot each figure
close all;
if is_one_page;
    figure();
    set(gcf, 'PaperOrientation', 'Portrait', 'PaperPositionMode', 'Manual', ...
        'PaperSize', [8.5 11], 'Color', 'White');
    set(gcf, 'PaperPosition', [0 0 8.5 11]);
    set(gcf, 'Position', [0, 0, 600, 800]);
end;

for i = 1:page_num_W
    for j = 1:page_num_H
        
        % extract y-axis information from cell
        d_temp = [];
        d_temp(1:size(d, 3), 1:size(d, 4)) = d(i, j, :, :);
        [xsel_pos_sub, xsel_txt_sub] = parse_xsel_label(band, j, 2);
        
        % add figure number to the top left and bottom right corner
        if is_one_page;
            fig_name = ['[', num2str(j), ']'];
        else
            fig_name = ['[', num2str(j), ', ', num2str(i), ']'];
        end;
        if is_page_no == 1;
            xsel_pos_sub(1) = round(num_sp / 2);
            xsel_txt_sub{1} = fig_name;
            name{i, w_length + 2} = fig_name;
        else
            xsel_pos_sub(1) = 1;
            xsel_txt_sub{1} = '';
            name{i, w_length + 2} = '';
        end;
        name(i) = fill_space_label(name(i), 0);
        xsel_txt_sub = fill_space_label(xsel_txt_sub, 1);
        
        % setup subplot panels for one_page
        if is_one_page;
            h = subplot(1, page_num_H, j);
            fig_num = j;
        else
            fig_num = (i - 1) * page_num_H + j;
            h = figure(fig_num);
        end;
        
        % adjustment for auto_title_length
        if ~is_one_page;
            fig_w = 600; fig_h = 800; 
            set(gcf, 'PaperOrientation', 'Portrait', 'PaperPositionMode', 'Manual', ...
                'PaperSize', [8.5 11], 'Color', 'White');
            set(gcf, 'PaperPosition', [0 0.5 8.5 10]);
            set(h, 'Position', [(fig_num - 1) * 20, 0, fig_w, fig_h]);
            title_h = 0.45;
        else
            title_h = 0.35;
        end;
        title_size_fc = 0.7;  title_offset = 0;
        
        % auto-contrast scale_factor adjustment
        scale_fc_h = scale_fc;
        if is_auto_contrast; scale_fc_h = scale_fc*1.25^(j-1); end;
        image(d_temp * scale_fc_h); hold on;
        colormap(1 - gray());
        
        % make lines
        % separating line set 1 and 2
        make_lines([], 'k', 0.5); hold on;
        if ~all([num_l1_s, num_l1_i, num_l1_e] == 0);
            make_lines((num_l1_s + 1):num_l1_i:(num_l1_e + 1), color_line_1, ln_wt_1);
        end;
        if ~all([num_l2_s, num_l2_i, num_l2_e] == 0);
            make_lines((num_l2_s + 1):num_l2_i:(num_l2_e + 1), color_line_2, ln_wt_2);
        end;
        hold on;
        
        % yellow border lines
        make_lines_horizontal(num_sp - 0.5, color_line_border, ln_wt_border);
        make_lines_horizontal(size(d_temp,1) - num_sp - 0.5, color_line_border, ln_wt_border);
        make_lines(1, color_line_border, ln_wt_border);
        make_lines(size(d_temp, 2)-1, color_line_border, ln_wt_border);
        hold on;
        
        % axis labeling
        % X-axis-tick, mutant names
        set(gca, 'FontSize', ft_sz_x_tick, 'FontWeight', ft_w_x_tick);
        set(gca, 'XTick', 1:size(name, 2), 'XTickLabel', char(name{i, :}), 'XColor', color_x_tick);
        xticklabel_rotate();
        hold on;
        
        % Y-axis-tick, band positions
        set(gca, 'FontSize', ft_sz_y_tick, 'FontWeight', ft_w_y_tick, 'YColor', color_y_tick);
        set(gca, 'YTick', xsel_pos_sub(:), 'YTickLabel', char(xsel_txt_sub(:)), 'YAxisLoc', 'Right');
        hold on;
        
        % make color lines
        for k = 2:length(xsel_pos_sub)
            res_lbl = char(xsel_txt_sub(k)); res_id = res_lbl(1);
            res_clr = get_residue_color(res_id);
            make_lines_horizontal(xsel_pos_sub(k) - 0.5, res_clr, 0.5);
        end;
        hold on;
        
        % make circles based from area_pred
        area_pred_sub = area_pred( ...
            (sum(h_l(1:j)) - h_l(j) + 1): sum(h_l(1:j)), ...
            (w_length * (i - 1) + 1): w_length * i);
        xsel_sub = xsel_pos_sub(xsel_pos_sub ~= 0);
        xsel_sub = xsel_sub(2:end);
        for m = 1:size(area_pred_sub, 1)
            mark_points = find( area_pred_sub(m, :) > 0.5 );
            if ~isempty(mark_points);
                plot(mark_points + 1, xsel_sub(m), 'ro');
            end;
            hold on;
        end;
        
        % version label of right-bottom corner
        if (is_version && j == page_num_H && i == page_num_W)
            ver_str = ['HiTRACE rev. ', num2str(ver_hitrace), ' / print\_xsel\_split ver. ', Script_VER, '   '];
            ylim = get(gca, 'YLim') - is_one_page*10; xlim = get(gca, 'XLim');
            text(xlim(2), ylim(2), ver_str, 'HorizontalAlignment', 'Right', 'VerticalALignment', 'Bottom', ...
                'FontWeight', ft_w_ver, 'FontSize', ft_sz_ver, 'Color', color_ver);
        end;

        % second y-axis labels
        if ~is_one_page;
            y_1 = gca; y_2 = copyobj(y_1, gcf);
            set(y_2, 'YAxisLoc', 'Left');
            hold on;
        end;
        
        % add title if asked
        if (is_title && i == page_num_W && ((~is_one_page && j == 1) || (is_one_page && j == page_num_H)) );
        
            comment = [author_str, date_str];
            xlim = get(gca, 'XLim');
            txt = text(xlim(2), 0, comment);
            set(txt, 'HorizontalAlignment', 'Right', 'VerticalAlignment', 'Bottom',...
                'FontWeight', ft_w_title_2, 'FontSize', ft_sz_title_2, 'FontName', 'Courier', 'Color', color_title_2);
            if is_auto_length; auto_font_size(txt, min((get(gcf,'PaperSize'))) * title_size_fc * 0.2, title_h / 2); end;
           
            title([' ' title_name ' '], 'HorizontalAlignment', 'Left', 'VerticalAlignment', 'Bottom',...
                'FontWeight', ft_w_title_1, 'FontSize', ft_sz_title_1, 'FontName', 'Courier', 'Color', color_title_1,'interp','none');
            tit = get(gca, 'Title'); pos = get(tit, 'Position');
            
            if ~is_one_page;
                set(tit, 'Position', [0 (pos(2) + title_offset) pos(3)]);
                tit_adjst = 0.8;
            else
                tit_x_adjst = [0 0.22 0.45 0.65 0.9];
                set(tit, 'Position', [-min((get(gcf,'PaperSize'))) * title_size_fc * tit_x_adjst(page_num_H) (pos(2) + title_offset) pos(3)]);
                tit_adjst = 0.6;
            end;
            if is_auto_length && ~isempty(title_name); auto_font_size(tit, min((get(gcf,'PaperSize'))) * title_size_fc * tit_adjst, title_h); end;
        end;
        hold off;
        
        % print to file if asked
        if is_print && ~is_one_page; print_save_figure(h, [file_name, num2str(fig_num)], dir_name, 0); end;
    end;
end;

if is_print; 
    if is_one_page;
        print_save_figure(gcf, [file_name, 'all'], dir_name, 0);
        fprintf(['1 page printed to folder "', dir_name,'".\n']);
    else
        fprintf([num2str(i*j),' pages printed to folder "', dir_name,'".\n']);
    end;
end;



%%%%%%
function str = print_yn (flag)

str = ['<strong>',num2yn(flag),'</strong>'];

function str = print_str (flag)

str = ['<strong>',num2str(flag),'</strong>'];