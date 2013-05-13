function print_CE_split(d_align, d_rdat, seqpos, mutpos, xsel, sequence, offset, num_flg, bol_flg, str_flg, ft_sz, clr)

% PRINT_CE_SPLIT(d_align, d_rdat, [seqpos], [mutpos], [xsel], [sequence], [offset], ...
%               [number_flag], [boolean_flag], [string_flag], [font_size], [color_code])
%
% Prints mutate-and-map electrophoregram in H x W splitted figures, with X-axis
%    labelled with mutants' name, Y-axis labelled with each nucleotides' position.
%
% Arguments
% =========
%   d_align         Required    Provides the electrophoregram matrix. Format in double
%                                   array.
%   d_rdat          Required    Provides d_rdat, at least for mutants names, e.g. G145C. 
%                                   Format in RDAT file. Provides SEQPOS, MUTPOS, XSEL, 
%                                   SEQUENCE, OFFSET if not included in arguments.                                   
%   [seqpos]        Optional    Provides bands annotation array. Format in string. If 
%                                   not provided, will be read out from D_RDAT.SEQPOS.  
%                                   Should be strictly monotonic decreasing, increasing  
%                                   array will be flipped automatically.
%   [mutpos]        Optional    Provides mutant position array. Format in double array. 
%                                   If not provided, will be read out from D_RDAT.MUTPOS. 
%                                   If not provided by D_RDAT (which is the case for RDAT
%                                   version 0.33+), will be generated from 
%                                   D_RDAT.DATA_ANNOTATIONS.
%   [xsel]          Optional    Provides annotation Y position. Format in double array.  
%                                   If not provided, willbe read out from D_RDAT.XSEL. 
%                                   Should be strict monotonicly increasing, decreasing  
%                                   array will be flipped automatically.
%   [sequence]      Optional    Provides sequence. Format in string. If not provided, 
%                                   will be read out from D_RDAT.SEQUENCE.
%   [offset]        Optional    Provides offset. Format in double. If not provided, will  
%                                   be read out from D_RDAT.OFFSET.
%   [number_flag]   Optional    Provides the layout format that will split into H x W 
%   [H  W  SC                       pages with vertical blank edge SP, using SC for 
%    SP LS SO                       intensity scaling, making lines every LS lanes. 
%    UO LO                          D_ALIGN will be auto-trimmed by UO and LO if asked,
%    YO XO]                         and axis titles will be placed with offset YO and
%                                   XO. Format in double array, default [0 0 0 50 10 ...
%                                   -1 75 100 0.01 0.19].
%                               'H' (height) denotes the number of pages vertically,
%                                   0 means use auto_size if AS in BOOLEAN_FLAG;
%                               'W' (width) denotes the number of pages horizontally,
%                                   0 means use auto_size if AS in BOOLEAN_FLAG;
%                               'SC' (scale_factor) denotes the scaling factor of image,
%                                   default is calculated to display optimally according 
%                                   to sc = 22.5 / mean(mean(d_align));
%                               'SP' (spacer) denotes the spacer width on both top and 
%                                   bottom of each sub-figure;
%                               'LS' (line_spacer) denotes the interval of blue thick 
%                                   lines;
%                               'SO' (line_offset) denotes the offset of blue thick 
%                                   lines. Default -1 means starts with MUTPOS 
%                                   divisible by 10.
%                               'UO' (upper_bound) denotes the excess D_ALIGN image 
%                                   shown above the first SEQPOS label. 
%                               'LO' (lower_bound) denotes the excess D_ALIGN image 
%                                   shown under the last SEQPOS label. 
%                               'YO' (y_title_offset) denotes the distance of y-axis
%                                   title to the axis.
%                               'XO' (x_title_offset) denotes the distance of x-axis
%                                   title to the axis.
%   [boolean_flag]  Optional    Provides the layout format that whether auto-trim 
%   [AT AS AL                       vertical boundaries, whether auto decide overall size, 
%    LN PR SQ                       whether auto decide text length for titles, whether draw 
%    TL VR PN]                      lines, whether print to file, whether squared page, 
%                                   whether add title, version label and page number. Format 
%                                   in double array, default [1 1 1 1 1 0 1 1 1].
%                               'AT' (is_auto_trim) denotes whether to auto trim top and 
%                                   bottom of D_ALIGN for optimal display. Excess image 
%                                   size is denoted by UO and LO in NUMBER_FLAG.
%                               'AS' (is_auto_size) denotes whether to auto decide number 
%                                   of pages based on D_ALIGN and SEQUENCE. This will 
%                                   override H and W in NUMBER_FLAG.
%                               'AL' (is_auto_length) denotes whether auto-decide font 
%                                   size of titles to fit in page margins. This will
%                                   override T1, T2, and T3 in FONT_SIZE.
%                               'LN' (is_line) denotes whether to make lines every LS of 
%                                   NUMBER_FLAG lanes or not;
%                               'PR' (if_print) denotes whether to print to .eps files, 
%                                   prints will be saved in folder DN of STRING_FLAG;
%                               'SQ' (is_square) denotes whether each page is squared;
%                               'TL' (is_title) denotes whether title is added to figure;
%                               'VR' (is_version) denotes whether version label is 
%                                   added to bottom right corner.
%                               'PN' (is_page_number) denotes whether page number is 
%                                   added to the corners of each figure;
%                               1 equals TRUE; 0 equals FALSE.
%   [string_flag]   Optional    Provides the string input for output .eps file name, 
%   [FN DN MF                       folder name, modifier label, authorship and date on
%    AU DT VR]                      printout. Format in string cell, default {'', ...
%                                   'print_CE_split_output', {modifier}, '', ...
%                                   'mmm yyyy', ''}.
%                               'FN' (file_name) denotes file name for print files.  
%                                   Numbers, underscore, and '.eps' extension will 
%                                   automatically append.
%                               'DN' (dir_name) denotes folder name for all files.
%                               'MF' (modifier_name) denotes modifier label. Default 
%                                   will use modifier entry from
%                                   D_RDAT.DATA_ANNOTATIONS. Useful when specifying 
%                                   SHAPE reagent with concentration used.
%                               'AU' (author_name) denotes authorship series. 'by' and 
%                                   '@' will automatically append.
%                               'DT' (date_string) denotes date string appearing on top 
%                                   right corner. Default will use current date.
%                               'VR' (hitrace_ver) denotes current HiTRACE subversion.
%   [font_size]     Optional    Provides font size values for figures. Format in double
%   [T1 T2 T3                       array, default [25 15 15 9 20 8 20 10].
%    VE YL YT                   'T1' font size of title (name) on page (1, 1);
%       XL XT]                  'T2' font size of title (conditions) on page (1, 2);
%                               'T3' font size of title (date) on page (1, 3);
%                               'VE' font size of version label on last page;
%                               'YL' font size of y-axis title (Sequence Position);
%                               'YT' font size of y-axis tick label;
%                               'XL' font size of x-axis title (Mutation Position);
%                               'XT' font size of x-axis tick label;
%   [color_code]    Optional    Provides color codes for figures. Format in string cell,
%   [T1 T2 T3                       default {'k', 'r', 'b', 'k', 'g', 'k', 'g', ...
%    VE YL YT                       'k', 'r', 'b', 'k'}.
%       XL XT                   'T1' font color of title (name) on page (1, 1);
%    L1 L2 L3]                  'T2' font color of title (conditions) on page (1, 2);
%                               'T3' font color of title (date) on page (1, 3);
%                               'YL' font color of y-axis title (Sequence Position);
%                               'VE' font color of version label on last page;
%                               'YT' font color of y-axis tick label;
%                               'XL' font color of x-axis title (Mutation Position);
%                               'XT' font color of x-axis tick label;
%                               'L1' line color of borders on each page;
%                               'L2' line color of all LS of NUMBER_FLAG thick lines;
%                               'L3' line color of border on each lane.
%
%
% e.g. PRINT_CE_SPLIT(d_align, d_rdat);
%      PRINT_CE_SPLIT(d_align, d_rdat, [], [], [], [], [], [], [1 1 1 0 1 1], ...
%                       {'', '', '1M7', 'T47', 'Apr 2013'}, [], {'k', 'b', 'r'});
%
%
% Notes
% =====
% Open all .eps files in single Preview window and print usin a color printer.
%   Deselect auto-rotate and auto-scale in the print dialog, make sure scale 
%   is 100%. Spaces on each border are included for easy splicing. Cut at red 
%   lines on each page. Assemble rows first, then put together the entire image.
% 
% ALL current opened figures will be LOST! Save before run.
%
% by T47, Apr 2013 - May 2013.
%

Script_VER = '2.1';

%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparation before splitting
% read offset, sequence, seqpos, mutpos, xsel from d_rdat if not exist
fprintf(['print_CE_split version ', Script_VER, '.\n']);
fprintf(['RDAT format version ', d_rdat.version,'.\n\n']);
fprintf(['d_align (',num2str(size(d_align,1)),' x ',num2str(size(d_align,2)),') read in.\n']);
d_align_size_H_org = size(d_align, 1);
d_align_size_W_org = size(d_align, 2);
fprintf(['d_rdat (', d_rdat.name, ') read in.\n']);

if ~exist('seqpos','var') || isempty(seqpos)
    seqpos = d_rdat.seqpos;
    fprintf(['seqpos (1 x ',num2str(length(seqpos)),') fetched from d_rdat.\n']); 
end;
if ~exist('mutpos','var') || isempty(mutpos)
    if exist('d.mutpos','var');
        mutpos = d_rdat.mutpos;
    else
        mutpos = [];
    end;
    fprintf(['mutpos (1 x ',num2str(length(mutpos)),') fetched from d_rdat.\n']); 
end;

% generate mutpos if not supplied by either argument or d_rdat
if isempty('mutpos') || length(mutpos) ~= size(d_align, 2)
    mutpos = generate_mutpos_from_rdat(d_rdat.data_annotations); 
end;

if ~exist('xsel','var') || isempty(xsel)
    xsel = d_rdat.xsel; 
    fprintf(['xsel (1 x ',num2str(length(xsel)),') fetched from d_rdat.\n']); 
end;
if ~exist('sequence','var') || isempty(sequence)
    sequence = d_rdat.sequence;
    fprintf(['sequence (1 x ',num2str(length(sequence)),') fetched from d_rdat.\n']); 
end;
if ~exist('offset','var') || isempty(offset)
    offset = d_rdat.offset; 
    fprintf(['offset (',num2str(offset),') fetched from d_rdat.\n']); 
end;
offset = round(offset);
fprintf('\n'); fprintf(['sequence  ',d_rdat.sequence,'\n']);
fprintf(['structure ',d_rdat.structure,'\n\n']);


% set all auxillary parameters
num_flg_org = [0 0 0 50 10 -1 75 100 0.01 0.19];
if ~exist('num_flg','var') || isempty(num_flg) ; 
    num_flg = num_flg_org;
else
    if length(num_flg) < length(num_flg_org);
        num_flg_org(1:length(num_flg)) = num_flg;
        num_flg = num_flg_org;
    end;
end;
num_flg([1:2, 4:7]) = round(num_flg([1:2, 4:7]));
num_flg(1) = max([num_flg(1), 0]);
num_flg(2) = max([num_flg(2), 0]);
num_flg(3) = max([num_flg(3), 0]);
if num_flg(4) < 1; num_flg(4) = 50; end;
if num_flg(5) < 1; num_flg(5) = 10; end;
num_flg(6) = max([num_flg(6), -1]);
page_num_H = num_flg(1); page_num_W = num_flg(2); scale_fc = num_flg(3);
num_sp = num_flg(4); num_line_sp = num_flg(5); num_line_offset = num_flg(6);
num_up_offset = num_flg(7); num_low_offset = num_flg(8);
num_y_offset = num_flg(9); num_x_offset = num_flg(10);

bol_flg_org = [1 1 1 1 1 0 1 1 1];
if ~exist('bol_flg','var') || isempty(bol_flg) ; 
    bol_flg = bol_flg_org;
else
    if length(bol_flg) < length(bol_flg_org);
        bol_flg_org(1:length(bol_flg)) = bol_flg;
        bol_flg = bol_flg_org;
    end;
end;
for i = length(bol_flg)
    bol_flg(i) = is_valid_boolean(bol_flg(i));
end;
is_auto_trim = bol_flg(1); is_auto_size = bol_flg(2); is_auto_length = bol_flg(3);
is_line = bol_flg(4); is_print = bol_flg(5); is_square = bol_flg(6);
is_title = bol_flg(7); is_version = bol_flg(8); is_page_no = bol_flg(9); 

str_flg_org = {'', 'print_CE_split_output', '', '', datestr(date, 'mmm yyyy'), ''};
if ~exist('str_flg','var') || isempty(str_flg)
    str_flg = str_flg_org;
else
    if length(str_flg) < length(str_flg_org);
        str_flg_org(1:length(str_flg)) = str_flg;
        str_flg = str_flg_org;
    end;
end;
file_name = str_flg{1}; dir_name = str_flg{2}; mdfr_str = str_flg{3};
author_str = str_flg{4}; date_str = [' @ ' str_flg{5}]; ver_hitrace = str_flg{6};
if ~isempty(file_name); file_name = [file_name '_']; end;
if isempty(dir_name); dir_name = 'print_CE_split_output'; end;
if ~isempty(author_str); author_str = [' by ' author_str]; end;
if isempty(date_str); date_str = datestr(date, 'mmm yyyy'); end;
if isempty(ver_hitrace); ver_hitrace = 'N/A'; end;

ft_sz_org = [25 15 15 9 20 8 20 10];
if ~exist('ft_sz','var') || isempty(ft_sz);
    ft_sz = ft_sz_org;
else
    if length(ft_sz) < length(ft_sz_org);
        ft_sz_org(1:length(ft_sz)) = ft_sz;
        ft_sz = ft_sz_org;
    end;
end;
ft_sz_title_1 = ft_sz(1); ft_sz_title_2 = ft_sz(2); ft_sz_title_3 = ft_sz(3);
ft_sz_ver = ft_sz(4);
ft_sz_y_title = ft_sz(5); ft_sz_y_tick = ft_sz(6);
ft_sz_x_title = ft_sz(7); ft_sz_x_tick = ft_sz(8);

clr_org = {'k', 'r', 'b', 'k', 'g', 'k', 'g', 'k', 'r', 'b', 'k'};
if ~exist('clr','var') || isempty(clr);
    clr = clr_org;
else
    if length(clr) < length(clr_org);
        clr_org(1:length(clr)) = clr;
        clr = clr_org;
    end;
end;
color_title_1 = clr{1}; color_title_2 = clr{2}; color_title_3 = clr{3};
color_ver = clr{4};
color_y_title = clr{5}; color_y_tick = clr{6};
color_x_title = clr{7}; color_x_tick = clr{8};
color_line_1 = clr{9}; color_line_2 = clr{10}; color_line_3 = clr{11};


% auto_size if asked
is_auto_size_force = (page_num_H == 0) && (page_num_W == 0);
[page_num_H, page_num_W] = auto_size(page_num_H, seqpos, page_num_W, mutpos, is_auto_size);

% print out summary
fprintf(['Divide into ', num2str(page_num_H), ' x ', num2str(page_num_W), ' pages, auto-size (', ...
    num2yn(is_auto_size || is_auto_size_force), ').\n']);
fprintf(['Spacer ', num2str(num_sp), ', make lines (', num2yn(is_line), ') every ', ...
    num2str(num_line_sp), ' lanes with offset (', num2str(num_line_offset), '), squared page (', num2yn(is_square), ').\n']);
fprintf(['Show title (', num2yn(is_title), '), show page number (', num2yn(is_page_no), ') show version (', ...
    num2yn(is_version), '), print to file (', num2yn(is_print), ').\n']);

% auto_trim if asked
fprintf(['Auto trim (', num2yn(is_auto_trim), ')']);
[d_align, xsel] = auto_trim(d_align, xsel, num_up_offset, num_low_offset, is_auto_trim);
fprintf('.\n');
fprintf(['Auto title size (', num2yn(is_auto_length), '), Y-axis title offset (', num2str(num_y_offset), ...
    ') and X-axis title offset (', num2str(num_x_offset), ').\n\n']);

% calculate auto_scale
auto_scale = 22.5 / mean(mean(d_align));
if scale_fc == 0; scale_fc = auto_scale; end;
fprintf( 'auto_scale_factor = %f\n', auto_scale);
fprintf( 'scale_used = %f\n\n', scale_fc);


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% split d_align into h*w sub matrices
% expand d_align to divisible size by h and w
[d_align, h_length, w_length] = d_expand_divisible(d_align, page_num_H, page_num_W);
fprintf(['In each page, there are ',num2str(w_length),' lanes (X-axis), and ',...
    num2str(h_length),' of traces (Y-axis).\n']);

% pause point
fprintf('\n'); fprintf('Press any key to continue...\n');
pause;


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% d(i, j, k, l):
% (i, j) denotes the subplot coordinate
% (k, l) denotes the band coordinate within one subplot

% add empty lane on both left and right of each figure
% add spacer on both top and bottom of each figure
d = d_split_add_spacer(d_align, page_num_H, h_length, num_sp, page_num_W, w_length, 1);

% read in mutants names (X-axis), trim to 'G145C' format
names = cell(1,size(d_align, 2)); 
for i = 1:length(mutpos)
    names{i} = strrep(d_rdat.data_annotations{i}{1}, 'mutation:', '');
end;

% split mutants names array into w*w_length array
name = cell(page_num_W, (w_length + 2)); 
for i = 1:page_num_W
    name(i, 2:(w_length + 1)) = names((w_length * (i - 1) + 1):(w_length * i));    
end;


% flip seqpos and xsel to correct order
% xsel be increasing, seqpos be decreasing, mutpos be increasing
[xsel_num_flag, xsel_str_flag] = check_monotone(xsel);
xsel_str_flag = lower(xsel_str_flag);
fprintf('\n');
fprintf(['Input xsel (1 x ',num2str(length(xsel)),') is ', xsel_str_flag]);
if xsel_num_flag == -2;
    fprintf(', FLIPPED for use.\n');
    xsel = fliplr(xsel);
elseif xsel_num_flag == 2;
    fprintf(', unchanged for use.\n');
else
    fprintf(', please check.\n');
    fprintf('** Strict monotonicity required! **\n');
end;

[seqpos_num_flag, seqpos_str_flag] = check_monotone(seqpos);
seqpos_str_flag = lower(seqpos_str_flag);
fprintf(['Input seqpos (1 x ',num2str(length(seqpos)),') is ', seqpos_str_flag]);
if seqpos_num_flag == 2;
    fprintf(', FLIPPED for use.\n');
    seqpos = fliplr(seqpos);
elseif seqpos_num_flag == -2;
    fprintf(', unchanged for use.\n');
else
    fprintf(', please check.\n');
    fprintf('** Strict monotonicity required! **\n');
end;

% read in band names (Y-axis)
% split band position array into h*h_length array
[bandpos, band] = xsel_label_split(xsel, seqpos, sequence, offset, page_num_H, h_length, num_sp);

% extend yticks on last column of figures if more than one blank lane there
extra_blank_lanes = w_length * page_num_W - d_align_size_W_org;
if extra_blank_lanes > 0;
    for i = 1:page_num_H
        xsel_pos_sub = []; ct = 1;
        for j = 1:size(band, 2)
            if ~isempty(band{i, j ,2});
                xsel_pos_sub(ct) = band{i, j, 2};
                ct = ct + 1;
            end;
        end;
        for j = 1:length(xsel_pos_sub)
            d(page_num_W, i, (round(xsel_pos_sub(j)) + [-1 0]), ...
                (end - extra_blank_lanes):end) = 22.5;
        end;
    end;
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot each figure
close all;

% make a folder to store print files
if is_print == 1; mkdir(dir_name); end;

for i = 1:page_num_W
    for j = 1:page_num_H
        
        % extract y-axis information from cell
        d_temp = []; 
        d_temp(1:size(d, 3), 1:size(d, 4)) = d(i, j, :, :);
        [xsel_pos_sub, xsel_txt_sub] = xsel_label_extract(band, j, 2);
        
        % add figure number to the top left and bottom right corner
        fig_name = ['[', num2str(j), ', ', num2str(i), '] of [',...
            num2str(page_num_H), ' ,', num2str(page_num_W), ']'];
        if is_page_no == 1;
            xsel_pos_sub(1) = round(num_sp / 2);
            xsel_txt_sub{1} = strrep(strcat(strtok(fig_name,']'), ']'), ' ', '');
            name{i, w_length + 2} = fig_name;
        end;
        
        fig_num = (i - 1) * page_num_H + j;
        h = figure(fig_num);
        
        % adjustment for auto_title_length
        fig_w = 600; title_size_fc = 0.775; title_h = 0.45;
        set(gcf, 'PaperOrientation', 'Portrait', 'PaperPositionMode', 'Manual', ...
            'PaperSize', [8.5 11], 'Color', 'White');
        if is_square == 1; 
            set(gcf, 'PaperPosition', [0 2 8.5 8.5]);
            fig_h = 600; title_offset = 0.025;
        else 
            set(gcf, 'PaperPosition', [0 1 8.5 10]);
            fig_h = 800; title_offset = 0;
        end;
        set(h, 'Position', [(fig_num - 1) * 20, 0, fig_w, fig_h]);
        image(d_temp * scale_fc);
        colormap(1 - gray());
        
        % make lines
        % blue thick lines per LS
        if is_line == 1
            
            % find mutpos divisible by 10
            if num_line_offset == -1;
                num_line_offset = mod(mutpos(2), num_line_sp);
            end;
            
            make_lines(0:1:w_length, color_line_3, 1);
            line_start = mod(num_line_offset + num_line_sp + 1 - mod(w_length * (i - 1), num_line_sp), num_line_sp);
            make_lines(line_start:num_line_sp:(w_length + 1), color_line_2, 2);
        end;
        
        % red border lines
        make_lines_horizontal(num_sp, color_line_1, 1);
        make_lines_horizontal(size(d_temp,1) - num_sp, color_line_1, 1);
        make_lines(1, color_line_1, 1);
        make_lines(size(d_temp, 2)-1, color_line_1, 1);
        
        % axis labeling
        % X-axis-tick, mutant names
        set(gca, 'FontSize', ft_sz_x_tick);
        set(gca, 'XTick', 1:size(name, 2), 'XTickLabel', char(name{i, :}), 'XColor', color_x_tick);
        xticklabel_rotate(); 
        
        % Y-axis-tick, band positions
        set(gca, 'FontSize', ft_sz_y_tick, 'YColor', color_y_tick);
        if i == page_num_W;
            set(gca, 'YTick', xsel_pos_sub(:), 'YTickLabel', char(xsel_txt_sub(:)), 'YAxisLoc', 'Right');
        else
            set(gca, 'YTick', xsel_pos_sub(:), 'YTickLabel', char(xsel_txt_sub(:)), 'YAxisLoc', 'Left');
        end;
        
        % X-axis-caption
        xlabel('Mutation Position', 'Color', color_x_title); 
        if j == page_num_H;
            xlabh = get(gca, 'XLabel'); set(xlabh, 'FontSize', ft_sz_x_title);
            set(xlabh, 'Position', get(xlabh, 'Position') + [0 (num_x_offset + title_offset) 0]);
        end;
        
        % Y-axis-caption
        ylabel('Sequence Position', 'Color', color_y_title);
        ylabh = get(gca, 'YLabel'); set(ylabh, 'FontSize', ft_sz_y_title);
        set(ylabh, 'Position', get(ylabh, 'Position') + [num_y_offset 0 0]);
        
        % add title if asked
        if is_title == 1;
            
            % modifier info from rdat if not specified
            if isempty(mdfr_str); mdfr_str = annotation_type_finder(d_rdat.annotations, 'modifier'); end;
            
            % title of left-top corner, for construct name
            if (i == 1 && j == 1);
                
                title([' ' d_rdat.name ' '], 'HorizontalAlignment', 'Left', 'VerticalAlignment', 'Bottom',...
                    'FontWeight', 'Bold', 'FontSize', ft_sz_title_1, 'FontName', 'Courier', 'Color', color_title_1);
                tit = get(gca, 'Title'); pos = get(tit, 'Position');
                set(tit, 'Position', [0 (pos(2) + title_offset) pos(3)]);
                if is_auto_length == 1; optimal_font_size(tit, min((get(gcf,'PaperSize'))) * title_size_fc, title_h); end;

            % title of second top page, for experiment details    
            elseif (j == 1 && i == 2);
                
                % two-line format, second line for all chemicals
                note{1} = [annotation_type_finder(d_rdat.annotations, 'experimentType'), '  ', ...
                    annotation_type_finder(d_rdat.annotations, 'modifier'), '  ',...
                    annotation_type_finder(d_rdat.annotations, 'temperature'), '  '];
                note{2} = [annotation_type_finder(d_rdat.annotations, 'chemical'), '  '];
                
                % include right-top corner text if only Nx2 pages
                if page_num_W == 2; note{3} = [mdfr_str author_str date_str '  ']; end;
                
                title(note, 'HorizontalAlignment', 'Right', 'VerticalAlignment', 'Bottom',...
                    'FontSize', ft_sz_title_2, 'FontName', 'Courier', 'Color', color_title_2);
                tit = get(gca, 'Title'); pos = get(tit, 'Position');
                set (tit, 'Position', [1 (pos(2) + title_offset) pos(3)]);
                if is_auto_length == 1; optimal_font_size(tit, min((get(gcf,'PaperSize'))) * title_size_fc, title_h); end;

            % title of right-top corner, for experiment date and authorship
            elseif (j == 1 && i == page_num_W)
                comment = [mdfr_str author_str date_str '  '];
                
                title(comment, 'HorizontalAlignment', 'Right', 'VerticalAlignment', 'Bottom',...
                    'FontWeight', 'Bold', 'FontSize', ft_sz_title_3, 'FontName', 'Courier', 'Color', color_title_3);
                tit = get(gca, 'Title'); pos = get(tit, 'Position');
                set (tit, 'Position', [1 (pos(2) + title_offset) pos(3)]);
                if is_auto_length == 1; optimal_font_size(tit, min((get(gcf,'PaperSize'))) * title_size_fc, title_h); end;
                
            end;
        end;
        
        % version label of right-bottom corner
        if (j == page_num_H && i == page_num_W && is_version)
            ver_str = ['HiTRACE rev. ', num2str(ver_hitrace), ' / RDAT ver. ', d_rdat.version, ...
                ' / print\_CE\_split ver. ', Script_VER, '  '];
            
            ylim = get(gca, 'YLim'); xlim = get(gca, 'XLim');
            if is_square == 1; 
                ver_offset_all = [220 110 90 70 50 40];
                ver_offset = ver_offset_all(page_num_H); 
            else
                ver_offset_all = [170 90 70 50 40 35];
                ver_offset = ver_offset_all(page_num_H); 
            end;
            tit = text(xlim(2), ylim(2) + ver_offset, ver_str);
            set(tit, 'HorizontalAlignment', 'Right', 'VerticalALignment', 'Top', 'FontSize', ft_sz_ver, 'Color', color_ver);
        end;
            
        % print to file if asked
        if is_print == 1; print_save_figure(h, [file_name, num2str(fig_num)], dir_name, 0); end;
    end;
end;

if is_print == 1; fprintf([num2str(i*j),' pages printed to folder "', dir_name,'".\n']); end;

