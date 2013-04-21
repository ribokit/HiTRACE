function print_CE_split(d_align, d_rdat, seqpos, mutpos, xsel, sequence, offset, flg, ft_sz, clr, file_name)

% PRINT_CE_SPLIT(d_align, d_rdat, [seqpos], [mutpos], [xsel], [sequence], [offset], ...
%               [boolean_flag], [font_size], [color_code], [filename])
%
% Prints mutate-and-map electrophoregram in H x W splitted figures, with X-axis
%    labelled with mutants' name, Y-axis labelled with each nucleotides' position.
%
% Arguments
% =========
%   d_align         Required    Provides the electrophoregram matrix. Format in double
%                                   array.
%   d_rdat          Required    Provides d_rdat, at least for mutants names, e.g. G145C. 
%                                   Format in RDAT file. Provides seqpos, mutpos, xsel, 
%                                   sequence, offset if not included in arguments.                                   
%   [seqpos]        Optional    Provides bands annotation array. Format in string. If 
%                                   not provided, will be read out from d_rdat.seqpos.  
%                                   Should be strictly monotonic decreasing, increasing  
%                                   array will be flipped automatically.
%   [mutpos]        Optional    Provides mutant position array. Format in double array. 
%                                   If not provided, will be read out from d_rdat.mutpos. 
%                                   If not provided by d_rdat, will be generated from 
%                                   d_rdat.data_annotations.
%   [xsel]          Optional    Provides annotation Y position. Format in double array.  
%                                   If not provided, willbe read out from d_rdat.xsel. 
%                                   Should be strict monotonicly increasing, decreasing  
%                                   array will be flipped automatically.
%   [sequence]      Optional    Provides sequence. Format in string. If not provided, 
%                                   will be read out from d_rdat.sequence.
%   [offset]        Optional    Provides offset. Format in double. If not provided, will  
%                                   be read out from d_rdat.offset.
%   [boolean_flag]  Optional    Provides the layout format that will split into H x W
%   [H W SP LN                      pages with vertical blank edge sp. Format in double
%    LS TL PR                       array. Default [3 3 50 1 10 1 1 1 0 1]: print in
%    SQ SC TC]                      3 x 3, spacer 50, with lines every 10 lanes, with 
%                                   title, print to file, squared page, image scaled
%                                   optimally, auto-trimmed vertical boundaries.
%                               'H' (height) denotes the number of pages vertically;
%                               'W' (width) denotes the number of pages horizontally;
%                               'SP' (spacer) denotes the spacer width on both top and 
%                                   bottom of each sub-figure;
%                               'LN' (if_line) denotes whether to make lines every ln_sp  
%                                   lanes or not;
%                               'LS' (line_spacer) denotes the interval of blue thick 
%                                   lines;
%                               'TL' (if_title) denotes whether title is added to figure;
%                               'PR' (if_print) denotes whether to print to .eps files, 
%                                   prints will be saved in folder "print_CE_split";
%                               'SQ' (if_square) denotes whether each page is squared;
%                               'SC' (scale_factor) denotes the scaling factor of image,
%                                   default is calculated to display optimally according 
%                                   to sc = 22.5 / mean(mean(d_align));
%                               'TC' (if_auto_cut) denotes whether to auto trim top and 
%                                   bottom of d_align for optimal display. 
%                               1 equals TRUE; 0 equals FALSE.
%   [font_size]     Optional    Provides font size values for figures. Format in double
%                                   array. Default [].
%                               'T1'
%                               'T2'
%                               'T3'
%                               'YL'
%                               'YT'
%                               'XL'
%                               'XT'
%   [color_code]    Optional    Provides color codes for figures. Format in string cell.
%                                   Default {}.
%                               'T1'
%                               'T2'
%                               'T3'
%                               'YL'
%                               'YT'
%                               'XL'
%                               'XT'
%                               'L1'
%                               'L2'
%                               'L3'
%   [filename]      Optional    Provides file name for print files. Format in string. 
%                                   Numbers, underscore, and '.eps' extension will 
%                                   automatically append. Default is blank.
%
% e.g. PRINT_CE_SPLIT(d_align, d_rdat);
%      PRINT_CE_SPLIT(d_align, d_rdat, [], [], 
%                                   
%
% Notes
% =====
% Print all the figures on US Letter, with "landscape fit":
%    [left -0.25; top 0.25; width 11.50; height 8.50]
%    For bottom figures, use top 0 to print X-axis fully; for top figures, use top 0.50
%        to print title fully.
% Spaces on each border are included for easy splicing. Cut at red lines on each page.
%
% ALL current opened figures will be LOST! Save before run.
%
% by T47, Apr 2013
%


% read offset, sequence, seqpos, mutpos, xsel from d_rdat if not exist
fprintf(['d_align (',num2str(size(d_align,1)),' x ',num2str(size(d_align,2)),') read in.\n']);
fprintf(['d_rdat (', d_rdat.name, ') created by ver ',d_rdat.version,' read in.\n']);
if ~exist('seqpos','var') || isempty(seqpos)
    seqpos = d_rdat.seqpos;
    fprintf(['seqpos (1 x ',num2str(length(seqpos)),') fetched from d_rdat.\n']); 
end;
if ~exist('mutpos','var') || isempty(mutpos)
    mutpos = d_rdat.mutpos; 
    fprintf(['mutpos (1 x ',num2str(length(mutpos)),') fetched from d_rdat.\n']); 
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
fprintf('\n'); fprintf(['sequence  ',d_rdat.sequence,'\n']);
fprintf(['structure ',d_rdat.structure,'\n']); fprintf('\n');

% set [h w sp ln ls tl pr sq sc tc] to default if not specified
flg_temp = [3 3 50 1 10 1 1 1 0 1];
if ~exist('flg','var') || isempty(flg) ; 
    flg = flg_temp;
else;
    if length(flg) < 10;
        flg_temp(1:length(flg)) = flg;
        flg = flg_temp;
    end;
end;
flg([1:8 10]) = round(flg([1:8 10])); offset = round(offset);
if flg(1) < 2; flg(1) = 1; end;
if flg(2) < 2; flg(2) = 1; end;
if flg(3) < 1; flg(3) = 50; end;
if flg(5) < 1; flg(5) = 10; end;
if flg(9) < 0; flg(9) = 0; end;

if ~exist('filename','var'); file_name = ''; end;
if ~isempty(file_name); file_name = [file_name '_']; end;

% print out summary
fprintf(['Divide into ', num2str(flg(1)), ' x ', num2str(flg(2)), ' pages.\n']);
fprintf(['Spacer ', num2str(flg(3)), ', make lines (', num2yn(flg(4)),...
    ') every ', num2str(flg(5)), ' lanes.\n']);
fprintf(['Show title (', num2yn(flg(6)), '), print to file (', num2yn(flg(7)),...
    '), squared page (', num2yn(flg(8)), ').\n']);
fprintf('\n');

% auto trim if asked


% calculate auto_scale
auto_scale = 22.5 / mean(mean(d_align));
if flg(9) == 0; flg(9) = auto_scale; end;
fprintf( 'auto_scale_factor = %f\n', auto_scale);
fprintf( 'scale_used = %f\n', flg(9));

% generate mutpos if not supplied by either argument or d_rdat
if isempty('mutpos') || length(mutpos) ~= size(d_align, 2)
    mutpos = generate_mutpos_from_rdat(d_rdat.data_annotations); 
end;

% split d_align into h*w sub matrices
% expand d_align to divisible size by h and w
h_length = floor(size(d_align, 1)/flg(1));
if h_length * flg(1) ~= size(d_align, 1); h_length = h_length + 1; end;
d_align((size(d_align, 1) + 1):(h_length * flg(1)), :) = 0;

w_length = floor(size(d_align, 2)/flg(2));
if w_length * flg(2) ~= size(d_align, 2); w_length = w_length + 1; end;
d_align(:, (size(d_align, 2) + 1):(w_length * flg(2))) = 0;

fprintf(['In each page, there are ',num2str(w_length),' lanes (X-axis), and ',...
    num2str(h_length),' of traces (Y-axis).\n']);

% pause point
fprintf('\n'); fprintf('Press any key to continue...\n');
pause;


% d(i, j, k, l):
% (i, j) denotes the subplot coordinate
% (k, l) denotes the band coordinate within one subplot

% add empty lane on both left and right of each figure
% add flgr on both top and bottom of each figure
d = zeros([flg(2) flg(1) (h_length + 2 * flg(3)) (w_length + 2)]); 
for i = 1:flg(2)
    for j = 1:flg(1)
        d(i, j, ((flg(3) + 1):(flg(3) + h_length)), 2:(size(d, 4) - 1)) = d_align((h_length * (j - 1) + 1):(h_length * j),...
            (w_length * (i - 1) + 1):(w_length * i));
    end;
end;


% read in mutants names (X-axis), trim to 'G145C' format
names = cell(1,size(d_align, 2)); 
for i = 1:length(mutpos)
    names{i} = strrep(d_rdat.data_annotations{i}{1}, 'mutation:', '');
end;

% split mutants names array into w*w_length array
name = cell(flg(2), (w_length + 2)); 
for i = 1:flg(2)
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
bandpos = cell(length(mutpos), 2);
for i = 1:length(seqpos)
    bandpos{i, 1} = [sequence(seqpos(i) - offset), num2str(seqpos(i))];
    bandpos{i, 2} = xsel(i);
end;

% split band position array into h*h_length array
band = cell(flg(1), size(bandpos, 1), 2);
for i = 1:flg(1)
    ymin = h_length * (i - 1) + 1; ymax = h_length * i;
    for j = 1:size(bandpos, 1)
        if bandpos{j, 2}>=ymin && bandpos{j, 2}<=ymax;
            band{i, j, 1} = bandpos{j, 1}; 
            band{i, j, 2} = bandpos{j, 2} + flg(3) - (ymin - 1);
        end;
    end;
end;


% plot each figure
close all;

% make a folder to store print files
if flg(7) == 1;
    dir_name = 'print_CE_split';
    mkdir(dir_name);
end;
for i = 1:flg(2)
    for j = 1:flg(1)
        
        % extract axis information from cell
        d_temp = []; band_temp = []; label_temp = {}; ct = 1;
        d_temp(1:size(d, 3), 1:size(d, 4)) = d(i, j, :, :);
        for k = 1:size(band, 2)
            if ~isempty(band{j, k, 2}); 
                band_temp(ct) = band{j, k, 2};
                label_temp{ct} = band{j, k, 1};
                ct = ct + 1;
            end;
        end;
        %band_temp = band_temp(end:-1:1);
        %label_temp = label_temp(end:-1:1);
        
        % add figure number to the bottom right corner
        name{i, w_length + 2} = ['[', num2str(j), ', ', num2str(i), '] of [',...
            num2str(flg(1)), ' ,', num2str(flg(2)), ']'];
        
        fig_num = (i - 1) * flg(1) + j;
        h = figure(fig_num);
        fig_h = 800; if flg(8) == 1; fig_w = 800; else; fig_w = 600; end;
        set(h,'Position',[(fig_num - 1) * 20, 0, fig_w, fig_h]);
        set(gcf, 'PaperOrientation', 'landscape', 'PaperPositionMode', 'auto', 'color', 'white');
        image(d_temp * flg(9));
        colormap(1 - gray());
        
        % make lines
        if flg(4) == 1
            make_lines(0:1:w_length, 'k', 1);
            line_start = flg(5) + 1 - mod(w_length * (i - 1), flg(5));
            make_lines(line_start:flg(5):(w_length + 1), 'b', 2);
        end;
        make_lines_horizontal(flg(3), 'r', 1);
        make_lines_horizontal(size(d_temp,1)-flg(3), 'r', 1);
        make_lines(1, 'r', 1);
        make_lines(size(d_temp,2)-1, 'r', 1);
        
        % axis labeling
        % X-axis-tick, mutant names
        set(gca, 'xtick', 1:size(name, 2), 'xticklabel', char(name{i, :}));
        xticklabel_rotate; 
        
        % Y-axis-tick, band positions
        set(gca, 'Fontsize', 8);
        if i == flg(2);
            set(gca, 'ytick', band_temp(:), 'yticklabel', char(label_temp(:)), 'yaxisloc', 'right');
        else
            set(gca, 'ytick', band_temp(:), 'yticklabel', char(label_temp(:)), 'yaxisloc', 'left');
        end;
        
        % X-axis-caption
        xlabel('Mutation Position', 'Color', [0 1 0]); 
        if j == flg(1);
            xlabh = get(gca, 'XLabel'); set(xlabh, 'fontsize', 20);
            set(xlabh, 'Position', get(xlabh, 'Position') + [0 .19 0])
        end;
        
        % Y-axis-caption
        ylabel('Sequence Position', 'Color', [0 1 0]);
        ylabh = get(gca, 'YLabel'); set(ylabh, 'Fontsize', 20);
        set(ylabh, 'Position', get(ylabh, 'Position') + [.01 0 0])
        
        % add title if asked
        if flg(6) == 1
            
            % title of left-top corner, for construct name
            if i == 1 && j == 1
                title(d_rdat.name, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom',...
                    'FontWeight', 'bold', 'FontSize', 25, 'FontName', 'Courier', 'Color', [0 0 0]);
                tit = get(gca, 'Title'); pos = get(tit, 'Position');
                pos(1) = 0; set (tit, 'Position', pos);
                
            % title of second top page, for experiment details    
            elseif (j == 1 && i == 2)
                note = strcat('` ', strrep(d_rdat.annotations(1), 'experimentType:', ''), strrep(d_rdat.annotations(5), 'modifier:', ' '), ...
                    strrep(d_rdat.annotations(3), 'chemical:', ' '), strrep(d_rdat.annotations(2), 'chemical:', ' '),...
                    strrep(d_rdat.annotations(4), 'temperature:', ' '), ' `');
                title(note, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom',...
                    'FontSize', 15, 'FontName', 'Courier', 'Color', [1 0 0]);
                tit = get(gca, 'Title'); pos = get(tit, 'Position');
                set (tit, 'Position', [1 pos(2) pos(3)]);
                
            % title of right-top corner, for experiment date
            elseif (j == 1 && i == flg(2))
                title(datestr(date,'mmm yyyy'), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom',...
                    'FontWeight', 'bold', 'FontSize', 15, 'FontName', 'Courier', 'Color', [0 0 1]);
                tit = get(gca, 'Title'); pos = get(tit, 'Position');
                set (tit, 'Position', [1 pos(2) pos(3)]);
            end;
        end;
        
        % print to file if asked
        if flg(7) == 1; print(h,'-depsc2',[dir_name, '/', file_name, num2str(fig_num),'.eps']); end;
    end;
end;

if flg(7) == 1; fprintf([num2str(i*j),' pages printed to folder "', dir_name,'".\n']); end;

