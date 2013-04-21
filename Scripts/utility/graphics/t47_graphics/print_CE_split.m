function print_CE_split(d_align, d_rdat, seqpos, mutpos, xsel, sequence, offset, space)

% PRINT_CE_SPLIT(d_align, d_rdat, [seqpos], [mutpos], [xsel], [sequence], [offset],...
%               [h w sp ln ln_sp tl pr sq sc])
%
% Prints the mutate-and-map electrophoregram in h*w splitted figures, with X-axis
%    labelled with mutants' name, Y-axis labelled with each nucleotides' position.
%
% Arguments
% =========
%   d_align         Required    Provides the electrophoregram matrix.
%   d_rdat          Required    Provides d_rdat, at least for mutants names, e.g. G145C. 
%                                   Provides seqpos, mutpos, xsel, sequence, offset if 
%                                   not included in arguments.                                   
%   [seqpos]        Optional    Provides bands annotation array. If not provided, will
%                                   be read out from d_rdat.seqpos. Should be strict 
%                                   monotonic decreasing, increasing array will be 
%                                   flipped automatically.
%   [mutpos]        Optional    Provides mutant position array. If not provided, will 
%                                   be read out from d_rdat.mutpos. If not provided by 
%                                   d_rdat, will be generated from 
%                                   d_rdat.data_annotations.
%   [xsel]          Optional    Provides annotation Y position. If not provided, will 
%                                   be read out from d_rdat.xsel. Should be strict 
%                                   monotonic increasing, decreasing array will be 
%                                   flipped automatically.
%   [sequence]      Optional    Provides sequence. If not provided, will be read out
%                                   from d_rdat.sequence.
%   [offset]        Optional    Provides offset. If not provided, will be read out from 
%                                   d_rdat.offset.
%   [h w sp ln      Optional    Provides the layout format that will split into h*w 
%    ln_sp tl                       pages with vertical blank edge sp. 
%    pr sq sc]                  h (height) denotes the number of pages vertically;
%                               w (width) denotes the number of pages horizontally;
%                               sp (spacer) denotes the spacer width on both top and 
%                                   bottom of each sub-figure;
%                               ln (if_line) denotes whether to make lines every ln_sp  
%                                   lanes or not;
%                               ln_sp (line_spacer) denotes the interval of blue thick 
%                                   lines;
%                               tl (if_title) denotes whether title is added to figure;
%                               pr (if_print) denotes whether to print to .eps files, 
%                                   prints will be saved in folder "print_CE_split";
%                               sq (if_square) denotes whether each page is squared;
%                                   1 equals true; 0 equals false.
%                               sc (scale_factor) denotes the scaling factor of image,
%                                   default is calculated to display optimally 
%                                   according to sc = 22.5 / mean(mean(d_align));
%                               Default [3 3 50 1 10 1 1 1 auto_scale]: print in 3*3, 
%                                   spacer 50, with lines every 10 lanes, with title,
%                                    print to file, squared page, image scaled optimally.
%                                    
%
% e.g. PRINT_CE_SPLIT(d_align, d_rdat, [], [], [], [], [], []);
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
% by T47, Feb 2013
%


% read offset, sequence, seqpos, mutpos, xsel from d_rdat if not exist
fprintf(['d_align (',num2str(size(d_align,1)),'x',num2str(size(d_align,2)),') read in.\n']);
fprintf(['d_rdat (', d_rdat.name, ') created by ver ',d_rdat.version,' read in.\n']);
fprintf(['sequence  ',d_rdat.sequence,'\n']);
fprintf(['structure ',d_rdat.structure,'\n']);
fprintf('\n');

if isempty(offset)
    offset = d_rdat.offset; 
    fprintf(['offset (',num2str(offset),') fetched from d_rdat.\n']); 
end;
if isempty(sequence)
    sequence = d_rdat.sequence;
    fprintf(['sequence (1x',num2str(length(sequence)),') fetched from d_rdat.\n']); 
end;
if isempty(seqpos)
    seqpos = (d_rdat.seqpos);
    fprintf(['seqpos (1x',num2str(length(seqpos)),') fetched from d_rdat.\n']); 
end;
if isempty(mutpos)
    mutpos = d_rdat.mutpos; 
    fprintf(['mutpos (1x',num2str(length(mutpos)),') fetched from d_rdat.\n']); 
end;
if isempty(xsel)
    xsel = d_rdat.xsel; 
    fprintf(['xsel (1x',num2str(length(xsel)),') fetched from d_rdat.\n']); 
end;

% set [h w sp ln ln_sp tl pr sq sc] to default if not specified
auto_scale = 22.5 / mean(mean(d_align));
space_temp = [3 3 50 1 10 1 1 1 auto_scale];
if isempty(space) || length(space) < 9; 
    space_temp(1:length(space)) = space;
    space = space_temp;
end;
space(1:8) = round(space(1:8)); offset = round(offset);
if space(1) < 2; space(1) = 1; end;
if space(2) < 2; space(2) = 1; end;
if space(3) < 1; space(3) = 0; end;
if space(5) < 1; space(5) = 1; end;
if space(9) < 0; space(9) = 0; end;

print_summary(space);
fprintf( 'auto_scale_factor = %f\n', auto_scale);
fprintf( 'scale_used = %f\n', space(9));

% generate mutpos if not supplied by either argument or d_rdat
if isempty(mutpos) || length(mutpos) ~= size(d_align, 2)
    mutpos = generate_mutpos_from_rdat(d_rdat.data_annotations); 
end;

% split d_align into h*w sub matrices
% expand d_align to divisible size by h and w
h_length = floor(size(d_align, 1)/space(1));
if h_length * space(1) ~= size(d_align, 1); h_length = h_length + 1; end;
d_align((size(d_align, 1) + 1):(h_length * space(1)), :) = 0;

w_length = floor(size(d_align, 2)/space(2));
if w_length * space(2) ~= size(d_align, 2); w_length = w_length + 1; end;
d_align(:, (size(d_align, 2) + 1):(w_length * space(2))) = 0;

fprintf(['In each page, there are ',num2str(w_length),' lanes (x-axis), and ',...
    num2str(h_length),' of traces (y-axis).\n']);

% d(i, j, k, l):
% (i, j) denotes the subplot coordinate
% (k, l) denotes the band coordinate within one subplot

% add empty lane on both left and right of each figure
% add spacer on both top and bottom of each figure
d = zeros([space(2) space(1) (h_length + 2 * space(3)) (w_length + 2)]); 
for i = 1:space(2)
    for j = 1:space(1)
        d(i, j, ((space(3) + 1):(space(3) + h_length)), 2:(size(d, 4) - 1)) = d_align((h_length * (j - 1) + 1):(h_length * j),...
            (w_length * (i - 1) + 1):(w_length * i));
    end;
end;


% read in mutants names (X-axis), trim to 'G145C' format
names = cell(1,size(d_align, 2)); 
for i = 1:length(mutpos)
    names(i) = strrep(d_rdat.data_annotations{i}, 'mutation:', '');
end;

% split mutants names array into w*w_length array
name = cell(space(2), (w_length + 2)); 
for i = 1:space(2)
    name(i, 2:(w_length + 1)) = names((w_length * (i - 1) + 1):(w_length * i));    
end;

% flip seqpos and xsel to correct order
% xsel be increasing, seqpos be decreasing, mutpos be increasing
[xsel_num_flag, xsel_str_flag] = check_monotone(xsel);
fprintf('\n');
fprintf(['Input xsel (1x',num2str(length(xsel)),') is ', xsel_str_flag]);
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
fprintf(['Input seqpos (1x',num2str(length(seqpos)),') is ', seqpos_str_flag]);
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
band = cell(space(1), size(bandpos, 1), 2);
for i = 1:space(1)
    ymin = h_length * (i - 1) + 1; ymax = h_length * i;
    for j = 1:size(bandpos, 1)
        if bandpos{j, 2}>=ymin && bandpos{j, 2}<=ymax;
            band{i, j, 1} = bandpos{j, 1}; 
            band{i, j, 2} = bandpos{j, 2} + space(3) - (ymin - 1);
        end;
    end;
end;

% plot each figure
close all;

% make a folder to store print files
if space(7) == 1;
    dir_name = 'print_CE_split';
    mkdir(dir_name);
end;
for i = 1:space(2)
    for j = 1:space(1)
        
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
            num2str(space(1)), ' ,', num2str(space(2)), ']'];
        
        fig_num = (i - 1) * space(1) + j;
        h = figure(fig_num);
        fig_h = 800; if space(8) == 1; fig_w = 800; else; fig_w = 600; end;
        set(h,'Position',[(fig_num - 1) * 20, 0, fig_w, fig_h]);
        set(gcf, 'PaperOrientation', 'landscape', 'PaperPositionMode', 'auto', 'color', 'white');
        image(d_temp * space(9));
        colormap(1 - gray());
        
        % make lines
        if space(4) == 1
            make_lines(0:1:w_length, 'k', 1);
            line_start = space(5) + 1 - mod(w_length * (i - 1), space(5));
            make_lines(line_start:space(5):(w_length + 1), 'b', 2);
        end;
        make_lines_horizontal(space(3), 'r', 1);
        make_lines_horizontal(size(d_temp,1)-space(3), 'r', 1);
        make_lines(1, 'r', 1);
        make_lines(size(d_temp,2)-1, 'r', 1);
        
        % axis labeling
        % x-axis-tick, mutant names
        set(gca, 'xtick', 1:size(name, 2), 'xticklabel', char(name{i, :}));
        xticklabel_rotate; 
        
        % y-axis-tick, band positions
        set(gca, 'Fontsize', 8);
        if i == space(2);
            set(gca, 'ytick', band_temp(:), 'yticklabel', char(label_temp(:)), 'yaxisloc', 'right');
        else
            set(gca, 'ytick', band_temp(:), 'yticklabel', char(label_temp(:)), 'yaxisloc', 'left');
        end;
        
        % x-axis-caption
        xlabel('Mutation Position', 'Color', [0 1 0]); 
        if j == space(1);
            xlabh = get(gca, 'XLabel'); set(xlabh, 'fontsize', 20);
            set(xlabh, 'Position', get(xlabh, 'Position') + [0 .19 0])
        end;
        
        % y-axis-caption
        ylabel('Sequence Position', 'Color', [0 1 0]);
        ylabh = get(gca, 'YLabel'); set(ylabh, 'Fontsize', 20);
        set(ylabh, 'Position', get(ylabh, 'Position') + [.01 0 0])
        
        % add title if asked
        if space(6) == 1
            
            % title of left-top corner, for construct name
            if i == 1 && j == 1
                title(d_rdat.name, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom',...
                    'FontWeight', 'bold', 'FontSize', 25, 'FontName', 'Courier', 'Color', [0 0 0]);
                tit = get(gca, 'Title'); pos = get(tit, 'Position');
                pos(1) = 0; set (tit, 'Position', pos);
                
            % title of second top page, for experiment details    
            elseif (j == 1 && i == 2)
                note = [strrep(d_rdat.annotations(1), 'experimentType:', ''), strrep(d_rdat.annotations(5), 'modifier:', ' '), ...
                    strrep(d_rdat.annotations(3), 'chemical:', ' '), strrep(d_rdat.annotations(2), 'chemical:', ' '),...
                    strrep(d_rdat.annotations(4), 'temperature:', ' '), ' `'];
                title(note, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom',...
                    'FontSize', 15, 'FontName', 'Courier', 'Color', [1 0 0]);
                tit = get(gca, 'Title'); pos = get(tit, 'Position');
                set (tit, 'Position', [1 pos(2) pos(3)]);
                
            % title of right-top corner, for experiment date
            elseif (j == 1 && i == space(2))
                title(datestr(date,'mmm yyyy'), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom',...
                    'FontWeight', 'bold', 'FontSize', 15, 'FontName', 'Courier', 'Color', [0 0 1]);
                tit = get(gca, 'Title'); pos = get(tit, 'Position');
                set (tit, 'Position', [1 pos(2) pos(3)]);
            end;
        end;
        
        % print to file if asked
        if space(7) == 1; print(h,'-depsc2',[dir_name, '/', num2str(fig_num),'.eps']); end;
    end;
end;
if space(7) == 1; fprintf([num2str(i*j),' pages printed to folder "', dir_name,'".\n']); end;


function print_summary(space)

fprintf('\n');
fprintf(['Divide into ', num2str(space(1)), 'x', num2str(space(2)), ' pages.\n']);
fprintf(['Spacer ', num2str(space(3)), ', make lines (', num2yn(space(4)),...
    ') every ', num2str(space(5)), ' lanes.\n']);
fprintf(['Show title (', num2yn(space(6)), '), print to file (', num2yn(space(7)),...
    '), squared page (', num2yn(space(8)), ').\n']);
fprintf('\n');




