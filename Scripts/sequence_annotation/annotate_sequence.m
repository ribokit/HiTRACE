function [xsel,seqpos,area_pred] = annotate_sequence(d_align, xsel, sequence_full, ...
    offset, data_types, first_RT_nucleotide, structure, font_size, JUST_PLOT)
% ANNOTATE_SEQUENCE - Tool for rapid manual assignment of bands in electropherograms.
%
%  [xsel,seqpos,area_pred] = annotate_sequence(d_align, xsel, sequence_full, offset, data_types, first_RT_nucleotide, structure, font_size, JUST_PLOT);
%
%
% Input:
%  d_align        = matrix of aligned electrophoretic traces.
%  xsel           = band positions if you already have them. Give [] if not initialized yet.
%  sequence_full  = sequence of RNA.
%  offset         = value that is added to sequence index to achieve 'historical'/favorite numbering. [default: 0]
%  data_type      = cell of tags of modification reactions in each trace, e.g.,
%               {'SHAPE','SHAPE','ddTTP'}
% first_RT_nucleotide = integer that gives nucleotide immediately 5' of primer binding site.
%                         [default: end of sequence]
% structure       = structure in dot/bracket notation [give as '' if unknown] [default: '']
% font_size       = font size of y-axis. Give 6 for printout, supply
%                         10 when you assign bands on the screen.  [optional, default 10]
% JUST_PLOT       = if 1, do not engage in interactive annotation; if 2, disasble window close button. [optional, default 0]
%
% Output:
% xsel      = positions of bands across all lanes.
% seqpos    = sequence numbers that go with each xsel
% area_pred = matrix of zeros and ones that mark band locations for entire assignable sequence.
%
% (C) R. Das, 2013
%

if nargin == 0; help(mfilename); return; end;
ver_flag = check_graphic_interface();

% initialize outputs.
if ~exist('xsel', 'var') || isempty(xsel);  xsel = []; end
seqpos = [];

if length(xsel) > 1 && (xsel(2) > xsel(1));
    fprintf('WARNING! WARNING! WARNING!\n');
    fprintf('You are using the old style of marking, where the bands are marked from 3'' to 5'' stop positions.\n');
    fprintf('This script is switching the order!\n');
    fprintf('Outputted seqpos will go from small to large values.\n');
    fprintf('Outputted xsel will go from large to small values\n');
    xsel = reverse_sort(xsel);
end;

if ~exist('offset', 'var') || isempty(offset);  offset = 0; end;
if ~exist('data_types', 'var') || isempty(data_types); data_types = []; end;
if ~exist('structure', 'var') || isempty(structure);  structure = ''; end;
if ~exist('font_size', 'var') || isempty(font_size); font_size = 10; end;

%
if exist('first_RT_nucleotide','var') || isempty(first_RT_nucleotide);
    sequence = sequence_full(1 : (first_RT_nucleotide-offset));
    if length(structure) > length(sequence); structure = structure(1: length(sequence)); end;
else
    sequence = sequence_full;
end;

% capitalize the sequence...
sequence = upper(sequence);


% fill out area_pred, based on data_types.
numlanes = size(d_align,2);
area_pred = generate_area_pred(sequence, structure, offset, data_types, numlanes);

if ~exist('JUST_PLOT', 'var') || isempty(JUST_PLOT); JUST_PLOT = 0; end

ylim = get(gca,'ylim');
ymin = ylim(1);
ymax = ylim(2);

xlim = get(gca,'xlim');
if (xlim(1) ~= 0.5 || xlim(2) ~= numlanes+0.5);
    ymin = 1;
    ymax = size(d_align,1);
    colormap(1 - gray(100));
end;
axis([ 0.5 numlanes+0.5 ymin ymax]);

scale_factor = 40/ mean(mean(abs(d_align)));
image(d_align * scale_factor);

set(gcf, 'WindowButtonMotionFcn', @mouseMove);

contrast_factor = 100.0/size(colormap,1);
numlanes = size(d_align,2);

stop_sel = 0;
xsel = reverse_sort(xsel);

annotation_handles = {};

update_plot = 1;
update_ylim = 1;
update_contrast = 1;

set(gcf, 'PaperPositionMode', 'auto', 'color', 'white');
if ver_flag; set(gcf, 'pointer', 'fullcross'); end;
%figure_full_screen();
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 0.3 0.9]);
if JUST_PLOT == 2; set(gcf, 'closerequestfcn', ''); end;

while ~stop_sel;
    
    % try to do 'lazy update'
    if update_plot; annotation_handles = make_plot(d_align, xsel, sequence, offset, area_pred, annotation_handles, font_size); end;
    if update_contrast; colormap(1 - gray(round(100/contrast_factor))); end;
    if (update_ylim || update_plot); do_ylim_update(annotation_handles, ymin, ymax); end;
    
    if JUST_PLOT == 1; break; end;
    
    update_plot = 0;
    update_contrast = 0;
    update_ylim = 0;
    
    title(['{\fontsize{14}{\bf{Keys: }}}',...
        '{\fontsize{14}{\color{blue}\bf{left click}}}: Select',...
        '; {\fontsize{14}{\color{red}\bf{middle click, e}}}: Erase Closest',...
        '; {\fontsize{14}{\color[rgb]{0,0.8,0}\bf{up/down, i/k}}}: Move View', ...
        '; {\fontsize{14}{\color{orange}\bf{-/+, j/l}}}: Zoom In/Out', ...
        '; {\fontsize{14}{\color{gray}\bf{1/2}}}: Contrast Dark/Light', char(10), ...
        '{\fontsize{14}{\color{magenta}\bf{q}}}: Save and Quit', ...
        '; {\fontsize{14}{\color{cyan}\bf{r}}}: Reset', ...
        '; {\fontsize{14}{\color[rgb]{0,0.5,0.5}\bf{x}}}: Autoassign (Old)', ...
        '; {\fontsize{14}{\color[rgb]{0.5,0.5,0}\bf{y}}}: Autoassign (New)', char(10), ...
        '{\fontsize{14}{\color[rgb]{0.5,0,0.5}\bf{# of Annotations}}: \bf{', num2str(length(xsel)), ' / ', num2str(length(sequence)), '}}']);

    if ver_flag;
        % This used to be the code, instead of the ginput(1) block above.
        keydown = waitforbuttonpress;
        button = get(gcf, 'SelectionType');
        mouse_pos = get(gca, 'CurrentPoint');
        xselpick = mouse_pos(1,2);
    else
        [~, xselpick, button ]  = ginput(1);
        keydown = (button > 3 );
        if ( button == 1 ) 
           button = 'normal';
        elseif ( button == 2 | button == 3 ) 
           button = 'extend';
        end;
    end;
    
    
    if (~keydown); % mousebutton pressed!
        switch(button);
            case 'normal' % left-click
                xsel = [xsel xselpick];
                xsel = reverse_sort(xsel);
                update_plot = 1;
            case 'extend' % middle click, or shift-click on Mac
                xsel = remove_pick(xsel, xselpick);
                update_plot = 1;
        end;
    else
        
        if ver_flag;
            % This used to be the code to get they keypress -- now, based on
            % button.
            keychar = get(gcf, 'CurrentCharacter');
        else
            keychar = char( button );
        end;

        if double(keychar) == 28; keychar = '+'; end; % left arrow
        if double(keychar) == 29; keychar = '-'; end; % right arrow
        if double(keychar) == 30; keychar = 'i'; end; % up arrow
        if double(keychar) == 31; keychar = 'k'; end; % down arrow
        
        %if double(keychar) == 33; keychar = 'w'; end; % page up
        %f double(keychar) == 34; keychar = 's'; end; % page down
        
        switch keychar;
%             case {'p', 'P'};
%                 [name, path] = uiputfile('xsel.mat', 'Write xsel to matlab file');
%                 if(name);
%                     save(strcat(path, name), 'xsel');
%                 end;
%             case {'o','O'};
%                 [name, path] = uigetfile('*.mat', 'Read xsel from matlab file');
%                 if(name)
%                     a = load(strcat(path, name));
%                     xsel = [xsel a.xsel];
%                     xsel = reverse_sort(xsel);
%                 end;
%                 update_plot = 1;
            case {'q','Q'};
                stop_sel = 1;
                if(~isempty(xsel) && (length(xsel) < length(sequence) || length(xsel) > length(sequence)));
                    %if isstruct(USE_GUI)
                    %btn = questdlg(sprintf('Warning: You have assigned only %d band position(s), which is less/more than the total number of bands (%d). How do you want to proceed?', length(xsel), length(sequence)), ...
                    %	       'Warning!', 'Return to image', 'Force to go','Force to go');
                    %if(~strcmp(btn, 'Force to go'))
                    %  stop_sel = 0;
                    %end
                    %else
                    %fprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n')
                    %fprintf('Warning: You have assigned only %d band position(s), which is less than the total number of bands (%d).\n', length(xsel), length(sequence));
                    %fprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n')
                    %end
                end;
                if JUST_PLOT == 2; set(gcf, 'closerequestfcn', 'closereq'); end;
            case {'e','E'}; % doesn't work -- use mouse middle-click to erase.
                xsel = remove_pick(xsel, xselpick);
                update_plot = 1;
            case {'j','J','+','='};
                %xselpick =  (ymax + ymin)/2;
                current_relative_pos =  (xselpick - ymin)/(ymax-ymin);
                yscale = (ymax - ymin) * 0.75;
                ymin = xselpick - yscale * (current_relative_pos);
                ymax = xselpick + yscale * (1- current_relative_pos);
                update_plot = 1; % for text labels
                update_ylim = 1;
            case {'l','L','-','_'};
                %xselpick =  (ymax + ymin)/2;
                current_relative_pos =  (xselpick - ymin)/(ymax-ymin);
                yscale = (ymax - ymin) / 0.75;
                ymin = xselpick - yscale * (current_relative_pos);
                ymax = xselpick + yscale * (1- current_relative_pos);
                update_plot = 1; % for text labels
                update_ylim = 1;
            case {'i','I'};
                yscale = (ymax - ymin);
                ymin = ymin - yscale * 0.05;
                ymax = ymax - yscale * 0.05;
                update_plot = 1; % for text labels
                update_ylim = 1;
            case {'k','K'};
                yscale = (ymax - ymin);
                ymin = ymin + yscale * 0.05;
                ymax = ymax + yscale * 0.05;
                update_plot = 1; % for text labels
                update_ylim = 1;
            case {'w','W'};
                yscale = (ymax - ymin);
                ymin = ymin - yscale;
                ymax = ymax - yscale;
                update_plot = 1; % for text labels
                update_ylim = 1;
            case {'s','S'};
                yscale = (ymax - ymin);
                ymin = ymin + yscale;
                ymax = ymax + yscale;
                update_plot = 1; % for text labels
                update_ylim = 1;
            case {'b', 'B'};
                yscale = (ymax - ymin);
                ymin = ymin + yscale;
                ymax = ymax + yscale;
                update_plot = 1; % for text labels
                update_ylim = 1;
            case {'t', 'T'};
                yscale = (ymax - ymin);
                ymin = ymin - yscale;
                ymax = ymax - yscale;
                update_plot = 1; % for text labels
                update_ylim = 1;
            case {'1'};
                contrast_factor = contrast_factor * sqrt(2);
                update_contrast = 1;
            case {'2'};
                contrast_factor = contrast_factor / sqrt(2);
                update_contrast = 1;
            case {'r', 'R'};
                xsel = []; % reset
                update_plot = 1;
            case {'x','X'};
                peak_spacing = 0; %size(d_align, 1)/ length( sequence);
                
                if ~exist('area_pred','var') || isempty(area_pred);
                    fprintf('You need to input data_types or area_pred if you want to use auto-assign!\n')
                else
                    input_bounds = [358 5755];
                    if length(xsel) == 2;
                        input_bounds = sort(xsel);
                        peak_spacing = (input_bounds(2) - input_bounds(1)) / length(sequence);
                    end;
                    area_pred_reverse = area_pred(end:-1:1,:); % should fix auto_assign to reverse.
                    fprintf('Running auto-assign. This might take a minute.\n');
                    
                    xsel = auto_assign_sequence(d_align, sequence, area_pred_reverse, peak_spacing, input_bounds, data_types);
                    %                     size(area_pred_reverse)
                    %                     size(d_align)
                    %                     size(sequence)
                    
                    xsel = reverse_sort(xsel); % should fix auto_assign to reverse.
                end;
                update_plot = 1;
            case {'y','Y'};
                param.WINDOW_SIZE = 'auto';
                param.window_jump = 'auto';
                param.begin = 293; %'auto';
                param.end = 7800; %'auto';
                param.peak_bonus_weight = 1; %1.5
                param.peak_range = 'auto';
                param.min_gamma = 'auto';
                param.gamma_add = 0.5;
                param.primary_weight = 'auto';
                param.exist_tailpeak = 1;
                
                if ~exist('area_pred','var') || isempty(area_pred);
                    fprintf('You need to input data_types or area_pred if you want to use auto-assign!\n');
                else
                    area_pred_reverse = area_pred(end:-1:1,:); % should fix auto_assign to reverse.
                    fprintf('Running auto-assign. This might take a minute.\n');
                    
                    [xsel, escore] = auto_assign_cubic (d_align, area_pred_reverse, reverse(sequence), param, data_types);
                    
                    %                     size(area_pred_reverse)
                    %                     size(d_align)
                    %                     size(sequence)
                    %                     escore
                    
                    xsel = reverse_sort(xsel); % should fix auto_assign to reverse.
                end;
                update_plot = 1;
        end;
    end;
    
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if length(xsel) > length(sequence) + 1;
    fprintf('WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!\n');
    fprintf(' Number of selected positions %d exceeds length of sequence %d!\n', length(xsel),length(sequence));
    fprintf('WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!\n');
end;
if length(xsel) < length(sequence);
    fprintf('WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!\n');
    fprintf(' Number of selected positions %d is less than length of sequence %d!\n', length(xsel),length(sequence));
    fprintf('WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!\n');
end;

title '';

seqpos = get_seqpos(sequence, offset, xsel);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xsel = remove_pick(xsel, xselpick)

if ~isempty(xsel);
    [~, closestpick] = min(abs(xsel - xselpick));
    
    xsel = xsel([1:(closestpick-1) ...
        (closestpick+1):length(xsel)] ...
        );
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = reverse_sort(x)
x = sort(x);
x = x(end:-1:1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function do_ylim_update(annotation_handles, ymin, ymax)

for i = 1:length(annotation_handles);
    handles = annotation_handles{i};
    
    m = 3; % text in handle
    if length(handles) < m; continue; end;
    h = handles{m};
    
    pos = get(h, 'Position');
    y = pos(2);
    if (y >= ymin && y <= ymax);
        set(h, 'visible', 'on');
    else
        set(h, 'visible', 'off');
    end;
end;

ylim([ymin ymax]);


function mouseMove(object, eventdata)
C = get(gca, 'CurrentPoint');
%title(gca, ['(X,Y) = (', num2str(C(1,1)), ', ',num2str(C(1,2)), ')']);
