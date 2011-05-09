function varargout = HiTRACE_figure(varargin)
% HITRACE_FIGURE MATLAB code for HiTRACE_figure.fig
%      HITRACE_FIGURE, by itself, creates a new HITRACE_FIGURE or raises the existing
%      singleton*.
%
%      H = HITRACE_FIGURE returns the handle to a new HITRACE_FIGURE or the handle to
%      the existing singleton*.
%
%      HITRACE_FIGURE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in HITRACE_FIGURE.M with the given input arguments.
%
%      HITRACE_FIGURE('Property','Value',...) creates a new HITRACE_FIGURE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before HiTRACE_figure_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to HiTRACE_figure_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help HiTRACE_figure

% Last Modified by GUIDE v2.5 09-May-2011 15:20:48

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @HiTRACE_figure_OpeningFcn, ...
                   'gui_OutputFcn',  @HiTRACE_figure_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before HiTRACE_figure is made visible.
function HiTRACE_figure_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to HiTRACE_figure (see VARARGIN)

% Choose default command line output for HiTRACE_figure
handles.output = hObject;

if(~isdeployed)
    addpath(genpath('../Scripts'));
    addpath(genpath('etc/'));
end

% initialization for settings
settings.refcol = 4;
settings.offset = -10;
settings.ymin = 2200;
settings.ymax = 4200;
settings.dist = 20;
settings.refine = 1;
settings.baseline = 1;

fineTune.slack = 10;
fineTune.shift = 100;
fineTune.windowsize = 100;
fineTune.dpAlign = 1;
fineTune.targetblock = '';

handles.settings = settings;
handles.fineTune = fineTune;
handles.displayComponents = [];

handles.d = [];
handles.da = [];
handles.d_bsub = [];
handles.d_align = [];


handles.filenames = {};
handles.sequence = [];
handles.mutposFile = [];
handles.xsel = [];

handles.step = 0;
handles.stages = {'initial', 'profile', 'baseline', 'dpalign', 'annotation', 'finalresult'};

handles.displayComponents = displaySetting('initial', handles);
menuSetting('initial', handles);
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes HiTRACE_figure wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = HiTRACE_figure_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --------------------------------------------------------------------
function settingMenu_Callback(hObject, eventdata, handles)
% hObject    handle to settingMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function dataMenu_Callback(hObject, eventdata, handles)
% hObject    handle to dataMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename pathname]= uigetfile('*.txt','Pick a list txt file');

if(filename)
    fid = fopen(strcat(pathname,filename), 'r');
    if(fid == -1)
        errordlg('File not found');
    else
        filenames = textscan(fid,'%s');
        handles.filenames = filenames{1};
        
        h = handles.displayComponents(1);
        set(h, 'String', filenames{1},'Max', length(filenames{1}), 'Value', 1);
        
        guidata(hObject, handles);
    end
end



% --------------------------------------------------------------------
function seqMenu_Callback(hObject, eventdata, handles)
% hObject    handle to seqMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
setStatus('Open a sequence',handles);
[filename pathname]= uigetfile('*.txt','Pick a sequence file');
if(filename)
    fid = fopen(strcat(pathname,filename));
    str = textscan(fid, '%s');
    handles.sequence = str{1}{1};
    fclose(fid);
    
    guidata(hObject, handles);
    
    setStatus('Sequence loaded',handles);
    
    for i = handles.displayComponents;
        delete(i);
    end
    
    handles.displayComponents = displaySetting('initial',handles);
    guidata(hObject, handles);
end
 
% --------------------------------------------------------------------
function markMenu_Callback(hObject, eventdata, handles)
% hObject    handle to markMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
setStatus('Open a sequence',handles);
[filename pathname]= uigetfile('*.txt','Pick a sequence file');
if(filename)
    handles.mutposFile = strcat(pathname, filename);
    
    guidata(hObject, handles);
    
    setStatus('Marks guide is set',handles);
    
    for i = handles.displayComponents;
        delete(i);
    end
    
    handles.displayComponents = displaySetting('initial',handles);
    guidata(hObject, handles);
end

% --------------------------------------------------------------------
function optionMenu_Callback(hObject, eventdata, handles)
% hObject    handle to optionMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.settings = GeneralOption(handles.settings);
guidata(hObject, handles);
for i = handles.displayComponents;
    delete(i);
end

handles.displayComponents = displaySetting('initial',handles);
guidata(hObject, handles);

% --------------------------------------------------------------------
function fineOptionMenu_Callback(hObject, eventdata, handles)
% hObject    handle to fineOptionMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.fineTune = FineTunning(handles.fineTune);
guidata(hObject, handles);
for i = handles.displayComponents;
    delete(i);
end

handles.displayComponents = displaySetting('initial',handles);
guidata(hObject, handles);

% --------------------------------------------------------------------
function exitMenu_Callback(hObject, eventdata, handles)
% hObject    handle to exitMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(handles.figure1);

function setStatus(string, handles)
str = sprintf('Status: %s', string);
set(handles.statusLabel,'String',str);


% --------------------------------------------------------------------
function testToolbar_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to testToolbar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp(handles);


function component = displaySetting(mode, handles)
figure(handles.figure1);
switch(mode)
    case 'initial'
        panel = uipanel('Parent', handles.figure1,'Units', 'normalized', 'Position', [0.03, 0.65, 0.94, 0.30], 'Title', 'Filenames');
        hList = uicontrol('Parent', panel, 'Style', 'listbox', 'Units', 'normalized', 'Position', [0.05, 0.05, 0.90, 0.90], 'String', handles.filenames);
        
        
        settings = handles.settings;
        fineTune = handles.fineTune;
        mutposFile = handles.mutposFile;
        sequence = handles.sequence;
        
        if(settings.refine)
            refine='true';
        else
            refine='false';
        end
        
        if(settings.baseline)
            baseline='true';
        else
            baseline='false';
        end
        
        if(fineTune.dpAlign)
            dpAlign='true';
        else
            dpAlign='false';
        end
        
        str = {sprintf('Reference column : %d', settings.refcol);
            sprintf('Offset : %d', settings.offset);
            sprintf('Time range : %d ~ %d', settings.ymin, settings.ymax);
            sprintf('Primer distance from end : %d', settings.dist);
            sprintf('Fitting without refinements : %s', refine);
            sprintf('Enabling baseline subtraction : %s', baseline);
            sprintf('Marks guide: %s', mutposFile);
            sprintf('Sequence: %s',sequence);
            };
            
        
        hOptions = uicontrol('Parent', handles.figure1, 'Style', 'text','HorizontalAlignment','left','Units', 'normalized', 'Position', [0.03, 0.08, 0.56, 0.4], 'String', str);
        
        str = {sprintf('DP alignment : %s', dpAlign);
            sprintf('Slack : %d', fineTune.slack);
            sprintf('Max shift: %d', fineTune.shift);
            sprintf('Window size : %d', fineTune.windowsize);
            sprintf('Target block : %s', fineTune.targetblock);
            };
        
        hDPOptions = uicontrol('Parent', handles.figure1, 'Style', 'text','HorizontalAlignment','left','Units', 'normalized', 'Position', [0.59, 0.08, 0.38, 0.4], 'String', str);
        
        % hAxes = axes('Parent', handles.figure1,'Units', 'normalized', 'Position', [0.35, 0.1, 0.6, 0.85]);
        component = [hList hOptions hDPOptions panel];
    case 'profile'
        h1 = axes('position', [0.1, 0.2, 0.35, 0.7]);
        image( 50 * handles.d);
        axis( [ 0.5 size(handles.da,2)+0.5 handles.settings.ymin handles.settings.ymax] );
        title( 'Signal')
        
        h2 = axes('position', [0.55, 0.2, 0.35, 0.7]);
        image( 50 * handles.da);
        axis( [ 0.5 size(handles.da,2)+0.5 handles.settings.ymin handles.settings.ymax] );
        title( 'Reference ladder')
        
        colormap( 1- gray(100));
        
        component = [h1, h2];
    case 'baseline'
        if(handles.settings.baseline)
            str = 'Signal (baseline subtraction enabled)';
        else
            str = 'Signal (baseline subtraction skipped)';
        end
        
        h1 = axes('Position', [0.1, 0.2, 0.8, 0.7]);
        
        image( 50 * handles.d_bsub);
        axis( [ 0.5 size(handles.d_bsub,2)+0.5 handles.settings.ymin handles.settings.ymax] );        
        title( str );
        
        component = h1;
    case 'dpalign'
        if(handles.fineTune.dpAlign)
            str = 'Signal (nonlinear alignment by DP)';
        else
            str = 'Signal (nonlinear alignment by DP skipped)';
        end
        
        h1 = axes('Position', [0.1, 0.2, 0.8, 0.7]);
        
        image( 50 * handles.d_align);
        axis( [ 0.5 size(handles.d_bsub,2)+0.5 1 size( handles.d_align,1)] )
        title( str );
        
        component = h1;
    case 'annotation'
        h1 = axes('Position', [0.1, 0.2, 0.8, 0.7]);
        
        d_align = handles.d_align / mean( mean( abs( handles.d_align ) ) );
        image( d_align * 20 );
        colormap( 1 - gray(100 ) );

        hold on
        for j = 1:handles.numpeaks
          plot( [0.5  size(handles.d_align,2)+0.5], [handles.xsel(j) handles.xsel(j)],'r');
        end
        if(~isempty(handles.marks))
            for m = 1:size(handles.d_align,2);
	      if ( m <= handles.mutpos )
		goodpos = find( handles.marks(:,1) ==  handles.mutpos(m)  );
		for k = goodpos'
		  which_xsel = handles.numpeaks - handles.marks( k, 2 ) + 1 + handles.settings.offset;
		  if ( which_xsel >= 1 & which_xsel <= length(handles.xsel) )
		    plot( m, handles.xsel(which_xsel),'ro' );
		  end
		end
	      end
            end
        end
        hold off
        
        component = h1;
    case 'finalresult'
        h1 = axes('Position', [0.08, 0.2, 0.27, 0.7]);
        scalefactor = 40/mean(mean(handles.d_align));
        image( scalefactor * handles.d_align );
        title( 'data' );

        h2 = axes('Position', [0.38, 0.2, 0.27, 0.7]);
        image( scalefactor * handles.prof_fit );
        title( 'fit' );

        h3 = axes('Position', [0.68, 0.2, 0.27, 0.7]);
        image( scalefactor * ( handles.d_align - handles.prof_fit) );
        title( 'residuals' );
        
        colormap( 1 - gray(100) );
        component = [h1 h2 h3];
end

function menuSetting(mode, handles)
v = get(handles.uitoolbar1, 'Children');
for i = v
    set(i,'Enable', 'on');
end

v = get(handles.settingMenu, 'Children');
for i = v
    set(i,'Enable', 'on');
end

v = get(handles.fileMenu, 'Children');
for i = v
    set(i,'Enable', 'on');
end

switch(mode)
    case 'initial'
        set(handles.savedataMenu, 'Enable', 'off');
        set(handles.saveresultMenu, 'Enable', 'off');
        set(handles.dataToolbar, 'Enable', 'off');
        set(handles.rdatToolbar, 'Enable', 'off');
        set(handles.nextToolbar, 'Enable', 'off');
        set(handles.prevToolbar, 'Enable', 'off');
    case 'running'
        v = get(handles.uitoolbar1, 'Children');
        for i = v
            set(i,'Enable', 'off');
        end

        v = get(handles.settingMenu, 'Children');
        for i = v
            set(i,'Enable', 'off');
        end
        
        v = get(handles.fileMenu, 'Children');
        for i = v
            set(i,'Enable', 'off');
        end
        
    case 'lastpage'
        set(handles.nextToolbar, 'Enable', 'off');
        set(handles.addlistToolbar, 'Enable', 'off');
        set(handles.removeToolbar, 'Enable', 'off');
    case 'midpage'
        set(handles.addlistToolbar, 'Enable', 'off');
        set(handles.removeToolbar, 'Enable', 'off');
    case 'firstpage'
        set(handles.prevToolbar, 'Enable', 'off');
end

% --------------------------------------------------------------------
function optionToolbar_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to optionToolbar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.settings = GeneralOption(handles.settings);
guidata(hObject, handles);
for i = handles.displayComponents
    delete(i);
end
handles.displayComponents = displaySetting('initial', handles);  
guidata(hObject, handles);

% --------------------------------------------------------------------
function finetuneToolbar_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to finetuneToolbar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% handles    structure with handles and user data (see GUIDATA)
handles.fineTune = FineTunning(handles.fineTune);
guidata(hObject, handles);

for i = handles.displayComponents
    delete(i);
end
handles.displayComponents = displaySetting('initial', handles);  
guidata(hObject, handles);

% --------------------------------------------------------------------
function opendataToolbar_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to opendataToolbar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[name path]= uigetfile('*.mat','Pick a setting file');

if(name)
    fullname = strcat(path, name);
    structure = load(fullname);
    if(isfield(structure, 'settings'))
        handles.settings = structure.settings;
    end
    
    if(isfield(structure, 'fineTune'))
        handles.fineTune = structure.fineTune;
    end
    
    if(isfield(structure, 'sequence'))
        handles.sequence = structure.sequence;
    end
    
    if(isfield(structure, 'mutpos'))
        handles.mutposFile = structure.mutposFile;
    end
    
    if(isfield(structure, 'filenames'))
        handles.filenames = structure.filenames;
    end
        
    guidata(hObject, handles);
    setStatus('Settings loaded', handles);
end

for i = handles.displayComponents
    delete(i);
end
handles.displayComponents = displaySetting('initial', handles);        
guidata(hObject, handles);



% --------------------------------------------------------------------
function saveToolbar_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to saveToolbar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
settings = handles.settings;
fineTune = handles.fineTune;
sequence = handles.sequence;
mutposFile =  handles.mutposFile;
filenames = handles.filenames;

uisave({'settings', 'fineTune', 'sequence', 'mutposFile', 'filenames'}, 'env.mat');

% --------------------------------------------------------------------
function addlistToolbar_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to addlistToolbar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 dir = uigetdir;
 if(dir)
     str = handles.filenames;
     if(isempty(str))
         str = {dir}; 
     else
         str = cat(1, str, {dir});
     end
     
     handles.filenames = str;
     set(handles.displayComponents(1),'String', str, 'Max', length(str), 'Value', 1);
     
     
     guidata(hObject, handles);
 end

% --------------------------------------------------------------------
function removeToolbar_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to removeToolbar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h = handles.displayComponents(1);
names = get(h,'String');
remove = get(h,'Value');

ind = zeros(length(names),1);
if remove > 0 & remove <= length(ind)
    ind(remove) = 1;
    ind = logical(ind);
    names = names(~ind);

    handles.filenames = names;

    set(h, 'String', names, 'Value', 1);

    guidata(hObject, handles);
end


% --------------------------------------------------------------------
function runToolbar_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to runToolbar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Run script
% initialization stage
menuSetting('running', handles);

timeline = tic;

if(isempty(handles.filenames))
    errordlg('Please add one or more data folders. You can select the location of a folder by clicking the ''Add'' button. Alternatively, you can load a text file that contains a list of folders (one path per line).','Error!');
    return;
end

if(isempty(handles.sequence))
    errordlg('Please specify your input sequence. You can load one from a plain text file.','Error!');
    return;
end

setStatus('Initailizing...', handles);


if matlabpool( 'size' ) == 0;
    res = findResource;
    matlabpool( res.ClusterSize );
end

step = handles.step;

filenames = handles.filenames;
ymin = handles.settings.ymin;
ymax = handles.settings.ymax;
dist = handles.settings.dist;
offset = handles.settings.offset;
refine = handles.settings.refine;
refcol = handles.settings.refcol;
baseline = handles.settings.baseline;

slack = handles.fineTune.slack;
shift = handles.fineTune.shift;
windowsize = handles.fineTune.windowsize;
dpAlign = handles.fineTune.dpAlign;

mutposFile = handles.mutposFile;
sequence = handles.sequence;

% target block
s = handles.fineTune.targetblock;
s_string = 'targetblock = {';
r = s;
while( ~isempty(r) )
  [t,r] = strtok( r, ';' );
  s_string = [ s_string,  '[',strrep( t, '-', ':' ),']' ];
  if ~isempty(r); s_string = [ s_string, ',' ]; end;
end
s_string = [s_string, '};' ];
eval( s_string );


for i = step:5
    switch(i)
        case 1
            setStatus('Loading and alignment...', handles);

            [ handles.d, handles.da] = quick_look( filenames, ymin, ymax, [], [], refcol, 0);
            
            for j = handles.displayComponents
                delete(j);
            end
            handles.displayComponents = displaySetting('profile', handles);
            
            guidata(hObject, handles);


        case 2
            if(baseline)
                str = 'Signal (baseline subtraction enabled)';
            else
                str = 'Signal (baseline subtraction skipped)';
            end
            
            setStatus(str, handles);
            
            if(baseline)
                handles.d_bsub = baseline_subtract_v2( handles.d, ymin, ymax, 2e6, 2e4, 0 );
            else
                handles.d_bsub = handles.d;
            end
           
            for j = handles.displayComponents
                delete(j);
            end
            handles.displayComponents = displaySetting('baseline', handles);

            guidata(hObject, handles);
            
        case 3

            setStatus('Nonlinear alignment by DP...', handles);

            if(~exist('targetblock')) targetblock = []; end;
            if( length( targetblock ) == 0 ); targetblock = [ 1: size( handles.d_bsub, 2 ) ]; end;
            if ~iscell( targetblock ) 
                targetblock = { targetblock };
            end

            % need to put in a warning here ... if ddATP ladders, etc., are in align_by_DP step they might be messed up.
            if(~isempty(mutposFile)) 
                data_type = textread(mutposFile,'%s'); % this probably should only occur once...
                for n = 1:length( targetblock )
                    for j = 1:length( targetblock{n} )
                        m = targetblock{n}(j);
                        if ( m <= length(data_type) &  ~isempty( strfind( data_type{m}, 'dd' ) ) )
                            errordlg('You might be applying dynamic programming fine alignment to a lane with a sequencing ladder, and could get weird results','OK'); % Sungroh, can we have an option "change align blocks" that returns to the setting window?
                        end
                    end
                end
            end
            
            if dpAlign
                handles.d_align= align_by_DP( handles.d_bsub(ymin:ymax, : ), targetblock, slack, shift, windowsize, 0);
            else
                handles.d_align = handles.d_bsub(ymin:ymax, : );
            end
            
            for j = handles.displayComponents
                delete(j);
            end
            handles.displayComponents = displaySetting('dpalign', handles);
            
            guidata(hObject,handles);
        case 4
            setStatus('Manual annotation...', handles);
            
            if ~exist( 'data_type' ) data_type = {}; end;
            
            
            seqpos = ( (length(sequence)-dist) : -1 :1 ) + offset;
            
            
            if ~exist( 'structure' );  % later fix this!
                structure = sequence;
                for m = 1:length(sequence); structure(m) = '.';end;
            end
            [ marks, area_pred, mutpos] = get_predicted_marks_SHAPE_DMS_CMCT( structure, sequence, offset, seqpos, data_type);	
            marks = [ marks; (1:length(sequence))'+offset, (1:length(sequence))'+offset ];
            
            xsel = handles.xsel;
	
            if length(xsel) == 0; cla; end;

            period = 1;
            
            for j = handles.displayComponents
                delete(j);
            end
            
            figure(handles.figure1);
            handles.displayComponents = axes('Position', [0.1, 0.2, 0.8, 0.7]);
            
            xsel = mark_sequence( handles.d_align, xsel, sequence(1:end-dist), 0, offset, period, marks, mutpos, handles.displayComponents);
            
            numpeaks = length(xsel);

            handles.xsel = xsel;
            handles.numpeaks = numpeaks;
            handles.mutpos = mutpos;
            handles.marks = marks;
            guidata(hObject, handles);
        case 5
            if(refine)
                str = 'Fitting without refinement';
            else
                str = 'Fitting with refinement';
            end

            setStatus(str, handles);
            
            if(refine)
                [area_peak, prof_fit] = do_the_fit_fast( handles.d_align, xsel', 0.0, 0 );            
            else
                [area_peak, prof_fit] = do_the_fit_GUI( handles.d_align, xsel' );
            end
            
            handles.area_peak = area_peak;
            handles.prof_fit = prof_fit;
            
            for j = handles.displayComponents
                delete(j);
            end
            handles.displayComponents = displaySetting('finalresult', handles);
            guidata(hObject, handles);
    end
end


setStatus('Ready...', handles);

fin=toc(timeline);
str = sprintf('Finished (running time = %f seconds). In order to return to the initial setting window, please close the viewer window.', fin);
h = questdlg(str, 'Finish!', 'Done', 'Done');

menuSetting('lastpage',handles);
handles.step = 5;
guidata(hObject, handles);


function saveResult(handles)
[name path] = uiputfile('Output.rdat', 'Save to RDAT file');
if(name)
    fullname = strcat(path, name);
    seqpos = length(handles.sequence)-handles.settings.dist - [0:(size(handles.area_peak,1)-1)] + handles.settings.offset;
    rdat = fill_rdat(name, handles.sequence, handles.settings.offset, seqpos, handles.area_peak);
    output_rdat_to_file(fullname,rdat);
end


% --------------------------------------------------------------------
function savedataMenu_Callback(hObject, eventdata, handles)
% hObject    handle to savedataMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
d = handles.d;
da = handles.da;
d_bsub = handles.d_bsub;
d_align = handles.d_align;

uisave({'d', 'da', 'd_bsub', 'd_align'}, 'data.mat');

% --------------------------------------------------------------------
function saveresultMenu_Callback(hObject, eventdata, handles)
% hObject    handle to saveresultMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
saveResult(handles);


% --------------------------------------------------------------------
function prevToolbar_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to prevToolbar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
for j = handles.displayComponents
    delete(j);
end

if(handles.step == 0)
    warndlg('This is the first step!');
    handles.displayComponents = displaySetting(handles.stages{handles.step + 1},handles);
    menuSetting('firstpage',handles);
else
    handles.step = handles.step -1;
    handles.displayComponents = displaySetting(handles.stages{handles.step + 1},handles);

    if(handles.step == 0)
        menuSetting('firstpage',handles);
    else
        menuSetting('midpage',handles);
    end
end

guidata(hObject,handles);

% --------------------------------------------------------------------
function nextToolbar_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to nextToolbar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
for j = handles.displayComponents
    delete(j);
end

if(handles.step == 5)
    warndlg('This is the last step!');
    handles.displayComponents = displaySetting(handles.stages{handles.step + 1},handles);
    menuSetting('lastpage',handles);
else
    handles.step = handles.step + 1;
    handles.displayComponents = displaySetting(handles.stages{handles.step + 1},handles);

    if(handles.step == 5)
        menuSetting('lastpage',handles);
    else
        menuSetting('midpage',handles);
    end
end

guidata(hObject,handles);


% --------------------------------------------------------------------
function fileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to fileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function listToolbar_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to listToolbar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename pathname]= uigetfile('*.txt','Pick a list txt file');

if(filename)
    fid = fopen(strcat(pathname,filename), 'r');
    if(fid == -1)
        errordlg('File not found');
    else
        filenames = textscan(fid,'%s');
        handles.filenames = filenames{1};
        
        h = handles.displayComponents(1);
        set(h, 'String', filenames{1},'Max', length(filenames{1}), 'Value', 1);
        
        guidata(hObject, handles);
    end
end

% --------------------------------------------------------------------
function seqToolbar_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to seqToolbar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
setStatus('Open a sequence',handles);
[filename pathname]= uigetfile('*.txt','Pick a sequence file');
if(filename)
    fid = fopen(strcat(pathname,filename));
    str = textscan(fid, '%s');
    handles.sequence = str{1}{1};
    fclose(fid);
    
    guidata(hObject, handles);
    
    setStatus('Sequence loaded',handles);
    
    for i = handles.displayComponents;
        delete(i);
    end
    
    handles.displayComponents = displaySetting('initial',handles);
    guidata(hObject, handles);
end

% --------------------------------------------------------------------
function marksToolbar_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to marksToolbar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
setStatus('Open a sequence',handles);
[filename pathname]= uigetfile('*.txt','Pick a sequence file');
if(filename)
    handles.mutposFile = strcat(pathname, filename);
    
    guidata(hObject, handles);
    
    setStatus('Marks guide is set',handles);
    
    for i = handles.displayComponents;
        delete(i);
    end
    
    handles.displayComponents = displaySetting('initial',handles);
    guidata(hObject, handles);
end

% --------------------------------------------------------------------
function rdatToolbar_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to rdatToolbar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
saveResult(handles);

% --------------------------------------------------------------------
function dataToolbar_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to dataToolbar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
d = handles.d;
da = handles.da;
d_bsub = handles.d_bsub;
d_align = handles.d_align;

uisave({'d', 'da', 'd_bsub', 'd_align'}, 'data.mat');
