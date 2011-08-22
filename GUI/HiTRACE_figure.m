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

% Last Modified by GUIDE v2.5 22-Aug-2011 11:10:04

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

axes(handles.logoAxes);
imshow('hitrace_logo.png');

% initialization for settings
settings.refcol = 4;
settings.offset = -10;
settings.ymin = 2200;
settings.ymax = 4200;
settings.dist = 20;
settings.refine = 1;
settings.baseline = 1;
settings.eternaStart = 6;
settings.eternaEnd = 6;
settings.eternaCheck = 0;
settings.peakspacing = 24;

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
handles.sequenceFile = [];
handles.structure = [];
handles.mutposFile = [];
handles.xsel = [];

handles.step = 0;
handles.stages = {'initial', 'profile', 'baseline', 'dpalign', 'annotation', 'finalresult', 'score'};

handles.displayComponents = displaySetting('initial', handles);
menuSetting('initial', handles);
refreshSetting(handles);
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
function seqMenu_Callback(hObject, eventdata, handles)
% hObject    handle to seqMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
setStatusLabel('Open a sequence',handles);
[filename pathname]= uigetfile('*.txt','Pick a sequence file');
if(filename)
    fid = fopen(strcat(pathname,filename));
    str = textscan(fid, '%s');
    handles.sequence = str{1}{1};
    fclose(fid);
    
    guidata(hObject, handles);
    
    setStatusLabel('Sequence loaded',handles);
    
    try
        for i = handles.displayComponents;
            delete(i);
        end
    catch err
        close(handles.figure1);
    end
    
    handles.displayComponents = displaySetting('initial',handles);
    guidata(hObject, handles);
end
 


function component = displaySetting(mode, handles)
figure(handles.figure1);
switch(mode)
    case 'initial'
        component = [];
    case 'profile'
        h1 = axes('position', [0.5 0.2 0.2 0.7]);
        image( 50 * handles.d);
        axis( [ 0.5 size(handles.da,2)+0.5 handles.settings.ymin handles.settings.ymax] );
        title( 'Signal')
        
        h2 = axes('position', [0.75 0.2 0.2 0.7]);
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
        
        h1 = axes('Position', [0.5 0.2 0.45 0.7]);
        
        if(handles.settings.eternaCheck)
            index = get(handles.profileCombo,'Value');
            image( 50 * handles.d_bsub{index});
            axis( [ 0.5 size(handles.d_bsub{index},2)+0.5 1 size( handles.d_bsub{index},1)] );        
            title( str );
        else
            image( 50 * handles.d_bsub);
            axis( [ 0.5 size(handles.d_bsub,2)+0.5 handles.settings.ymin handles.settings.ymax] );        
            title( str );
        end
        
        component = h1;
    case 'dpalign'
        if(handles.fineTune.dpAlign)
            str = 'Signal (nonlinear alignment by DP)';
        else
            str = 'Signal (nonlinear alignment by DP skipped)';
        end
        
        h1 = axes('Position', [0.5 0.2 0.45 0.7]);
        
        image( 50 * handles.d_align);
        axis( [ 0.5 size(handles.d_align,2)+0.5 1 size( handles.d_align,1)] )
        title( str );
        
        component = h1;
    case 'annotation'
        if(handles.settings.eternaCheck)
            h1 = axes('Position', [0.5 0.2 0.45 0.7]);

            index = get(handles.profileCombo,'Value');
            
            d_bsub{index} = handles.d_bsub{index} / mean( mean( abs( handles.d_bsub{index}) ) );
            image( d_bsub{index} * 20 );
            colormap( 1 - gray(100 ) );

            hold on
            for j = 1:handles.numpeaks{index}
              plot( [0.5  size(handles.d_bsub{index},2)+0.5], [handles.xsel{index}(j) handles.xsel{index}(j)],'r');
            end
            if(~isempty(handles.marks{index}))
                for m = 1:size(handles.d_bsub{index},2);
              if ( m <= handles.mutpos{index} )
            goodpos = find( handles.marks{index}(:,1) ==  handles.mutpos{index}(m)  );
            for k = goodpos'
              which_xsel = handles.numpeaks{index} - handles.marks{index}( k, 2 ) + 1 + handles.settings.offset;
              if ( which_xsel >= 1 & which_xsel <= length(handles.xsel{index}) )
                plot( m, handles.xsel{index}(which_xsel),'ro' );
              end
            end
              end
                end
            end
            hold off
        else
            h1 = axes('Position', [0.5 0.2 0.45 0.7]);

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
        end
        
        component = h1;
    case 'finalresult'
        if(handles.settings.eternaCheck)
            index = get(handles.profileCombo,'Value');
            
            h1 = axes('Position', [0.5 0.2 0.15 0.7]);
            scalefactor = 40/mean(mean(handles.d_bsub{index}));
            image( scalefactor * handles.d_bsub{index} );
            title( 'data' );

            h2 = axes('Position', [0.67 0.2 0.15 0.7]);
            image( scalefactor * handles.prof_fit{index} );
            title( 'fit' );

            h3 = axes('Position', [0.84 0.2 0.15 0.7]);
            image( scalefactor * ( handles.d_bsub{index} - handles.prof_fit{index}) );
            title( 'residuals' );
        else
            h1 = axes('Position', [0.5 0.2 0.15 0.7]);
            scalefactor = 40/mean(mean(handles.d_align));
            image( scalefactor * handles.d_align );
            title( 'data' );

            h2 = axes('Position', [0.67 0.2 0.15 0.7]);
            image( scalefactor * handles.prof_fit );
            title( 'fit' );

            h3 = axes('Position', [0.84 0.2 0.15 0.7]);
            image( scalefactor * ( handles.d_align - handles.prof_fit) );
            title( 'residuals' );
        end
        
        colormap( 1 - gray(100) );
        component = [h1 h2 h3];
    case 'score'
        index = get(handles.profileCombo,'Value');
        index2 = get(handles.dataCombo, 'Value');
        area_bsub = handles.area_bsub;
        nres = size( area_bsub{index}, 1);

        goodbins = [(nres-handles.settings.eternaEnd):-1:handles.settings.eternaStart+1]; % have to go backwards. Silly convention switch.

        % average over appropriate columns...
        data = mean(area_bsub{index}( goodbins, index2 ), 2)';

        % normalize.
        [data_norm, ~] = SHAPE_normalize( data );
        pred = handles.all_area_pred{index}( goodbins, 2 ); 

        n = length( data_norm );

        f  = [   0,           0,    ones(1,n),  ones(1,n), zeros( 1, n )           ];

        LB(1) = 0.60 / mean( max(data_norm,0.0) );

        data_sort=sort(data_norm);

        UB(1) = inf;

        LB(2) = -0.1;
        UB(2) =  0.1;

        LB( 2 + [1:n] ) = zeros(1,n);
        UB( 2 + [1:n] ) = inf * ones(1,n);

        LB( 2 + n + [1:n] ) = zeros(1,n);
        UB( 2 + n + [1:n] ) = inf * ones(1,n);

        pred_high = find(  pred );

        LB( 2 + 2*n + pred_high) =  0.0;
        UB( 2 + 2*n + pred_high) =  0.0;

        pred_low  = find( ~pred);
        LB( 2 + 2*n + pred_low) = -0.0;
        UB( 2 + 2*n + pred_low) = +0.0;

        Aeq = [ data_norm', ones(n,1), eye(n,n), -eye(n,n), eye(n,n) ]; 

        beq = pred;

        options = optimset('Display','off');
        params = linprog( f, [], [], Aeq, beq, LB, UB,[],options);


        h2 = axes('position', [0.5 0.6 0.4 0.3]);
        plot( pred, 'k','linewidth',2 );
        hold on;
        plot( params(1)*data_norm + params(2), 'kx-' ); 
        plot( params(1)*data_norm + params(2) + params( 2+2*n+[1:n])', 'rx-' ); 
        hold off;

        h1 = axes('position', [0.5 0.2 0.4 0.3]);
        scale_factor =  params(1);
        baseline     =  params(2);
        max_SHAPE  = (1.0 - baseline)/scale_factor;
        min_SHAPE =  (0.0 - baseline)/scale_factor;
        threshold_SHAPE = ( 0.5 - baseline)/scale_factor;

        bar( data_norm,'facecolor',[0.7 0.7 0.7],'edgecolor',[0 0 0],'barwidth',0.5 );
        hold on
        plot( [1 length(data_norm)], [max_SHAPE max_SHAPE], 'k' );
        plot( [1 length(data_norm)], [min_SHAPE min_SHAPE], 'k' );
        plot( [1 length(data_norm)], [threshold_SHAPE threshold_SHAPE], 'k' );
        plot( (pred - baseline)/scale_factor,'color',[0.5 0 0] );

        ETERNA_score = 0.0;
        badpt = [];
        for k = 1:length( data_norm )
          if ( pred(k) ) % should be unpaired
            if ( data_norm(k) > (0.25*threshold_SHAPE + 0.75*min_SHAPE ) );    
              ETERNA_score  = ETERNA_score + 1;
            else
              badpt = [badpt, k ];
            %  [k,pred(k)]
            end
          else
            if ( data_norm(k) < threshold_SHAPE )
              ETERNA_score  = ETERNA_score + 1;
            else
              badpt = [badpt, k ];
            %  [k,pred(k)]
            end
          end
        end
        plot( badpt, (pred(badpt)-baseline)/scale_factor, 'x','color',[0. 0.5 0],'linewidth',2 );
        ylim([ (-0.2-baseline)/scale_factor (1.5-baseline)/scale_factor ]);
        hold off;

        component = [h1, h2];
end

function menuSetting(mode, handles)
v = get(handles.uitoolbar1, 'Children');
for i = v
    set(i,'Enable', 'on');
end

v = get(handles.listPanel, 'Children');
for i = v
    set(i,'Enable', 'on');
end

v = get(handles.seqPanel, 'Children');
for i = v
    set(i,'Enable', 'on');
end

v = get(handles.optionPanel, 'Children');
for i = v
    set(i,'Enable', 'on');
end

v = get(handles.dpPanel, 'Children');
for i = v
    set(i,'Enable', 'on');
end

v = get(handles.ioPanel, 'Children');
for i = v
    set(i,'Enable', 'on');
end

set(handles.runBtn, 'Enable', 'on');
set(handles.resetBtn, 'Enable', 'on');
set(handles.prevBtn, 'Enable', 'on');
set(handles.nextBtn, 'Enable', 'on');
set(handles.profileCombo, 'Enable', 'on');
set(handles.dataCombo, 'Enable', 'on');


switch(mode)
    case 'initial'
        set(handles.nextToolbar, 'Enable', 'off');
        set(handles.prevToolbar, 'Enable', 'off');
        set(handles.prevBtn, 'Enable', 'off');
        set(handles.nextBtn, 'Enable', 'off');
        set(handles.rdatBtn, 'Enable', 'off');
        set(handles.saveworkspaceBtn, 'Enable', 'off');
        set(handles.profileCombo, 'Enable', 'off');
        set(handles.dataCombo, 'Enable', 'off');
        set(handles.startEdit, 'Enable', 'off');
        set(handles.endEdit, 'Enable', 'off');
	case 'running'
        v = get(handles.uitoolbar1, 'Children');
        for i = v
            set(i,'Enable', 'off');
        end
        
        v = get(handles.listPanel, 'Children');
        for i = v
            set(i,'Enable', 'off');
        end
        
        v = get(handles.seqPanel, 'Children');
        for i = v
            set(i,'Enable', 'off');
        end
        
        v = get(handles.optionPanel, 'Children');
        for i = v
            set(i,'Enable', 'off');
        end
        
        v = get(handles.dpPanel, 'Children');
        for i = v
            set(i,'Enable', 'off');
        end
        
        v = get(handles.ioPanel, 'Children');
        for i = v
            set(i,'Enable', 'off');
        end
        
        set(handles.runBtn, 'Enable', 'off');
        set(handles.prevBtn, 'Enable', 'off');
        set(handles.nextBtn, 'Enable', 'off');
        set(handles.profileCombo, 'Enable', 'off');
        set(handles.dataCombo, 'Enable', 'off');
        set(handles.resetBtn, 'Enable', 'off');
    case 'lastpage'
        set(handles.nextToolbar, 'Enable', 'off');
        set(handles.nextBtn, 'Enable', 'off');
        
        if(handles.settings.eternaCheck)
            set(handles.profileCombo, 'Enable', 'on');
            set(handles.dataCombo, 'Enable', 'on');
        else
            set(handles.profileCombo, 'Enable', 'off');
            set(handles.dataCombo, 'Enable', 'off');
        end
    case 'midpage'
        if(handles.settings.eternaCheck)
            set(handles.profileCombo, 'Enable', 'on');
            set(handles.dataCombo, 'Enable', 'on');
        else
            set(handles.profileCombo, 'Enable', 'off');
            set(handles.dataCombo, 'Enable', 'off');
        end
    case 'firstpage'
        set(handles.prevToolbar, 'Enable', 'off');
        set(handles.prevBtn, 'Enable', 'off');
        if(handles.settings.eternaCheck)
            set(handles.profileCombo, 'Enable', 'on');
            set(handles.dataCombo, 'Enable', 'on');
        else
            set(handles.profileCombo, 'Enable', 'off');
            set(handles.dataCombo, 'Enable', 'off');
        end
end

% --------------------------------------------------------------------
function optionToolbar_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to optionToolbar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.settings = GeneralOption(handles.settings);
guidata(hObject, handles);
try
    for i = handles.displayComponents
        delete(i);
    end
catch err
    close(handles.figure1);
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

try
    for i = handles.displayComponents
        delete(i);
    end
catch err
    close(handles.figure1);
end
handles.displayComponents = displaySetting('initial', handles);  
guidata(hObject, handles);

% --------------------------------------------------------------------
function openenvToolbar_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to openenvToolbar (see GCBO)
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
    setStatusLabel('Settings loaded', handles);
end


try
    for i = handles.displayComponents
        delete(i);
    end
catch err
    close(handles.figure1);
end
guidata(hObject, handles);
refreshSetting(handles);


% --------------------------------------------------------------------
function saveenvToolbar_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to saveenvToolbar (see GCBO)
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


function saveResult(handles)
[name path] = uiputfile('Output.rdat', 'Save to RDAT file');
if(name)
    if(handles.settings.eternaCheck)
        seqpos = [length(handles.sequence{1})-20: -1: 1 ];
        nres = size( handles.area_bsub{1}, 1);
        goodbins = [(nres-handles.settings.eternaEnd):-1:handles.settings.eternaStart+1];
        data_type = {'SHAPE','SHAPE','nomod','ddTTP'}; 
        added_salt = {'chemical:MgCl2:10mM','chemical:NaCl:1M','',''};
        bad_lanes = [];
        comments = {'Chemical mapping data for the EteRNA design project.'};
        eterna_create_rdat_files( strcat(path,name), handles.target_names{1}, handles.structure, handles.sequence, ...
            handles.area_bsub_norm, {get(handles.sequenceEdit, 'String')}, seqpos, goodbins, data_type, added_salt, ...
            bad_lanes, handles.min_SHAPE, handles.max_SHAPE, handles.threshold_SHAPE, handles.ETERNA_score, comments, ...
            handles.darea_bsub_norm );
    else
        fullname = strcat(path, name);
        seqpos = length(handles.sequence)-handles.settings.dist - [0:(size(handles.area_peak,1)-1)] + handles.settings.offset;
        rdat = fill_rdat(name, handles.sequence, handles.settings.offset, seqpos, handles.area_peak);
        output_rdat_to_file(fullname,rdat);
    end
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
try
    for j = handles.displayComponents
        delete(j);
    end
catch err
    close(handles.figure1);
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
try
    for j = handles.displayComponents
        delete(j);
    end
catch err
    close(handles.figure1);
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


% --------------------------------------------------------------------
function marksToolbar_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to marksToolbar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
setStatusLabel('Open a mark guide',handles);
[filename pathname]= uigetfile('*.txt','Pick a mark guide file');
if(filename)
    handles.mutposFile = strcat(pathname, filename);
    
    guidata(hObject, handles);
    
    setStatusLabel('Marks guide is set',handles);
    try
        for i = handles.displayComponents;
            delete(i);
        end
    catch err
        close(handles.figure1);
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


% --- Executes on button press in runBtn.
function runBtn_Callback(hObject, eventdata, handles)
% hObject    handle to runBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = setupSetting(handles);
guidata(hObject, handles);

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

setStatusLabel('Initailizing...', handles);


if exist('mablabpool')
    if matlabpool( 'size' ) == 0;
        res = findResource;
        matlabpool( res.ClusterSize );
    end
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

eternaCheck = handles.settings.eternaCheck;
eternaStart = handles.settings.eternaStart;
eternaEnd = handles.settings.eternaEnd;

peak_spacing = handles.settings.peakspacing;

slack = handles.fineTune.slack;
shift = handles.fineTune.shift;
windowsize = handles.fineTune.windowsize;
dpAlign = handles.fineTune.dpAlign;

mutposFile = handles.mutposFile;
sequence = handles.sequence;
structure = handles.structure;

num_sequences = length( sequence );
which_sets = 1 : num_sequences;

period= 1;

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


for i = step:handles.max
    switch(i)
        case 1
            setStatusLabel('Loading and alignment...', handles);

            [ handles.d, handles.da] = quick_look( filenames, ymin, ymax, [], [], refcol, 0);
            
            for j = handles.displayComponents
                delete(j);
            end
            handles.displayComponents = displaySetting('profile', handles);
            
            guidata(hObject, handles);


        case 2
            if(eternaCheck)
                str = 'DP Alignment using reference column';
                
                if(~exist('targetblock')) targetblock = []; end;
                if( length( targetblock ) == 0 ); targetblock = [ 1: size( handles.d_bsub, 2 ) ]; end;
                if ~iscell( targetblock ) 
                    targetblock = { targetblock };
                end
                
                setStatusLabel(str, handles);
                d_align = align_by_DP_using_ref( handles.d(ymin:ymax,:), handles.da(ymin:ymax,:), [], slack, shift, windowsize, 0);
                
                handles.d_align = d_align;
                
                for j = handles.displayComponents
                    delete(j);
                end
                
                handles.displayComponents = displaySetting('dpalign', handles);
            else
                if(baseline)
                    str = 'Signal (baseline subtraction enabled)';
                else
                    str = 'Signal (baseline subtraction skipped)';
                end

                setStatusLabel(str, handles);

                if(baseline)
                    handles.d_bsub = baseline_subtract_v2( handles.d, ymin, ymax, 2e6, 2e4, 0 );
                else
                    handles.d_bsub = handles.d;
                end

                for j = handles.displayComponents
                    delete(j);
                end
                
                handles.displayComponents = displaySetting('baseline', handles);
            end
            
            guidata(hObject, handles);
            
        case 3
            if(eternaCheck)
                if(baseline)
                    for j = which_sets;
                      str = sprintf('Signal (baseline subtraction enabled) ( %d / %d )', j, which_sets(end));
                      setStatusLabel(str, handles);
                      d_r{j} = handles.d_align( : ,j + ([1:4] - 1)*num_sequences );
                      handles.d_bsub{j} = baseline_subtract_v2( d_r{j} );
                    end
                else
                    which_sets = 1 : num_sequences; 
                    for j = which_sets;
                      str = sprintf('Signal (baseline subtraction skipped) ( %d / %d )', j, which_sets(end));
                      setStatusLabel(str, handles);
                      d_r{j} = handles.d_align( : ,j + ([1:4] - 1)*num_sequences );
                      handles.d_bsub{j} = d_r{j};
                    end
                    
                end

                for j = handles.displayComponents
                    delete(j);
                end
                
                handles.displayComponents = displaySetting('baseline', handles);
            else
                setStatusLabel('Nonlinear alignment by DP...', handles);

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
                                % errordlg('You might be applying dynamic programming fine alignment to a lane with a sequencing ladder, and could get weird results','OK'); % Sungroh, can we have an option "change align blocks" that returns to the setting window?
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
            end
            guidata(hObject,handles);
        case 4
            setStatusLabel('Manual annotation...', handles);
            
            if ~exist( 'data_type' ) data_type = {}; end;
            
            if(~isempty(mutposFile)) 
                data_type = textread(mutposFile,'%s'); % this probably should only occur once...
                for n = 1:length( targetblock )
                    for j = 1:length( targetblock{n} )
                        m = targetblock{n}(j);
                        if ( m <= length(data_type) &  ~isempty( strfind( data_type{m}, 'dd' ) ) )
                            % errordlg('You might be applying dynamic programming fine alignment to a lane with a sequencing ladder, and could get weird results','OK'); % Sungroh, can we have an option "change align blocks" that returns to the setting window?
                        end
                    end
                end
            end
            
            seqpos = ( (length(sequence)-dist) : -1 :1 ) + offset;
            
            if(eternaCheck)
                data_types = {'SHAPE','SHAPE','nomod','ddTTP'};
                for j = which_sets  
                  seqpos = length(sequence{j})-20 - [1:(length(sequence{j})-20)] + 1 + offset;
                  [ marks{j}, all_area_pred{j}, mutpos{j} ] = get_predicted_marks_SHAPE_DMS_CMCT( structure, sequence{j}, offset , seqpos, data_types );
                end

                if(verLessThan('matlab', '7.10.0'))
                    for j = which_sets;
                      setStatusLabel(sprintf('Auto assign annotation... ( %d / %d )',j,which_sets(end)), handles);
                      xsel{j} = auto_assign_sequence( handles.d_bsub{j}, sequence{j}(1:end-dist), seqpos, offset, all_area_pred{j}, peak_spacing, [], 0 );
                    end
                else
                    setStatusLabel('Auto assign annotation parallely...', handles);
                    parfor j = which_sets;
                      xsel{j} = auto_assign_sequence( handles.d_bsub{j}, sequence{j}(1:end-dist), seqpos, offset, all_area_pred{j}, peak_spacing, [], 0 );
                    end
                end
                
                
                for j = handles.displayComponents
                    delete(j);
                end
                
                figure(handles.figure1);
                handles.displayComponents = axes('Position', [0.5 0.2 0.45 0.7]);
                
                for j = which_sets;
                    setStatusLabel(sprintf('Manual annotation... ( %d / %d )', j, which_sets(end)), handles);
                    set(handles.profileCombo, 'Value', j);
                    xsel{j} = mark_sequence( handles.d_bsub{j}, xsel{j}, sequence{j}(1:end-dist), 0, offset, period, marks{j}, mutpos{j}, all_area_pred{j},handles.displayComponents);
                    numpeaks{j} = length(xsel{j});
                end
                
                handles.xsel = xsel;
                handles.numpeaks = numpeaks;
                handles.mutpos = mutpos;
                handles.marks = marks;
                
                handles.all_area_pred = all_area_pred;
                
            else
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
                handles.displayComponents = axes('Position', [0.5 0.2 0.45 0.7]);

                xsel = mark_sequence( handles.d_align, xsel, sequence(1:end-dist), 0, offset, period, marks, mutpos, area_pred, handles.displayComponents);

                numpeaks = length(xsel);

                handles.xsel = xsel;
                handles.numpeaks = numpeaks;
                handles.mutpos = mutpos;
                handles.marks = marks;
            end
            guidata(hObject, handles);
        case 5
            if(refine)
                str = 'Fitting without refinement';
            else
                str = 'Fitting with refinement';
            end

            setStatusLabel(str, handles);
            
            if(eternaCheck)
                for j = which_sets
                    if(refine)
                        [area_peak{j}, prof_fit{j}] = do_the_fit_FAST( handles.d_bsub{j}, handles.xsel{j}', 0.0, 0);
                    else
                        [area_peak{j}, prof_fit{j}] = do_the_fit_GUI( handles.d_bsub{j}, handles.xsel{j}' );
                    end
                end
            else
                if(refine)
                    [area_peak, prof_fit] = do_the_fit_fast( handles.d_align, xsel', 0.0, 0 );            
                else
                    [area_peak, prof_fit] = do_the_fit_GUI( handles.d_align, xsel' );
                end
            end
            
            handles.area_peak = area_peak;
            handles.prof_fit = prof_fit;
            
            for j = handles.displayComponents
                delete(j);
            end
            handles.displayComponents = displaySetting('finalresult', handles);
            guidata(hObject, handles);
        case 6
            backgd_sub_col = 3;
            penalize_negative_weight = 2.0;
            for j = which_sets
              [ area_bsub{j}, darea_bsub{j}] = overmod_and_background_correct_logL( handles.area_peak{j}, backgd_sub_col, [4:size(handles.area_peak{j},1)-4], handles.all_area_pred{j});
             % pause;
            end

            % Average data over Ann/DC replicates, and output to screen in nice format that can be copy/pasted into EteRNA.
            % Also, automated EteRNA scoring.
            START = eternaStart; % where to start data for eterna input
            END   = eternaEnd;   % where to end data for eterna input.
            min_SHAPE = {}; max_SHAPE={};threshold_SHAPE={};ETERNA_score={};
            area_bsub_norm  = area_bsub;
            for j = which_sets
              for a = 1:2  

              nres = size( area_bsub{j}, 1);

              goodbins = [(nres-END):-1:START+1]; % have to go backwards. Silly convention switch.
              data_cols = [a];

              % average over appropriate columns...
              data = mean(area_bsub{j}( goodbins, data_cols ), 2)';

              % normalize.
              [data_norm, scalefactor] = SHAPE_normalize( data );
              area_bsub_norm{j}(:,data_cols)  = area_bsub{j}(:,data_cols)/scalefactor;
              darea_bsub_norm{j}(:,data_cols) = darea_bsub{j}(:,data_cols)/scalefactor;

              % To compute EteRNA score, need predicted paired/unpaired for 'perfect' design.
              % this was computed above to aid in sequence annotation, background subtraction, etc.
              pred = handles.all_area_pred{j}( goodbins, 2 ); 
              [min_SHAPE{j,a}, max_SHAPE{j,a}, threshold_SHAPE{j,a}, ETERNA_score{j,a} ] = determine_thresholds_and_ETERNA_score( data_norm, pred );
              %pause;
              end
            end
            
            handles.area_bsub = area_bsub;
            handles.darea_bsub = darea_bsub;
            
            handles.area_bsub_norm = area_bsub_norm;
            handles.darea_bsub_norm = darea_bsub_norm;
            
            handles.min_SHAPE = min_SHAPE;
            handles.max_SHAPE = max_SHAPE;
            handles.threshold_SHAPE = threshold_SHAPE;
            handles.ETERNA_score = ETERNA_score;
            
            for j = handles.displayComponents
                delete(j);
            end
            handles.displayComponents = displaySetting('score', handles);
            guidata(hObject, handles);
            close(3);
    end
end


setStatusLabel('Ready...', handles);

fin=toc(timeline);
str = sprintf('Finished (running time = %f seconds). In order to return to the initial setting window, please close the viewer window.', fin);
h = questdlg(str, 'Finish!', 'Done', 'Done');

menuSetting('lastpage',handles);
handles.step = handles.max;
guidata(hObject, handles);


% --- Executes on button press in rdatBtn.
function rdatBtn_Callback(hObject, eventdata, handles)
% hObject    handle to rdatBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
saveResult(handles);

% --- Executes on button press in saveworkspaceBtn.
function saveworkspaceBtn_Callback(hObject, eventdata, handles)
% hObject    handle to saveworkspaceBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
settings = handles.settings;
fineTune = handles.fineTune;
sequence = handles.sequence;
sequenceFile = handles.sequenceFile;
mutposFile =  handles.mutposFile;
filenames = handles.filenames;
structure = handles.structure;


d = handles.d;
da = handles.da;
d_bsub = handles.d_bsub;
d_align = handles.d_align;
xsel = handles.xsel;
numpeaks = handles.numpeaks;
marks = handles.marks;        
mutpos = handles.mutpos;
prof_fit = handles.prof_fit;

if(handles.settings.eternaCheck)
    area_bsub = handles.area_bsub;
    darea_bsub = handles.darea_bsub;
    area_bsub_norm = handles.area_bsub_norm;
    darea_bsub_norm = handles.darea_bsub_norm;
    min_SHAPE = handles.min_SHAPE;
    max_SHAPE = handles.max_SHAPE;
    threshold_SHAPE = handles.threshold_SHAPE;
    ETERNA_score = handles.ETERNA_score;

    ids = handles.ids;
    target_names = handles.target_names;
    subrounds = handles.subrounds;
    design_names = handles.design_names;
    
    all_area_pred = handles.all_area_pred;
    uisave({'settings', 'fineTune', 'sequence', 'sequenceFile', 'structure', 'mutposFile', ...
        'filenames', 'd', 'da' , 'd_bsub', 'd_align', 'xsel', 'numpeaks', 'marks', 'mutpos', ...
        'prof_fit','area_bsub', 'darea_bsub', 'area_bsub_norm', 'darea_bsub_norm', 'min_SHAPE', ...
        'max_SHAPE', 'threshold_SHAPE', 'ETERNA_score', 'ids', 'target_names', 'subrounds', ...
        'design_names', 'all_area_pred'}, 'workspace.mat');
else
    uisave({'settings', 'fineTune', 'sequence', 'sequenceFile', 'structure', 'mutposFile', 'filenames', 'd', 'da' , 'd_bsub', 'd_align', 'xsel', 'numpeaks', 'marks', 'mutpos', 'prof_fit'}, 'workspace.mat');    
end

% --- Executes on button press in loadworkspaceBtn.
function loadworkspaceBtn_Callback(hObject, eventdata, handles)
% hObject    handle to loadworkspaceBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[name path]= uigetfile('*.mat','Pick a workspace file');

if(name)
    fullname = strcat(path, name);
    structure = load(fullname);
    if(prod(double(isfield(structure, {'settings', 'fineTune', 'sequence', 'sequenceFile', 'structure', 'mutposFile', 'filenames', 'd', 'da' , 'd_bsub', 'd_align', 'xsel', 'numpeaks', 'marks', 'mutpos', 'prof_fit'}))))
        handles.settings = structure.settings;
        handles.fineTune = structure.fineTune;
        handles.sequence = structure.sequence;
        handles.sequenceFile = structure.sequenceFile;
        handles.mutposFile = structure.mutposFile;
        handles.filenames = structure.filenames;
        handles.structure = structure.structure;
        
        handles.d = structure.d;
        handles.da = structure.da;
        handles.d_bsub = structure.d_bsub;
        handles.d_align = structure.d_align;
        
        handles.xsel = structure.xsel;
        handles.numpeaks = structure.numpeaks;
        handles.marks = structure.marks;        
        handles.mutpos = structure.mutpos;
        handles.prof_fit = structure.prof_fit;
        
        if(handles.settings.eternaCheck)
            handles.area_bsub = structure.area_bsub;
            handles.darea_bsub = structure.darea_bsub;
            handles.area_bsub_norm = structure.area_bsub_norm;
            handles.darea_bsub_norm = structure.darea_bsub_norm;
            handles.min_SHAPE = structure.min_SHAPE;
            handles.max_SHAPE = structure.max_SHAPE;
            handles.threshold_SHAPE = structure.threshold_SHAPE;
            handles.ETERNA_score = structure.ETERNA_score;

            handles.ids = structure.ids;
            handles.target_names = structure.target_names;
            handles.subrounds = structure.subrounds;
            handles.design_names = structure.design_names;

            handles.all_area_pred = structure.all_area_pred;
            
            set(handles.profileCombo,'String', handles.design_names);
            set(handles.dataCombo,'String', {'1', '2'});
        end
        
        handles.step = 0;
        
        for i = handles.displayComponents
            delete(i);
        end
        handles.displayComponents = displaySetting('initial',handles);
        menuSetting('firstpage',handles);
        refreshSetting(handles);
        
        if(handles.settings.eternaCheck)
            handles.max = 6;
        else
            handles.max = 5;
        end
        
        guidata(hObject, handles);
        setStatusLabel('Workspace loaded', handles);
    else
        errordlg('Invalid workspace file.', 'Error!');
    end
end

function guideEdit_Callback(hObject, eventdata, handles)
% hObject    handle to guideEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of guideEdit as text
%        str2double(get(hObject,'String')) returns contents of guideEdit as a double


% --- Executes during object creation, after setting all properties.
function guideEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to guideEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in guideBtn.
function guideBtn_Callback(hObject, eventdata, handles)
% hObject    handle to guideBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
setStatusLabel('Open a mark guides file',handles);
[filename pathname]= uigetfile('*.txt','Pick a guide file');
if(filename)
    handles.mutposFile = strcat(pathname, filename);
    set(handles.guideEdit, 'String', handles.mutposFile);
    setStatusLabel('Marks guide is set',handles);

    guidata(hObject, handles);
end


function shiftEdit_Callback(hObject, eventdata, handles)
% hObject    handle to shiftEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of shiftEdit as text
%        str2double(get(hObject,'String')) returns contents of shiftEdit as a double


% --- Executes during object creation, after setting all properties.
function shiftEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to shiftEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function windowEdit_Callback(hObject, eventdata, handles)
% hObject    handle to windowEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of windowEdit as text
%        str2double(get(hObject,'String')) returns contents of windowEdit as a double


% --- Executes during object creation, after setting all properties.
function windowEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to windowEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in dpCheck.
function dpCheck_Callback(hObject, eventdata, handles)
% hObject    handle to dpCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of dpCheck



function targetblockEdit_Callback(hObject, eventdata, handles)
% hObject    handle to targetblockEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of targetblockEdit as text
%        str2double(get(hObject,'String')) returns contents of targetblockEdit as a double


% --- Executes during object creation, after setting all properties.
function targetblockEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to targetblockEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function slackEdit_Callback(hObject, eventdata, handles)
% hObject    handle to slackEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of slackEdit as text
%        str2double(get(hObject,'String')) returns contents of slackEdit as a double


% --- Executes during object creation, after setting all properties.
function slackEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slackEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sequenceEdit_Callback(hObject, eventdata, handles)
% hObject    handle to sequenceEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sequenceEdit as text
%        str2double(get(hObject,'String')) returns contents of sequenceEdit as a double


% --- Executes during object creation, after setting all properties.
function sequenceEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sequenceEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in seqBtn.
function seqBtn_Callback(hObject, eventdata, handles)
% hObject    handle to seqBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
setStatusLabel('Open a sequence',handles);
[filename pathname]= uigetfile('*.txt','Pick a sequence file');
if(filename)
    set(handles.sequenceEdit, 'String', strcat(pathname,filename));
end

% --- Executes on button press in loadlistBtn.
function loadlistBtn_Callback(hObject, eventdata, handles)
% hObject    handle to loadlistBtn (see GCBO)
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
        
        h = handles.fileListbox;
        set(h, 'String', filenames{1},'Max', length(filenames{1}), 'Value', 1);
        
        guidata(hObject, handles);
    end
end

% --- Executes on selection change in fileListbox.
function fileListbox_Callback(hObject, eventdata, handles)
% hObject    handle to fileListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns fileListbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from fileListbox


% --- Executes during object creation, after setting all properties.
function fileListbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fileListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in addBtn.
function addBtn_Callback(hObject, eventdata, handles)
% hObject    handle to addBtn (see GCBO)
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
     set(handles.fileListbox,'String', str, 'Max', length(str), 'Value', 1);
     
     
     guidata(hObject, handles);
 end

% --- Executes on button press in removeBtn.
function removeBtn_Callback(hObject, eventdata, handles)
% hObject    handle to removeBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h = handles.fileListbox;
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


function refcolEdit_Callback(hObject, eventdata, handles)
% hObject    handle to refcolEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of refcolEdit as text
%        str2double(get(hObject,'String')) returns contents of refcolEdit as a double


% --- Executes during object creation, after setting all properties.
function refcolEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to refcolEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function yminEdit_Callback(hObject, eventdata, handles)
% hObject    handle to yminEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of yminEdit as text
%        str2double(get(hObject,'String')) returns contents of yminEdit as a double


% --- Executes during object creation, after setting all properties.
function yminEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to yminEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ymaxEdit_Callback(hObject, eventdata, handles)
% hObject    handle to ymaxEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ymaxEdit as text
%        str2double(get(hObject,'String')) returns contents of ymaxEdit as a double


% --- Executes during object creation, after setting all properties.
function ymaxEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ymaxEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function offsetEdit_Callback(hObject, eventdata, handles)
% hObject    handle to offsetEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of offsetEdit as text
%        str2double(get(hObject,'String')) returns contents of offsetEdit as a double


% --- Executes during object creation, after setting all properties.
function offsetEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to offsetEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in baselineCheck.
function baselineCheck_Callback(hObject, eventdata, handles)
% hObject    handle to baselineCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of baselineCheck


% --- Executes on button press in refineCheck.
function refineCheck_Callback(hObject, eventdata, handles)
% hObject    handle to refineCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of refineCheck



function distEdit_Callback(hObject, eventdata, handles)
% hObject    handle to distEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of distEdit as text
%        str2double(get(hObject,'String')) returns contents of distEdit as a double


% --- Executes during object creation, after setting all properties.
function distEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to distEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function refreshSetting(handles)
set(handles.refcolEdit, 'String', handles.settings.refcol);
set(handles.offsetEdit, 'String', handles.settings.offset);
set(handles.yminEdit, 'String', handles.settings.ymin);
set(handles.ymaxEdit, 'String', handles.settings.ymax);
set(handles.distEdit, 'String', handles.settings.dist);
set(handles.refineCheck, 'Value', handles.settings.refine);
set(handles.baselineCheck, 'Value', handles.settings.baseline);

set(handles.startEdit, 'String', handles.settings.eternaStart);
set(handles.endEdit, 'String', handles.settings.eternaEnd);
set(handles.eternaCheck, 'Value', handles.settings.eternaCheck);

set(handles.peakspacingEdit, 'String', handles.settings.peakspacing);

set(handles.slackEdit, 'String', handles.fineTune.slack);
set(handles.shiftEdit, 'String', handles.fineTune.shift);
set(handles.windowEdit, 'String', handles.fineTune.windowsize);
set(handles.targetblockEdit, 'String', handles.fineTune.targetblock);
set(handles.dpCheck, 'Value', handles.fineTune.dpAlign);

set(handles.fileListbox, 'String', handles.filenames, 'Value', 1);
set(handles.sequenceEdit, 'String', handles.sequenceFile);
set(handles.strEdit, 'String', handles.structure);
set(handles.guideEdit, 'String', handles.mutposFile);
set(handles.strEdit, 'String', handles.structure);

function handles = setupSetting(handles)
handles.settings.refcol = str2double(get(handles.refcolEdit, 'String'));
handles.settings.offset = str2double(get(handles.offsetEdit, 'String'));
handles.settings.ymin = str2double(get(handles.yminEdit, 'String'));
handles.settings.ymax = str2double(get(handles.ymaxEdit, 'String'));
handles.settings.dist = str2double(get(handles.distEdit, 'String'));
handles.settings.refine = get(handles.refineCheck, 'Value');
handles.settings.baseline = get(handles.baselineCheck, 'Value');

handles.fineTune.slack = str2double(get(handles.slackEdit, 'String'));
handles.fineTune.shift = str2double(get(handles.shiftEdit, 'String'));
handles.fineTune.windowsize = str2double(get(handles.windowEdit, 'String'));
handles.fineTune.targetblock = get(handles.targetblockEdit, 'String');
handles.fineTune.dpAlign = get(handles.dpCheck, 'Value');

handles.settings.eternaStart = str2double(get(handles.startEdit, 'String'));
handles.settings.eternaEnd = str2double(get(handles.startEdit, 'String'));
handles.settings.eternaCheck = get(handles.eternaCheck, 'Value');

handles.settings.peakspacing = str2double(get(handles.peakspacingEdit, 'String'));

if(handles.settings.eternaCheck)
    file = get(handles.sequenceEdit, 'String');
    [handles.ids, handles.target_names, handles.subrounds, handles.sequence, handles.design_names ] = parse_EteRNA_sequence_file( file );
    set(handles.profileCombo,'String', handles.design_names);
    set(handles.dataCombo,'String', {'1', '2'});
    handles.max = 6;
else
    file = get(handles.sequenceEdit, 'String');
    fid = fopen(file);
    str = textscan(fid, '%s');
    handles.sequence = str{1}{1};
    fclose(fid);
    handles.max = 5;
end

handles.filenames = get(handles.fileListbox, 'String');
handles.sequenceFile = get(handles.sequenceEdit, 'String');
handles.structure = get(handles.strEdit ,'String');
handles.mutposFile = get(handles.guideEdit, 'String');



function edit18_Callback(hObject, eventdata, handles)
% hObject    handle to edit18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit18 as text
%        str2double(get(hObject,'String')) returns contents of edit18 as a double


% --- Executes during object creation, after setting all properties.
function edit18_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function setStatusLabel(str,handles)
set(handles.statusLabel, 'String', sprintf('Status: %s',str));


% --- Executes on button press in prevBtn.
function prevBtn_Callback(hObject, eventdata, handles)
% hObject    handle to prevBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
    for j = handles.displayComponents
        delete(j);
    end
catch err
    close(handles.figure1);
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

% --- Executes on button press in nextBtn.
function nextBtn_Callback(hObject, eventdata, handles)
% hObject    handle to nextBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
    for j = handles.displayComponents
        delete(j);
    end
catch err
    close(handles.figure1);
end

if(handles.step == handles.max)
    warndlg('This is the last step!');
    handles.displayComponents = displaySetting(handles.stages{handles.step + 1},handles);
    menuSetting('lastpage',handles);
else
    handles.step = handles.step + 1;
    handles.displayComponents = displaySetting(handles.stages{handles.step + 1},handles);

    if(handles.step == handles.max)
        menuSetting('lastpage',handles);
    else
        menuSetting('midpage',handles);
    end
end

guidata(hObject,handles);



function strEdit_Callback(hObject, eventdata, handles)
% hObject    handle to strEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of strEdit as text
%        str2double(get(hObject,'String')) returns contents of strEdit as a double


% --- Executes during object creation, after setting all properties.
function strEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to strEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton23.
function pushbutton23_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename pathname]= uigetfile('*.txt','Pick a structure file');
if(filename)
    fid = fopen(strcat(pathname,filename));
    str = textscan(fid, '%s');
    fclose(fid);
    set(handles.strEdit, 'String', str{1}{1});
end


function startEdit_Callback(hObject, eventdata, handles)
% hObject    handle to startEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of startEdit as text
%        str2double(get(hObject,'String')) returns contents of startEdit as a double


% --- Executes during object creation, after setting all properties.
function startEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to startEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function endEdit_Callback(hObject, eventdata, handles)
% hObject    handle to endEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of endEdit as text
%        str2double(get(hObject,'String')) returns contents of endEdit as a double


% --- Executes during object creation, after setting all properties.
function endEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to endEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in eternaCheck.
function eternaCheck_Callback(hObject, eventdata, handles)
% hObject    handle to eternaCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
eterna = get(hObject, 'Value');
if(eterna)
    set(handles.startEdit, 'Enable', 'On');
    set(handles.endEdit, 'Enable', 'On');
else
    set(handles.startEdit, 'Enable', 'Off');
    set(handles.endEdit, 'Enable', 'Off');
end
    
% Hint: get(hObject,'Value') returns toggle state of eternaCheck


% --- Executes on selection change in profileCombo.
function profileCombo_Callback(hObject, eventdata, handles)
% hObject    handle to profileCombo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns profileCombo contents as cell array
%        contents{get(hObject,'Value')} returns selected item from profileCombo
for j = handles.displayComponents
    delete(j);
end
handles.displayComponents = displaySetting(handles.stages{handles.step + 1}, handles);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function profileCombo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to profileCombo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in dataCombo.
function dataCombo_Callback(hObject, eventdata, handles)
% hObject    handle to dataCombo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns dataCombo contents as cell array
%        contents{get(hObject,'Value')} returns selected item from dataCombo
for j = handles.displayComponents
    delete(j);
end
handles.displayComponents = displaySetting(handles.stages{handles.step + 1}, handles);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function dataCombo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dataCombo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in resetBtn.
function resetBtn_Callback(hObject, eventdata, handles)
% hObject    handle to resetBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(handles.figure1);
HiTRACE_figure;



function peakspacingEdit_Callback(hObject, eventdata, handles)
% hObject    handle to peakspacingEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of peakspacingEdit as text
%        str2double(get(hObject,'String')) returns contents of peakspacingEdit as a double


% --- Executes during object creation, after setting all properties.
function peakspacingEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to peakspacingEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
