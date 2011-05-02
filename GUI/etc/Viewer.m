function varargout = Viewer(varargin)
% VIEWER M-file for Viewer.fig
%      VIEWER, by itself, creates a new VIEWER or raises the existing
%      singleton*.
%
%      H = VIEWER returns the handle to a new VIEWER or the handle to
%      the existing singleton*.
%
%      VIEWER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VIEWER.M with the given input arguments.
%
%      VIEWER('Property','Value',...) creates a new VIEWER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Viewer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Viewer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Viewer

% Last Modified by GUIDE v2.5 28-Mar-2011 20:40:15

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Viewer_OpeningFcn, ...
                   'gui_OutputFcn',  @Viewer_OutputFcn, ...
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


% --- Executes just before Viewer is made visible.
function Viewer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Viewer (see VARARGIN)

axis(handles.axes3, 'off');
[img map] = imread('loading48.gif', 'frame', 'all');


% Choose default command line output for Viewer
handles.output = hObject;
handles.step = 0;
handles.loadingImage = img;
handles.loadingImageColor = map;
guidata(hObject, handles);

global loadingTimer;
loadingTimer = 1;

global t;
t = timer('TimerFcn',{@loading_timer, handles}, 'ExecutionMode', 'fixedRate' ,'Period', 0.13);


start(t);


% Update handles structure

% UIWAIT makes Viewer wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Viewer_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles;

% --- Executes on button press in prevBtn.
function prevBtn_Callback(hObject, eventdata, handles)
% hObject    handle to prevBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if(handles.step < 2)
    warndlg('This is the first step!');
else
    handles.step = handles.step - 1;
    guidata(handles.output, handles);
    
    deleteImage(handles.axesHandles);
    createImage(handles.step, handles);
end



% --- Executes on button press in restartBtn.
function restartBtn_Callback(hObject, eventdata, handles)
% hObject    handle to restartBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global stepNo;
stepNo = handles.step - 1;

deleteImage(handles.axesHandles);
set(handles.prevBtn, 'Enable', 'Off');
set(handles.nextBtn, 'Enable', 'Off');
set(handles.restartBtn, 'Enable', 'Off');
uiresume;

% --- Executes on button press in nextBtn.
function nextBtn_Callback(hObject, eventdata, handles)
% hObject    handle to nextBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if(handles.step > 4)
    warndlg('This is the last step!');
else
    handles.step = handles.step + 1;
    guidata(handles.output, handles);
    
    deleteImage(handles.axesHandles);
    createImage(handles.step, handles);
end

function deleteImage(handles)
for i = handles
    delete(i);
end
        

function createImage(currentStep, handles)
global SETTING;
switch(currentStep)
    case 1
        h1 = axes('position', [0.1, 0.2, 0.35, 0.7]);
        image( 50 * handles.d);
        axis( [ 0.5 size(handles.da,2)+0.5 SETTING.ymin SETTING.ymax] );
        title( 'Signal')
        
        h2 = axes('position', [0.55, 0.2, 0.35, 0.7]);
        image( 50 * handles.da);
        axis( [ 0.5 size(handles.da,2)+0.5 SETTING.ymin SETTING.ymax] );
        title( 'Reference ladder')
        
        colormap( 1- gray(100));
        
        set(handles.stepLabel, 'String', 'Loading and alignment (step 1/5)');
        
        handles.axesHandles = [h1, h2];        
        guidata(handles.output, handles);
    case 2
        if(SETTING.enableBaseline)
            str = 'Signal (baseline subtraction enabled)';
        else
            str = 'Signal (baseline subtraction skipped)';
        end
        
        h1 = axes('Position', [0.1, 0.2, 0.8, 0.7]);
        
        image( 50 * handles.d_bsub);
        axis( [ 0.5 size(handles.d_bsub,2)+0.5 SETTING.ymin SETTING.ymax] );        
        title( str );
        
        colormap( 1- gray(100));
        
        set(handles.stepLabel, 'String', 'Baseline subtraction (step 2/5)');
        
        handles.axesHandles = h1;        
        guidata(handles.output, handles);
    case 3
        if SETTING.enableDP
            str = 'Signal (nonlinear alignment by DP)';
        else
            str = 'Signal (nonlinear alignment by DP skipped)';
        end
        
        h1 = axes('Position', [0.1, 0.2, 0.8, 0.7]);
        
        image( 50 * handles.d_align);
        axis( [ 0.5 size(handles.d_bsub,2)+0.5 1 size( handles.d_align,1)] )
        title( str );
        
        colormap( 1- gray(100));
        
        set(handles.stepLabel, 'String', 'Nonlinear alignment by DP (step 3/5)');
        handles.axesHandles = h1;        
        guidata(handles.output, handles);
    case 4
        
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
		  which_xsel = handles.numpeaks - handles.marks( k, 2 ) + 1 + SETTING.offset;
		  if ( which_xsel >= 1 & which_xsel <= length(handles.xsel) )
		    plot( m, handles.xsel(which_xsel),'ro' );
		  end
		end
	      end
            end
        end
        hold off      
        
        set(handles.stepLabel, 'String', 'Manual annotation (step 4/5)');
        handles.axesHandles = h1;        
        guidata(handles.output, handles);
    case 5
        if(SETTING.enableRefine)
            str = 'Fitting without refinement (step 5/5)';
        else
            str = 'Fitting with refinement (step 5/5)';
        end
        
        h1 = axes('Position', [0.08, 0.2, 0.27, 0.7]);
        scalefactor = 40/mean(mean(handles.d_align));
        image( scalefactor * handles.d_align );
        if length( handles.xsel_start ) > 0; ylim( [ min(handles.xsel_start)-100, max(handles.xsel_start)+100 ] ); end
        title( 'data' );

        h2 = axes('Position', [0.38, 0.2, 0.27, 0.7]);
        image( scalefactor * handles.prof_fit );
        if length( handles.xsel_start ) > 0; ylim( [ min(handles.xsel_start)-100, max(handles.xsel_start)+100 ] ); end
        title( 'fit' );

        h3 = axes('Position', [0.68, 0.2, 0.27, 0.7]);
        image( scalefactor * ( handles.d_align - handles.prof_fit) );
        if length( handles.xsel_start ) > 0; ylim( [ min(handles.xsel_start)-100, max(handles.xsel_start)+100 ] ); end;
        title( 'residuals' );
        
        colormap( 1 - gray(100) );
        
        set(handles.stepLabel, 'String', str);
        handles.axesHandles = [h1 h2 h3];        
        guidata(handles.output, handles);
end


% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global stepNo;
global setter;
stepNo = -1;
set(setter.startBtn, 'String', 'RUN');


function loading_timer(obj, event, handles)
global loadingTimer;
global t;
if loadingTimer == 0
    stop(t);
else
    loadingTimer = mod(loadingTimer,size(handles.loadingImage,4)) + 1;

    axes(handles.axes3);
        
    image(handles.loadingImage(401:end,:,1,loadingTimer));
    colormap(handles.loadingImageColor);
    axis off;
end


% --- Executes when figure1 is resized.
function figure1_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pos = get(hObject, 'Position');
posPrev = get(handles.prevBtn, 'Position');
posNext = get(handles.nextBtn, 'Position');
posRestart = get(handles.restartBtn, 'Position');

posPrev(3) = 129 / pos(3);
posPrev(4) = 30 / pos(4);

posNext(3) = 101 / pos(3);
posNext(4) = 30 / pos(4);

posRestart(3) = 166 / pos(3);
posRestart(4) = 30 / pos(4);


set(handles.prevBtn, 'Position', posPrev);
set(handles.nextBtn, 'Position', posNext);
set(handles.restartBtn, 'Position', posRestart);


% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
