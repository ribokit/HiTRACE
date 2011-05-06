function varargout = FineTunning(varargin)
% FINETUNNING MATLAB code for FineTunning.fig
%      FINETUNNING, by itself, creates a new FINETUNNING or raises the existing
%      singleton*.
%
%      H = FINETUNNING returns the handle to a new FINETUNNING or the handle to
%      the existing singleton*.
%
%      FINETUNNING('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FINETUNNING.M with the given input arguments.
%
%      FINETUNNING('Property','Value',...) creates a new FINETUNNING or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before FineTunning_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to FineTunning_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help FineTunning

% Last Modified by GUIDE v2.5 03-May-2011 15:22:27

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @FineTunning_OpeningFcn, ...
                   'gui_OutputFcn',  @FineTunning_OutputFcn, ...
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


% --- Executes just before FineTunning is made visible.
function FineTunning_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to FineTunning (see VARARGIN)

handles.settings = varargin{1};
initialization(handles.settings, handles);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GeneralOption wait for user response (see UIRESUME)
uiwait(hObject);


% --- Outputs from this function are returned to the command line.
function varargout = FineTunning_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.settings;
delete(hObject);

% --- Executes on button press in dpCheck.
function dpCheck_Callback(hObject, eventdata, handles)
% hObject    handle to dpCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 h1 = handles.slackText;
 h2 = handles.shiftText;
 h3 = handles.windowsizeText;
 h4 = handles.targetblockText;

 if(get(hObject,'Value'))
     set(h1, 'Enable', 'on');
     set(h2, 'Enable', 'on');
     set(h3, 'Enable', 'on');
     set(h4, 'Enable', 'on');
 else
     set(h1, 'Enable', 'off');
     set(h2, 'Enable', 'off');
     set(h3, 'Enable', 'off');
     set(h4, 'Enable', 'off');
 end



function targetblockText_Callback(hObject, eventdata, handles)
% hObject    handle to targetblockText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of targetblockText as text
%        str2double(get(hObject,'String')) returns contents of targetblockText as a double


% --- Executes during object creation, after setting all properties.
function targetblockText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to targetblockText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function shiftText_Callback(hObject, eventdata, handles)
% hObject    handle to shiftText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of shiftText as text
%        str2double(get(hObject,'String')) returns contents of shiftText as a double


% --- Executes during object creation, after setting all properties.
function shiftText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to shiftText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function windowsizeText_Callback(hObject, eventdata, handles)
% hObject    handle to windowsizeText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of windowsizeText as text
%        str2double(get(hObject,'String')) returns contents of windowsizeText as a double


% --- Executes during object creation, after setting all properties.
function windowsizeText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to windowsizeText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function slackText_Callback(hObject, eventdata, handles)
% hObject    handle to slackText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of slackText as text
%        str2double(get(hObject,'String')) returns contents of slackText as a double


% --- Executes during object creation, after setting all properties.
function slackText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slackText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in okButton.
function okButton_Callback(hObject, eventdata, handles)
% hObject    handle to okButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
settings.slack = str2double(get(handles.slackText, 'String'));
settings.shift = str2double(get(handles.shiftText, 'String'));
settings.windowsize = str2double(get(handles.windowsizeText, 'String'));
settings.targetblock = get(handles.targetblockText, 'String');
settings.dpAlign = get(handles.dpCheck, 'Value');


% input error check
if(sum(isnan([settings.slack settings.shift settings.windowsize settings.dpAlign])))
    errordlg('Input values must be number!', 'Error');
    return;
end

% save setting status
handles.settings = settings;
guidata(hObject, handles);
uiresume(handles.figure1);

% --- Executes on button press in cancelButton.
function cancelButton_Callback(hObject, eventdata, handles)
% hObject    handle to cancelButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiresume(handles.figure1);


function initialization(settings, handles)
if(isfield(settings, 'slack'))
    set(handles.slackText, 'String', int2str(settings.slack));
end

if(isfield(settings, 'shift'))
    set(handles.shiftText, 'String', int2str(settings.shift));
end

if(isfield(settings, 'windowsize'))
    set(handles.windowsizeText, 'String', int2str(settings.windowsize));
end

if(isfield(settings, 'dpAlign'))
    h1 = handles.slackText;
    h2 = handles.shiftText;
    h3 = handles.windowsizeText;
    h4 = handles.targetblockText;
    
    if(settings.dpAlign)
        set(h1, 'Enable', 'on');
        set(h2, 'Enable', 'on');
        set(h3, 'Enable', 'on');
        set(h4, 'Enable', 'on');
    else
        set(h1, 'Enable', 'off');
        set(h2, 'Enable', 'off');
        set(h3, 'Enable', 'off');
        set(h4, 'Enable', 'off');
    end
    
    set(handles.dpCheck, 'Value', settings.dpAlign);
end

if(isfield(settings, 'targetblock'))
    set(handles.targetblockText, 'String', settings.targetblock);
end


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiresume(hObject);
