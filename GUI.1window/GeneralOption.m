function varargout = GeneralOption(varargin)
% GENERALOPTION MATLAB code for GeneralOption.fig
%      GENERALOPTION, by itself, creates a new GENERALOPTION or raises the existing
%      singleton*.
%
%      H = GENERALOPTION returns the handle to a new GENERALOPTION or the handle to
%      the existing singleton*.
%
%      GENERALOPTION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GENERALOPTION.M with the given input arguments.
%
%      GENERALOPTION('Property','Value',...) creates a new GENERALOPTION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GeneralOption_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GeneralOption_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GeneralOption

% Last Modified by GUIDE v2.5 03-May-2011 14:27:53

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GeneralOption_OpeningFcn, ...
                   'gui_OutputFcn',  @GeneralOption_OutputFcn, ...
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


% --- Executes just before GeneralOption is made visible.
function GeneralOption_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GeneralOption (see VARARGIN)

handles.settings = varargin{1};
initialization(handles.settings, handles);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GeneralOption wait for user response (see UIRESUME)
uiwait(hObject);


% --- Outputs from this function are returned to the command line.
function varargout = GeneralOption_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.settings;
delete(hObject);


function refcolText_Callback(hObject, eventdata, handles)
% hObject    handle to refcolText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of refcolText as text
%        str2double(get(hObject,'String')) returns contents of refcolText as a double


% --- Executes during object creation, after setting all properties.
function refcolText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to refcolText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function yminText_Callback(hObject, eventdata, handles)
% hObject    handle to yminText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of yminText as text
%        str2double(get(hObject,'String')) returns contents of yminText as a double


% --- Executes during object creation, after setting all properties.
function yminText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to yminText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ymaxText_Callback(hObject, eventdata, handles)
% hObject    handle to ymaxText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ymaxText as text
%        str2double(get(hObject,'String')) returns contents of ymaxText as a double


% --- Executes during object creation, after setting all properties.
function ymaxText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ymaxText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function offsetText_Callback(hObject, eventdata, handles)
% hObject    handle to offsetText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of offsetText as text
%        str2double(get(hObject,'String')) returns contents of offsetText as a double


% --- Executes during object creation, after setting all properties.
function offsetText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to offsetText (see GCBO)
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



function distText_Callback(hObject, eventdata, handles)
% hObject    handle to distText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of distText as text
%        str2double(get(hObject,'String')) returns contents of distText as a double


% --- Executes during object creation, after setting all properties.
function distText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to distText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



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


% --- Executes on button press in okButton.
function okButton_Callback(hObject, eventdata, handles)
% hObject    handle to okButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
settings.refcol = str2double(get(handles.refcolText, 'String'));
settings.offset = str2double(get(handles.offsetText, 'String'));
settings.ymin = str2double(get(handles.yminText, 'String'));
settings.ymax = str2double(get(handles.ymaxText, 'String'));
settings.dist = str2double(get(handles.distText, 'String'));
settings.refine= get(handles.refineCheck, 'Value');
settings.baseline = get(handles.baselineCheck, 'Value');


% input error check
if(sum(isnan([settings.refcol settings.offset settings.ymin settings.ymax settings.dist settings.refine settings.baseline])))
    errordlg('An input value must be number!', 'Error');
    return;
end

if(settings.ymin > settings.ymax)
    errordlg('Ymin value cannot exceed Ymax value!', 'Error');
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

if(isfield(settings, 'refcol'))
    set(handles.refcolText, 'String', int2str(settings.refcol));
end

if(isfield(settings, 'offset'))
    set(handles.offsetText, 'String', int2str(settings.offset));
end

if(isfield(settings, 'ymin'))
    set(handles.yminText, 'String', int2str(settings.ymin));
end

if(isfield(settings, 'ymax'))
    set(handles.ymaxText, 'String', int2str(settings.ymax));
end

if(isfield(settings, 'dist'))
    set(handles.distText, 'String', int2str(settings.dist));
end

if(isfield(settings, 'refine'))
    set(handles.refineCheck, 'Value', settings.refine);
end

if(isfield(settings, 'baseline'))
    set(handles.baselineCheck, 'Value', settings.baseline);
end


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiresume(handles.figure1);
