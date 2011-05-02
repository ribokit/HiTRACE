function varargout = SettingWindow(varargin)
% SettingWindow M-file for SettingWindow.fig
%      SettingWindow, by itself, creates a new SettingWindow or raises the existing
%      singleton*.
%
%      H = SettingWindow returns the handle to a new SettingWindow or the handle to
%      the existing singleton*.
%
%      SettingWindow('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SettingWindow.M with the given input arguments.
%
%      SettingWindow('Property','Value',...) creates a new SettingWindow or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SettingWindow_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SettingWindow_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

 % Edit the above text to modify the response to help SettingWindow

 % Last Modified by GUIDE v2.5 07-Apr-2011 13:08:02

 % Begin initialization code - DO NOT EDIT
 gui_Singleton = 1;
 gui_State = struct('gui_Name',       mfilename, ...
		    'gui_Singleton',  gui_Singleton, ...
		    'gui_OpeningFcn', @SettingWindow_OpeningFcn, ...
		    'gui_OutputFcn',  @SettingWindow_OutputFcn, ...
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


 % --- Executes just before SettingWindow is made visible.
 function SettingWindow_OpeningFcn(hObject, eventdata, handles, varargin)
 % This function has no output args, see OutputFcn.
 % hObject    handle to figure
 % eventdata  reserved - to be defined in a future version of MATLAB
 % handles    structure with handles and user data (see GUIDATA)
 % varargin   command line arguments to SettingWindow (see VARARGIN)

 % Choose default command line output for SettingWindow
 handles.output = hObject;

 %axes(handles.imageAxes);
 %imshow('HiTRACE.png');

 axes(handles.logo);
 imshow('hitrace_logo.png');


 global res;
 str = sprintf('%d thread(s) available',res.ClusterSize);
 %set(handles.recommendTxt, 'String', str);
 set(handles.numOfThreads, 'String', res.ClusterSize);
 % Update handles structure
 guidata(hObject, handles);

 % UIWAIT makes SettingWindow wait for user response (see UIRESUME)
 % uiwait(handles.figure1);


 % --- Outputs from this function are returned to the command line.
 function varargout = SettingWindow_OutputFcn(hObject, eventdata, handles) 
 % varargout  cell array for returning output args (see VARARGOUT);
 % hObject    handle to figure
 % eventdata  reserved - to be defined in a future version of MATLAB
 % handles    structure with handles and user data (see GUIDATA)

 % Get default command line output from handles structure
 varargout{1} = handles;


 % --- Executes on button press in startBtn.
 function startBtn_Callback(hObject, eventdata, handles)
 % hObject    handle to startBtn (see GCBO)
 % eventdata  reserved - to be defined in a future version of MATLAB
 % handles    structure with handles and user data (see GUIDATA)

 global SETTING;
 global res;

 handles_to_SETTINGS( handles )


 if(isempty(SETTING.list))
     errordlg('Please add one or more data folders. You can select the location of a folder by clicking the ''Add'' button. Alternatively, you can load a text file that contains a list of folders (one path per line).','Error!');
     return;
 end

 if(isempty(SETTING.sequence))
     errordlg('Please specify your input sequence. You can either type a sequence string directly in the window or load one from a plain text file.','Error!');
     return;
 end

 %if(get(hParallel, 'Value'))

 % Parallelize by default. 
 % Currently, internal scripts don't actually use res.ClusterSize do determine whether 
 % to use parfor or for, do they? Can fix this later...
 
 hNumPar = handles.numOfThreads;

 if( true )
     nums = str2double(get(hNumPar, 'String'));
     if(nums > res.ClusterSize)
	 str = sprintf('You specified more threads (%d) than currently available (%d). Only %d threads will be activated.', nums, res.ClusterSize, res.ClusterSize);
	 reply = questdlg(str, 'Confirm!', 'OK', 'Cancel', 'OK');
	 if(strcmp( reply,'OK'))
	     set(hNumPar, 'String', res.ClusterSize);
	     SETTING.numOfJob = res.ClusterSize;
	     uiresume;
	 end
     else
	 SETTING.numOfJob = nums;
     end
 else
     SETTING.numOfJob = 0;

 end

 str = get(hObject, 'String');
 if(strcmp(str, 'RUN'))
     set(hObject, 'String', 'Update parameters');
     uiresume;
 else
     helpdlg('Your parameter change has been saved. Click "Restart from this step" button in the viewer window to rerun analysis');
 end

 function listFileName_Callback(hObject, eventdata, handles)
 % hObject    handle to listFileName (see GCBO)
 % eventdata  reserved - to be defined in a future version of MATLAB
 % handles    structure with handles and user data (see GUIDATA)

 % Hints: get(hObject,'String') returns contents of listFileName as text
 %        str2double(get(hObject,'String')) returns contents of listFileName as a double


 % --- Executes during object creation, after settingwindow all properties.
 function listFileName_CreateFcn(hObject, eventdata, handles)
 % hObject    handle to listFileName (see GCBO)
 % eventdata  reserved - to be defined in a future version of MATLAB
 % handles    empty - handles not created until after all CreateFcns called

 % Hint: edit controls usually have a white background on Windows.
 %       See ISPC and COMPUTER.
 if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
     set(hObject,'BackgroundColor','white');
 end


 % --- Executes on button press in openListDialog.
 function openListDialog_Callback(hObject, eventdata, handles)
 % hObject    handle to openListDialog (see GCBO)
 % eventdata  reserved - to be defined in a future version of MATLAB
 % handles    structure with handles and user data (see GUIDATA)
 [filename pathname]= uigetfile('*.txt','Pick a list txt file');
 if(filename)
     fid = fopen(strcat(pathname,filename), 'r');
     if(fid == -1)
	 errordlg('File not found');
     else
	 filenames = textscan(fid,'%s');
	 h = handles.abiListbox;
	 set(h, 'String', filenames{1},'Max', length(filenames{1}), 'Value', 1);
     end
 end

 % --- If Enable == 'on', executes on mouse press in 5 pixel border.
 % --- Otherwise, executes on mouse press in 5 pixel border or over openListDialog.
 function openListDialog_ButtonDownFcn(hObject, eventdata, handles)
 % hObject    handle to openListDialog (see GCBO)
 % eventdata  reserved - to be defined in a future version of MATLAB
 % handles    structure with handles and user data (see GUIDATA)



 function ymin_Callback(hObject, eventdata, handles)
 % hObject    handle to ymin (see GCBO)
 % eventdata  reserved - to be defined in a future version of MATLAB
 % handles    structure with handles and user data (see GUIDATA)

 % Hints: get(hObject,'String') returns contents of ymin as text
 %        str2double(get(hObject,'String')) returns contents of ymin as a double


 % --- Executes during object creation, after settingwindow all properties.
 function ymin_CreateFcn(hObject, eventdata, handles)
 % hObject    handle to ymin (see GCBO)
 % eventdata  reserved - to be defined in a future version of MATLAB
 % handles    empty - handles not created until after all CreateFcns called

 % Hint: edit controls usually have a white background on Windows.
 %       See ISPC and COMPUTER.
 if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
     set(hObject,'BackgroundColor','white');
 end



 function ymax_Callback(hObject, eventdata, handles)
 % hObject    handle to ymax (see GCBO)
 % eventdata  reserved - to be defined in a future version of MATLAB
 % handles    structure with handles and user data (see GUIDATA)

 % Hints: get(hObject,'String') returns contents of ymax as text
 %        str2double(get(hObject,'String')) returns contents of ymax as a double


 % --- Executes during object creation, after settingwindow all properties.
 function ymax_CreateFcn(hObject, eventdata, handles)
 % hObject    handle to ymax (see GCBO)
 % eventdata  reserved - to be defined in a future version of MATLAB
 % handles    empty - handles not created until after all CreateFcns called

 % Hint: edit controls usually have a white background on Windows.
 %       See ISPC and COMPUTER.
 if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
     set(hObject,'BackgroundColor','white');
 end


 % --- Executes on button press in checkbox3.
 function checkbox3_Callback(hObject, eventdata, handles)
 % hObject    handle to checkbox3 (see GCBO)
 % eventdata  reserved - to be defined in a future version of MATLAB
 % handles    structure with handles and user data (see GUIDATA)

 % Hint: get(hObject,'Value') returns toggle state of checkbox3


 % --- Executes on button press in parallelCheck.
 function parallelCheck_Callback(hObject, eventdata, handles)
 % hObject    handle to parallelCheck (see GCBO)
 % eventdata  reserved - to be defined in a future version of MATLAB
 % handles    structure with handles and user data (see GUIDATA)

 % Hint: get(hObject,'Value') returns toggle state of parallelCheck

 h = handles.numOfThreads;
 if(get(hObject,'Value'))
     set(h, 'Enable', 'on');
 else
     set(h, 'Enable', 'off');
 end


 function numOfThreads_Callback(hObject, eventdata, handles)
 % hObject    handle to numOfThreads (see GCBO)
 % eventdata  reserved - to be defined in a future version of MATLAB
 % handles    structure with handles and user data (see GUIDATA)

 % Hints: get(hObject,'String') returns contents of numOfThreads as text
 %        str2double(get(hObject,'String')) returns contents of numOfThreads as a double


 % --- Executes during object creation, after settingwindow all properties.
 function numOfThreads_CreateFcn(hObject, eventdata, handles)
 % hObject    handle to numOfThreads (see GCBO)
 % eventdata  reserved - to be defined in a future version of MATLAB
 % handles    empty - handles not created until after all CreateFcns called

 % Hint: edit controls usually have a white background on Windows.
 %       See ISPC and COMPUTER.
 if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
     set(hObject,'BackgroundColor','white');
 end



 function refCol_Callback(hObject, eventdata, handles)
 % hObject    handle to refCol (see GCBO)
 % eventdata  reserved - to be defined in a future version of MATLAB
 % handles    structure with handles and user data (see GUIDATA)

 % Hints: get(hObject,'String') returns contents of refCol as text
 %        str2double(get(hObject,'String')) returns contents of refCol as a double


 % --- Executes during object creation, after settingwindow all properties.
 function refCol_CreateFcn(hObject, eventdata, handles)
 % hObject    handle to refCol (see GCBO)
 % eventdata  reserved - to be defined in a future version of MATLAB
 % handles    empty - handles not created until after all CreateFcns called

 % Hint: edit controls usually have a white background on Windows.
 %       See ISPC and COMPUTER.
 if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
     set(hObject,'BackgroundColor','white');
 end


 % --- Executes on button press in exitBtn.
 function exitBtn_Callback(hObject, eventdata, handles)
 % hObject    handle to exitBtn (see GCBO)
 % eventdata  reserved - to be defined in a future version of MATLAB
 % handles    structure with handles and user data (see GUIDATA)
 quit;



 function sequence_Callback(hObject, eventdata, handles)
 % hObject    handle to sequence (see GCBO)
 % eventdata  reserved - to be defined in a future version of MATLAB
 % handles    structure with handles and user data (see GUIDATA)

 % Hints: get(hObject,'String') returns contents of sequence as text
 %        str2double(get(hObject,'String')) returns contents of sequence as a double


 % --- Executes during object creation, after settingwindow all properties.
 function sequence_CreateFcn(hObject, eventdata, handles)
 % hObject    handle to sequence (see GCBO)
 % eventdata  reserved - to be defined in a future version of MATLAB
 % handles    empty - handles not created until after all CreateFcns called

 % Hint: edit controls usually have a white background on Windows.
 %       See ISPC and COMPUTER.
 if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
     set(hObject,'BackgroundColor','white');
 end



 function offset_Callback(hObject, eventdata, handles)
 % hObject    handle to offset (see GCBO)
 % eventdata  reserved - to be defined in a future version of MATLAB
 % handles    structure with handles and user data (see GUIDATA)

 % Hints: get(hObject,'String') returns contents of offset as text
 %        str2double(get(hObject,'String')) returns contents of offset as a double


 % --- Executes during object creation, after settingwindow all properties.
 function offset_CreateFcn(hObject, eventdata, handles)
 % hObject    handle to offset (see GCBO)
 % eventdata  reserved - to be defined in a future version of MATLAB
 % handles    empty - handles not created until after all CreateFcns called

 % Hint: edit controls usually have a white background on Windows.
 %       See ISPC and COMPUTER.
 if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
     set(hObject,'BackgroundColor','white');
 end



 function parConfig_Callback(hObject, eventdata, handles)
 % hObject    handle to parConfig (see GCBO)
 % eventdata  reserved - to be defined in a future version of MATLAB
 % handles    structure with handles and user data (see GUIDATA)

 % Hints: get(hObject,'String') returns contents of parConfig as text
 %        str2double(get(hObject,'String')) returns contents of parConfig as a double


 % --- Executes during object creation, after settingwindow all properties.
 function parConfig_CreateFcn(hObject, eventdata, handles)
 % hObject    handle to parConfig (see GCBO)
 % eventdata  reserved - to be defined in a future version of MATLAB
 % handles    empty - handles not created until after all CreateFcns called

 % Hint: edit controls usually have a white background on Windows.
 %       See ISPC and COMPUTER.
 if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
     set(hObject,'BackgroundColor','white');
 end


 % --- Executes on button press in parConfigBtn.
 function parConfigBtn_Callback(hObject, eventdata, handles)
 % hObject    handle to parConfigBtn (see GCBO)
 % eventdata  reserved - to be defined in a future version of MATLAB
 % handles    structure with handles and user data (see GUIDATA)
 [filename pathname]= uigetfile('*.mat','Pick a parallel configuraion mat file');
 if(filename)
     h = handles.parConfig;
     set(h, 'String', strcat(pathname, filename));
 end


 % --- Executes on selection change in abiListbox.
 function abiListbox_Callback(hObject, eventdata, handles)
 % hObject    handle to abiListbox (see GCBO)
 % eventdata  reserved - to be defined in a future version of MATLAB
 % handles    structure with handles and user data (see GUIDATA)

 % Hints: contents = cellstr(get(hObject,'String')) returns abiListbox contents as cell array
 %        contents{get(hObject,'Value')} returns selected item from
 %        abiListbox


 % --- Executes during object creation, after settingwindow all properties.
 function abiListbox_CreateFcn(hObject, eventdata, handles)
 % hObject    handle to abiListbox (see GCBO)
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
     h = handles.abiListbox;
     str = get(h,'String');
     if(isempty(str))
	 str = {dir}; 
     else
	 str = cat(1, str, {dir});
     end
     set(h,'String', str, 'Max', length(str), 'Value', 1);
 end


 % --- Executes on button press in removeBtn.
 function removeBtn_Callback(hObject, eventdata, handles)
 % hObject    handle to removeBtn (see GCBO)
 % eventdata  reserved - to be defined in a future version of MATLAB
 % handles    structure with handles and user data (see GUIDATA)
 h = handles.abiListbox;
 names = get(h,'String');
 remove = get(h,'Value');

 ind = zeros(length(names),1);

 if remove > 0 & remove <= length(ind)
   ind(remove) = 1;
   ind = logical(ind);
   names = names(~ind);

   set(h, 'String', names, 'Max', length(names));

   if(length(names) < remove)
     set(h, 'Value', length(names));
   end
 end


 function slack_Callback(hObject, eventdata, handles)
 % hObject    handle to slack (see GCBO)
 % eventdata  reserved - to be defined in a future version of MATLAB
 % handles    structure with handles and user data (see GUIDATA)

 % Hints: get(hObject,'String') returns contents of slack as text
 %        str2double(get(hObject,'String')) returns contents of slack as a double


 % --- Executes during object creation, after settingwindow all properties.
 function slack_CreateFcn(hObject, eventdata, handles)
 % hObject    handle to slack (see GCBO)
 % eventdata  reserved - to be defined in a future version of MATLAB
 % handles    empty - handles not created until after all CreateFcns called

 % Hint: edit controls usually have a white background on Windows.
 %       See ISPC and COMPUTER.
 if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
     set(hObject,'BackgroundColor','white');
 end



 function maxShift_Callback(hObject, eventdata, handles)
 % hObject    handle to maxShift (see GCBO)
 % eventdata  reserved - to be defined in a future version of MATLAB
 % handles    structure with handles and user data (see GUIDATA)

 % Hints: get(hObject,'String') returns contents of maxShift as text
 %        str2double(get(hObject,'String')) returns contents of maxShift as a double


 % --- Executes during object creation, after settingwindow all properties.
 function maxShift_CreateFcn(hObject, eventdata, handles)
 % hObject    handle to maxShift (see GCBO)
 % eventdata  reserved - to be defined in a future version of MATLAB
 % handles    empty - handles not created until after all CreateFcns called

 % Hint: edit controls usually have a white background on Windows.
 %       See ISPC and COMPUTER.
 if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
     set(hObject,'BackgroundColor','white');
 end



 function windowSize_Callback(hObject, eventdata, handles)
 % hObject    handle to windowSize (see GCBO)
 % eventdata  reserved - to be defined in a future version of MATLAB
 % handles    structure with handles and user data (see GUIDATA)

 % Hints: get(hObject,'String') returns contents of windowSize as text
 %        str2double(get(hObject,'String')) returns contents of windowSize as a double


 % --- Executes during object creation, after settingwindow all properties.
 function windowSize_CreateFcn(hObject, eventdata, handles)
 % hObject    handle to windowSize (see GCBO)
 % eventdata  reserved - to be defined in a future version of MATLAB
 % handles    empty - handles not created until after all CreateFcns called

 % Hint: edit controls usually have a white background on Windows.
 %       See ISPC and COMPUTER.
 if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
     set(hObject,'BackgroundColor','white');
 end


 % --- Executes on button press in loadSeqBtn.
 function loadSeqBtn_Callback(hObject, eventdata, handles)
 % hObject    handle to loadSeqBtn (see GCBO)
 % eventdata  reserved - to be defined in a future version of MATLAB
 % handles    structure with handles and user data (see GUIDATA)
 [filename pathname]= uigetfile('*.txt','Pick a sequence file');
 if(filename)
     h = handles.sequence;
     a = textread(strcat(pathname,filename), '%s');
     set(h, 'String', a{1});
 end

 % --- Executes when user attempts to close figure1.
 function figure1_CloseRequestFcn(hObject, eventdata, handles)
 % hObject    handle to figure1 (see GCBO)
 % eventdata  reserved - to be defined in a future version of MATLAB
 % handles    structure with handles and user data (see GUIDATA)
 global SETTING;
 SETTING = [];
 % Hint: delete(hObject) closes the figure
 delete(hObject);


 % --- Executes on button press in enableDP.
 function enableDP_Callback(hObject, eventdata, handles)
 % hObject    handle to enableDP (see GCBO)
 % eventdata  reserved - to be defined in a future version of MATLAB
 % handles    structure with handles and user data (see GUIDATA)

 h1 = handles.slack;
 h2 = handles.maxShift;
 h3 = handles.windowSize;
 h4 = handles.targetBlock;

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

 % Hint: get(hObject,'Value') returns toggle state of enableDP



 function targetBlock_Callback(hObject, eventdata, handles)
 % hObject    handle to targetBlock (see GCBO)
 % eventdata  reserved - to be defined in a future version of MATLAB
 % handles    structure with handles and user data (see GUIDATA)

 % Hints: get(hObject,'String') returns contents of targetBlock as text
 %        str2double(get(hObject,'String')) returns contents of targetBlock as a double


 % --- Executes during object creation, after settingwindow all properties.
 function targetBlock_CreateFcn(hObject, eventdata, handles)
 % hObject    handle to targetBlock (see GCBO)
 % eventdata  reserved - to be defined in a future version of MATLAB
 % handles    empty - handles not created until after all CreateFcns called

 % Hint: edit controls usually have a white background on Windows.
 %       See ISPC and COMPUTER.
 if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
     set(hObject,'BackgroundColor','white');
 end


 % --- Executes on button press in applyBtn.
 function applyBtn_Callback(hObject, eventdata, handles)
 % hObject    handle to applyBtn (see GCBO)
 % eventdata  reserved - to be defined in a future version of MATLAB
 % handles    structure with handles and user data (see GUIDATA)


 % --- Executes on button press in enableBaseline.
 function enableBaseline_Callback(hObject, eventdata, handles)
 % hObject    handle to enableBaseline (see GCBO)
 % eventdata  reserved - to be defined in a future version of MATLAB
 % handles    structure with handles and user data (see GUIDATA)

 % Hint: get(hObject,'Value') returns toggle state of enableBaseline



 function edit14_Callback(hObject, eventdata, handles)
 % hObject    handle to edit14 (see GCBO)
 % eventdata  reserved - to be defined in a future version of MATLAB
 % handles    structure with handles and user data (see GUIDATA)

 % Hints: get(hObject,'String') returns contents of edit14 as text
 %        str2double(get(hObject,'String')) returns contents of edit14 as a double


 % --- Executes during object creation, after setting all properties.
 function edit14_CreateFcn(hObject, eventdata, handles)
 % hObject    handle to edit14 (see GCBO)
 % eventdata  reserved - to be defined in a future version of MATLAB
 % handles    empty - handles not created until after all CreateFcns called

 % Hint: edit controls usually have a white background on Windows.
 %       See ISPC and COMPUTER.
 if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
     set(hObject,'BackgroundColor','white');
 end



 function edit15_Callback(hObject, eventdata, handles)
 % hObject    handle to edit15 (see GCBO)
 % eventdata  reserved - to be defined in a future version of MATLAB
 % handles    structure with handles and user data (see GUIDATA)

 % Hints: get(hObject,'String') returns contents of edit15 as text
 %        str2double(get(hObject,'String')) returns contents of edit15 as a double


 % --- Executes during object creation, after setting all properties.
 function edit15_CreateFcn(hObject, eventdata, handles)
 % hObject    handle to edit15 (see GCBO)
 % eventdata  reserved - to be defined in a future version of MATLAB
 % handles    empty - handles not created until after all CreateFcns called

 % Hint: edit controls usually have a white background on Windows.
 %       See ISPC and COMPUTER.
 if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
     set(hObject,'BackgroundColor','white');
 end



 function edit16_Callback(hObject, eventdata, handles)
 % hObject    handle to edit16 (see GCBO)
 % eventdata  reserved - to be defined in a future version of MATLAB
 % handles    structure with handles and user data (see GUIDATA)

 % Hints: get(hObject,'String') returns contents of edit16 as text
 %        str2double(get(hObject,'String')) returns contents of edit16 as a double


 % --- Executes during object creation, after setting all properties.
 function edit16_CreateFcn(hObject, eventdata, handles)
 % hObject    handle to edit16 (see GCBO)
 % eventdata  reserved - to be defined in a future version of MATLAB
 % handles    empty - handles not created until after all CreateFcns called

 % Hint: edit controls usually have a white background on Windows.
 %       See ISPC and COMPUTER.
 if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
     set(hObject,'BackgroundColor','white');
 end


 % --- Executes on button press in checkbox7.
 function checkbox7_Callback(hObject, eventdata, handles)
 % hObject    handle to checkbox7 (see GCBO)
 % eventdata  reserved - to be defined in a future version of MATLAB
 % handles    structure with handles and user data (see GUIDATA)

 % Hint: get(hObject,'Value') returns toggle state of checkbox7



 function edit17_Callback(hObject, eventdata, handles)
 % hObject    handle to edit17 (see GCBO)
 % eventdata  reserved - to be defined in a future version of MATLAB
 % handles    structure with handles and user data (see GUIDATA)

 % Hints: get(hObject,'String') returns contents of edit17 as text
 %        str2double(get(hObject,'String')) returns contents of edit17 as a double


 % --- Executes during object creation, after setting all properties.
 function edit17_CreateFcn(hObject, eventdata, handles)
 % hObject    handle to edit17 (see GCBO)
 % eventdata  reserved - to be defined in a future version of MATLAB
 % handles    empty - handles not created until after all CreateFcns called

 % Hint: edit controls usually have a white background on Windows.
 %       See ISPC and COMPUTER.
 if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
     set(hObject,'BackgroundColor','white');
 end


 % --- Executes on button press in enableRefine.
 function enableRefine_Callback(hObject, eventdata, handles)
 % hObject    handle to enableRefine (see GCBO)
 % eventdata  reserved - to be defined in a future version of MATLAB
 % handles    structure with handles and user data (see GUIDATA)

 % Hint: get(hObject,'Value') returns toggle state of enableRefine



 function mutposFile_Callback(hObject, eventdata, handles)
 % hObject    handle to mutposFile (see GCBO)
 % eventdata  reserved - to be defined in a future version of MATLAB
 % handles    structure with handles and user data (see GUIDATA)

 [filename pathname]= uigetfile('*.txt','Pick a mutpos file');
 if(filename)
     h = handles.mutpos;
     a = textread(strcat(pathname,filename), '%d');
     set(h, 'String', a);
 end


 % --- Executes during object creation, after setting all properties.
 function mutposFile_CreateFcn(hObject, eventdata, handles)
 % hObject    handle to mutposFile (see GCBO)
 % eventdata  reserved - to be defined in a future version of MATLAB
 % handles    empty - handles not created until after all CreateFcns called

 % Hint: edit controls usually have a white background on Windows.
 %       See ISPC and COMPUTER.
 if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
     set(hObject,'BackgroundColor','white');
 end


 % --- Executes on button press in pushbutton11.
 function pushbutton11_Callback(hObject, eventdata, handles)
 % hObject    handle to pushbutton11 (see GCBO)
 % eventdata  reserved - to be defined in a future version of MATLAB
 % handles    structure with handles and user data (see GUIDATA)
 [filename pathname]= uigetfile('*.txt','Pick a file for mutation position');
 if(filename)
     set(handles.mutposFile, 'String', strcat(pathname, filename));
 end

 % --- Executes on button press in pushbutton12.
 function pushbutton12_Callback(hObject, eventdata, handles)
 % hObject    handle to pushbutton12 (see GCBO)
 % eventdata  reserved - to be defined in a future version of MATLAB
 % handles    structure with handles and user data (see GUIDATA)
 [filename pathname]= uigetfile('*.txt','Pick a file for marks');
 if(filename)
     set(handles.marksFile, 'String', strcat(pathname, filename));
 end


 function marksFile_Callback(hObject, eventdata, handles)
 % hObject    handle to marksFile (see GCBO)
 % eventdata  reserved - to be defined in a future version of MATLAB
 % handles    structure with handles and user data (see GUIDATA)

 % Hints: get(hObject,'String') returns contents of marksFile as text
 %        str2double(get(hObject,'String')) returns contents of marksFile as a double


 % --- Executes during object creation, after setting all properties.
 function marksFile_CreateFcn(hObject, eventdata, handles)
 % hObject    handle to marksFile (see GCBO)
 % eventdata  reserved - to be defined in a future version of MATLAB
 % handles    empty - handles not created until after all CreateFcns called

 % Hint: edit controls usually have a white background on Windows.
 %       See ISPC and COMPUTER.
 if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
     set(hObject,'BackgroundColor','white');
 end


 % --- Executes on button press in saveBtn.
 function saveBtn_Callback(hObject, eventdata, handles)
 % hObject    handle to saveBtn (see GCBO)
 % eventdata  reserved - to be defined in a future version of MATLAB
 % handles    structure with handles and user data (see GUIDATA)
 global SETTING;
 sequence = SETTING.sequence;
 primer_distance_from_end = SETTING.primer_distance_from_end;
 area_peak = handles.area_peak;
 offset = SETTING.offset;
 [name path] = uiputfile('Output.rdat', 'Save to RDAT file');
 if(name)
   fullname = strcat(path, name);
   seqpos = length(sequence)-primer_distance_from_end - [0:(size(area_peak,1)-1)] + offset;
   rdat = fill_rdat(name, sequence, offset, seqpos,area_peak);
   output_rdat_to_file(fullname,rdat);
 end


 % --- Executes on button press in saveBtn.
 function saveWorkspace_Callback(hObject, eventdata, handles)
 % hObject    handle to saveBtn (see GCBO)
 % eventdata  reserved - to be defined in a future version of MATLAB
 % handles    structure with handles and user data (see GUIDATA)
 global SETTING;
 global stepNo;

 handles_to_SETTINGS( handles );
 
 uisave({'handles', 'SETTING', 'stepNo'}, 'workspace.mat');

% --- Executes on button press in loadWorkspace
function loadWorkspace_Callback(hObject, eventdata, handles)
% hObject    handle to saveBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global SETTING;
global stepNo;

 [name path] = uigetfile('*.mat', 'Load MATLAB workspace file');
if(name)
  fullname = strcat(path, name);
  
  vars = load( fullname );
  SETTING = vars.SETTING;
  
  % need to fill out 'set' commands
  set( handles.refCol, 'String', num2str( SETTING.refcol ) );
  set( handles.ymin, 'String', num2str( SETTING.ymin ) );
  set( handles.ymax, 'String', num2str( SETTING.ymax ) );
  set( handles.offset, 'String', num2str( SETTING.offset ) );
  set( handles.primer_distance_from_end, 'String', num2str( SETTING.primer_distance_from_end ) );
  set( handles.sequence, 'String',  SETTING.sequence );
  set( handles.abiListbox, 'String',  SETTING.list, 'Max', length(SETTING.list), 'Value', 1 );
  set( handles.enableDP, 'Value',  SETTING.enableDP );
  set( handles.slack, 'String',  num2str(SETTING.slack) );
  set( handles.maxShift, 'String',  num2str(SETTING.maxShift) );
  set( handles.windowSize, 'String',  num2str(SETTING.windowSize) );
  set( handles.targetBlock, 'String',  SETTING.targetBlock );
  set( handles.enableBaseline, 'Value',  SETTING.enableBaseline );
  set( handles.enableRefine, 'Value',  SETTING.enableRefine );
  set( handles.mutposFile, 'String',  SETTING.mutposFile );

end


 % --- Needed for saving, and starting run.
function handles_to_SETTINGS( handles )

global SETTING;

SETTING.refcol = str2double(get(handles.refCol, 'String'));
SETTING.ymin = str2double(get(handles.ymin, 'String'));
SETTING.ymax = str2double(get(handles.ymax, 'String'));
SETTING.offset = str2double(get(handles.offset, 'String'));
SETTING.primer_distance_from_end = str2double(get(handles.primer_distance_from_end, 'String'));
SETTING.sequence = get(handles.sequence, 'String');    
SETTING.list = get(handles.abiListbox, 'String');
SETTING.enableDP = get(handles.enableDP, 'Value');
SETTING.slack = str2double(get(handles.slack, 'String'));
SETTING.maxShift = str2double(get(handles.maxShift, 'String'));
SETTING.windowSize = str2double(get(handles.windowSize, 'String'));
SETTING.targetBlock = get(handles.targetBlock, 'String');
SETTING.enableBaseline = get(handles.enableBaseline, 'Value');
SETTING.enableRefine = get(handles.enableRefine, 'Value');
SETTING.mutposFile = get(handles.mutposFile, 'String');



function primer_distance_from_end_Callback(hObject, eventdata, handles)
% hObject    handle to primer_distance_from_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of primer_distance_from_end as text
%        str2double(get(hObject,'String')) returns contents of primer_distance_from_end as a double


% --- Executes during object creation, after setting all properties.
function primer_distance_from_end_CreateFcn(hObject, eventdata, handles)
% hObject    handle to primer_distance_from_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
