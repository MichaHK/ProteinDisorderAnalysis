function varargout = GUI_FEBS(varargin)
% GUI_FEBS MATLAB code for GUI_FEBS.fig
%      GUI_FEBS, by itself, creates a new GUI_FEBS or raises the existing
%      singleton*.
%
%      H = GUI_FEBS returns the handle to a new GUI_FEBS or the handle to
%      the existing singleton*.
%
%      GUI_FEBS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_FEBS.M with the given input arguments.
%
%      GUI_FEBS('Property','Value',...) creates a new GUI_FEBS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_FEBS_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_FEBS_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_FEBS

% Last Modified by GUIDE v2.5 29-Apr-2015 15:01:51

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_FEBS_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_FEBS_OutputFcn, ...
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


% --- Executes just before GUI_FEBS is made visible.
function GUI_FEBS_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_FEBS (see VARARGIN)

% Choose default command line output for GUI_FEBS
handles.output = hObject;
handles.CustomCommandsHistoryFilename = ['CustomCommandHistory.' getenv('COMPUTERNAME') '.txt'];
set(handles.CheckBoxType1,'Value',1)
set(handles.TailBox,'Value',1);
set(handles.PhosCheck,'Value',1);

% Don't let this interrupt the program from loading
try
    handles.CustomCommandHistory = {};
    if (exist(handles.CustomCommandsHistoryFilename, 'file'))
        data = ReadTextFileLines(handles.CustomCommandsHistoryFilename);
        handles.CustomCommandHistory = data;
    end
catch
end

% Update handles structure
guidata(hObject, handles);
set(handles.PopUpY,'Value',2);
main_FEBS(handles)

set(handles.filterPopup,'Value',3);
set(handles.filterRangeEditbox, 'String', '[10, 1e4]');

% UIWAIT makes GUI_FEBS wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_FEBS_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in PopUpX.
function PopUpX_Callback(hObject, eventdata, handles)
% hObject    handle to PopUpX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
main_FEBS(handles)
% Hints: contents = cellstr(get(hObject,'String')) returns PopUpX contents as cell array
%        contents{get(hObject,'Value')} returns selected item from PopUpX


% --- Executes during object creation, after setting all properties.
function PopUpX_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PopUpX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in PopUpY.
function PopUpY_Callback(hObject, eventdata, handles)
% hObject    handle to PopUpY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
main_FEBS(handles)
% Hints: contents = cellstr(get(hObject,'String')) returns PopUpY contents as cell array
%        contents{get(hObject,'Value')} returns selected item from PopUpY


% --- Executes during object creation, after setting all properties.
function PopUpY_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PopUpY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in YesPhosButton.
function YesPhosButton_Callback(hObject, eventdata, handles)
% hObject    handle to YesPhosButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
main_FEBS(handles)
% Hint: get(hObject,'Value') returns toggle state of YesPhosButton


% --- Executes on button press in NoPhosButton.
function NoPhosButton_Callback(hObject, eventdata, handles)
% hObject    handle to NoPhosButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
main_FEBS(handles)
% Hint: get(hObject,'Value') returns toggle state of NoPhosButton


% --- Executes on button press in CheckBoxType1.
function CheckBoxType1_Callback(hObject, eventdata, handles)
% hObject    handle to CheckBoxType1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
main_FEBS(handles)
% Hint: get(hObject,'Value') returns toggle state of CheckBoxType1


% --- Executes on button press in CheckBoxType2.
function CheckBoxType2_Callback(hObject, eventdata, handles)
% hObject    handle to CheckBoxType2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
main_FEBS(handles)
% Hint: get(hObject,'Value') returns toggle state of CheckBoxType2


% --- Executes on button press in CheckBoxType3.
function CheckBoxType3_Callback(hObject, eventdata, handles)
% hObject    handle to CheckBoxType3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
main_FEBS(handles)
% Hint: get(hObject,'Value') returns toggle state of CheckBoxType3


% --- Executes on button press in CheckBoxType5.
function CheckBoxType5_Callback(hObject, eventdata, handles)
% hObject    handle to CheckBoxType5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
main_FEBS(handles)
% Hint: get(hObject,'Value') returns toggle state of CheckBoxType5


% --- Executes on button press in CheckBoxType4.
function CheckBoxType4_Callback(hObject, eventdata, handles)
% hObject    handle to CheckBoxType4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
main_FEBS(handles)
% Hint: get(hObject,'Value') returns toggle state of CheckBoxType4


% --- Executes on button press in CheckBoxType6.
function CheckBoxType6_Callback(hObject, eventdata, handles)
% hObject    handle to CheckBoxType6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
main_FEBS(handles)
% Hint: get(hObject,'Value') returns toggle state of CheckBoxType6


% --- Executes on button press in radiobutton4.
function radiobutton4_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
main_FEBS(handles)
% Hint: get(hObject,'Value') returns toggle state of radiobutton4


% --- Executes on button press in PushSwitchXY.
function PushSwitchXY_Callback(hObject, eventdata, handles)
% hObject    handle to PushSwitchXY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Xind=get(handles.PopUpX,'Value');
Yind=get(handles.PopUpY,'Value');
set(handles.PopUpY,'Value',Xind);
set(handles.PopUpX,'Value',Yind);
main_FEBS(handles)

% --- Executes on button press in PlotPush.
function PlotPush_Callback(hObject, eventdata, handles)
% hObject    handle to PlotPush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
main_FEBS(handles)


% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1


% --- Executes on button press in PhosCheck.
function PhosCheck_Callback(hObject, eventdata, handles)
% hObject    handle to PhosCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
main_FEBS(handles)
% Hint: get(hObject,'Value') returns toggle state of PhosCheck


% --- Executes on button press in UnPhosCheck.
function UnPhosCheck_Callback(hObject, eventdata, handles)
% hObject    handle to UnPhosCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
main_FEBS(handles)
% Hint: get(hObject,'Value') returns toggle state of UnPhosCheck


% --- Executes on button press in LogY_Check.
function LogY_Check_Callback(hObject, eventdata, handles)
% hObject    handle to LogY_Check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
main_FEBS(handles)
% Hint: get(hObject,'Value') returns toggle state of LogY_Check


% --- Executes on button press in LogX_Check.
function LogX_Check_Callback(hObject, eventdata, handles)
% hObject    handle to LogX_Check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
main_FEBS(handles)
% Hint: get(hObject,'Value') returns toggle state of LogX_Check


% --- Executes on button press in ABS_X_check.
function ABS_X_check_Callback(hObject, eventdata, handles)
% hObject    handle to ABS_X_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
main_FEBS(handles)
% Hint: get(hObject,'Value') returns toggle state of ABS_X_check


% --- Executes on button press in ABS_Y_check.
function ABS_Y_check_Callback(hObject, eventdata, handles)
% hObject    handle to ABS_Y_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
main_FEBS(handles)
% Hint: get(hObject,'Value') returns toggle state of ABS_Y_check


% --- Executes on button press in PushClipboardFigure.
function PushClipboardFigure_Callback(hObject, eventdata, handles)
% hObject    handle to PushClipboardFigure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
1;
hgexport(handles.figure1,'-clipboard')

%set(handles.axes1, 'Units', 'pixels');
%axesCoordinates = get(handles.axes1, 'Position');

% saveas(handles.figure1, 'temp.jpg');
% print(handles.figure1, '-dbmp16m', 'temp.bmp');
% I = imread('temp.bmp');
% imclipboard('copy', I)


% --- Executes during object creation, after setting all properties.
function TextPanel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TextPanel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in HeadBox.
function HeadBox_Callback(hObject, eventdata, handles)
% hObject    handle to HeadBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
main_FEBS(handles)

% Hint: get(hObject,'Value') returns toggle state of HeadBox


% --- Executes on button press in TailBox.
function TailBox_Callback(hObject, eventdata, handles)
% hObject    handle to TailBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
main_FEBS(handles)

% Hint: get(hObject,'Value') returns toggle state of TailBox


% --- Executes on button press in RodBox.
function RodBox_Callback(hObject, eventdata, handles)
% hObject    handle to RodBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
main_FEBS(handles)

% Hint: get(hObject,'Value') returns toggle state of RodBox


% --- Executes on button press in CheckPlotExternal.
function CheckPlotExternal_Callback(hObject, eventdata, handles)
% hObject    handle to CheckPlotExternal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
main_FEBS(handles)

% Hint: get(hObject,'Value') returns toggle state of CheckPlotExternal


function [scatterHandles] = GetScatterPlots(hObject, handles)
childrenHandles = allchild(handles.axes1);
stringEndsWith = @(s, sub)(numel(s) >= numel(sub) && strcmp(s(end-numel(sub)+1:end), sub));
whichScatter = [arrayfun(@(a)stringEndsWith(class(a), 'Scatter'), childrenHandles)] ~= 0;
scatterHandles = childrenHandles(whichScatter);

% --- Executes on button press in colorAccordingToKeratinTypeButton.
function colorAccordingToKeratinTypeButton_Callback(hObject, eventdata, handles)
% hObject    handle to colorAccordingToKeratinTypeButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axesChildrenHandles = allchild(handles.axes1);
whichScatter = [arrayfun(@(a)~isempty(get(a, 'UserData')), axesChildrenHandles)] ~= 0;

%scatterHandles = GetScatterPlots(hObject, handles);
scatterHandles = axesChildrenHandles(whichScatter);

newColors = hot(6);
for i = 1:numel(scatterHandles)
    sh = scatterHandles(i);
    if (~isstruct(sh)) % For compatibility before/after 2014b
        s = get(sh);
    else
        s = sh;
    end
    
    if (isempty(s.UserData) || ~isstruct(s.UserData))
        continue;
    end
    
    proteinsData = s.UserData.ProteinData;
    whichKeratin = ([proteinsData.IsKeratin] ~= 0);
    keratinTypes = [proteinsData(whichKeratin).KeratinType];
    keratinTissues = {proteinsData(whichKeratin).KeratinTissue};
    
    % 1,2 -> Epithelial type 1/2
    % 3,4 -> Hair type 1/2
    whichHair = [cellfun(@(c)c(1)=='H', keratinTissues)];
    
    cdata = zeros(numel(proteinsData), 3);
    cdata(whichKeratin, :) = [newColors(whichHair * 2 + keratinTypes, :)];
    set(sh, 'CData', cdata);
    1;
end


% --- Executes on button press in redrawButton.
function redrawButton_Callback(hObject, eventdata, handles)
% hObject    handle to redrawButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
main_FEBS(handles)


% --- Executes on button press in curveFitButton.
function curveFitButton_Callback(hObject, eventdata, handles)
% hObject    handle to curveFitButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
scatterHandles = GetScatterPlots(hObject, handles)

x = [];
y = [];

for i = 1:numel(scatterHandles)
    s = scatterHandles(i);
    x = [x; s.XData(:)];
    y = [y; s.YData(:)];
end

cftool(x, y);

function [lines] = ReadTextFileLines(filename)
lines = {};

f = fopen(filename, 'r');
if (f >= 0)
    while (~feof(f))
        l = strtrim(fgetl(f));
        if (~isempty(l))
            lines{end+1} = l;
        end
    end
    
    fclose(f);
end

function [] = WriteTextFileLines(filename, lines)
f = fopen(filename, 'w');
if (f >= 0)
    for i = 1:numel(lines)
        l = strtrim(lines{i});
        if (~isempty(l))
            fprintf(f, '%s\n', l);
        end
    end
    fclose(f);
end

function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
str = get(hObject,'String');
%set(hObject,'String', '');

try
    whichIdentical = cellfun(@(c)strcmp(c, str), handles.CustomCommandHistory);
    
    if (nnz(whichIdentical) == 0)
        handles.CustomCommandHistory = horzcat(str, handles.CustomCommandHistory);
        
        if (numel(handles.CustomCommandHistory) > 10)
            handles.CustomCommandHistory = handles.CustomCommandHistory(1:10);
        end
        
        guidata(hObject, handles);
        
        WriteTextFileLines(handles.CustomCommandsHistoryFilename, handles.CustomCommandHistory);
    end
catch
end

if (iscell(str))
    cellfun(@(c)eval(c), str);
else
    eval(str);
end

% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on key press with focus on edit1 and none of its controls.
function edit1_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

str = get(handles.edit1, 'String');
whichIdentical = cellfun(@(c)strcmp(c, str), handles.CustomCommandHistory);
whichIdentical = find(whichIdentical);

if (isempty(whichIdentical))
    whichIdentical = 0;
end

index = whichIdentical;

if (numel(handles.CustomCommandHistory) > 0)
switch (eventdata.Key)
    case 'uparrow'
        index = index - 1;
        if (index < 1)
            index = numel(handles.CustomCommandHistory);
        end
        set(handles.edit1, 'String', handles.CustomCommandHistory(index));
        
    case 'downarrow'
        index = index + 1;
        if (index > numel(handles.CustomCommandHistory))
            index = 1;
        end
        set(handles.edit1, 'String', handles.CustomCommandHistory(index));
        
    otherwise
        1;
end
end


% --- Executes on button press in displayLabelsCheckbox.
function displayLabelsCheckbox_Callback(hObject, eventdata, handles)
% hObject    handle to displayLabelsCheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of displayLabelsCheckbox
main_FEBS(handles)


% --- Executes on selection change in filterPopup.
function filterPopup_Callback(hObject, eventdata, handles)
% hObject    handle to filterPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns filterPopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from filterPopup
main_FEBS(handles)


% --- Executes during object creation, after setting all properties.
function filterPopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filterPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function filterRangeEditbox_Callback(hObject, eventdata, handles)
% hObject    handle to filterRangeEditbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of filterRangeEditbox as text
%        str2double(get(hObject,'String')) returns contents of filterRangeEditbox as a double
main_FEBS(handles)


% --- Executes during object creation, after setting all properties.
function filterRangeEditbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filterRangeEditbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in filterCheckbox.
function filterCheckbox_Callback(hObject, eventdata, handles)
% hObject    handle to filterCheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of filterCheckbox
main_FEBS(handles)
