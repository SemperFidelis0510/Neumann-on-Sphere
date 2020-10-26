function varargout = resolution_function(varargin)
% RESOLUTION_FUNCTION MATLAB code for resolution_function.fig
%      RESOLUTION_FUNCTION, by itself, creates a new RESOLUTION_FUNCTION or raises the existing
%      singleton*.
%
%      H = RESOLUTION_FUNCTION returns the handle to a new RESOLUTION_FUNCTION or the handle to
%      the existing singleton*.
%
%      RESOLUTION_FUNCTION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RESOLUTION_FUNCTION.M with the given input arguments.
%
%      RESOLUTION_FUNCTION('Property','Value',...) creates a new RESOLUTION_FUNCTION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before resolution_function_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to resolution_function_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help resolution_function

% Last Modified by GUIDE v2.5 13-Mar-2019 20:21:27

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @resolution_function_OpeningFcn, ...
                   'gui_OutputFcn',  @resolution_function_OutputFcn, ...
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


% --- Executes just before resolution_function is made visible.
function resolution_function_OpeningFcn(hObject, eventdata, handles, varargin)
    handles.res = varargin{1};
    set(handles.reg_res,'string',handles.res.resolution);
    set(handles.ang_res,'string',handles.res.angular_resolution);
    guidata(hObject, handles);
    uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = resolution_function_OutputFcn(hObject, eventdata, handles) 
    varargout{1} = handles.res;
    close(handles.figure1);


% --- Executes during object creation, after setting all properties.
function reg_res_CreateFcn(hObject, eventdata, handles)
% hObject    handle to reg_res (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



function ang_res_Callback(hObject, eventdata, handles)
% hObject    handle to ang_res (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ang_res as text
%        str2double(get(hObject,'String')) returns contents of ang_res as a double


% --- Executes during object creation, after setting all properties.
function ang_res_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ang_res (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in enter.
function enter_Callback(hObject, eventdata, handles)
    handles.res.resolution = str2num(get(handles.reg_res,'string'));
    handles.res.angular_resolution = str2num(get(handles.ang_res,'string'));
    guidata(hObject,handles);
    uiresume(handles.figure1);


% --- Executes on button press in cancel.
function cancel_Callback(hObject, eventdata, handles)
    guidata(hObject,handles);
    uiresume(handles.figure1);
