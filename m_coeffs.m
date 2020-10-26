function varargout = m_coeffs(varargin)
% M_COEFFS MATLAB code for m_coeffs.fig
%      M_COEFFS, by itself, creates a new M_COEFFS or raises the existing
%      singleton*.
%
%      H = M_COEFFS returns the handle to a new M_COEFFS or the handle to
%      the existing singleton*.
%
%      M_COEFFS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in M_COEFFS.M with the given input arguments.
%
%      M_COEFFS('Property','Value',...) creates a new M_COEFFS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before m_coeffs_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to m_coeffs_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help m_coeffs

% Last Modified by GUIDE v2.5 15-Mar-2019 16:32:52

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @m_coeffs_OpeningFcn, ...
                   'gui_OutputFcn',  @m_coeffs_OutputFcn, ...
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


% --- Executes just before m_coeffs is made visible.
function m_coeffs_OpeningFcn(hObject, eventdata, handles, varargin)
    handles.coeffs = varargin{1};    
    handles.coeffs_ = handles.coeffs;
    handles.flag = 0;
    
    
    set(handles.numbers,'ColumnName',(-handles.coeffs_.l:1:handles.coeffs_.l));
    set(handles.numbers,'data',handles.coeffs.vals);
    set(handles.l_value,'string',handles.coeffs_.l);
    set(handles.gauss_mu,'string',(-handles.coeffs_.l:1:handles.coeffs_.l));
    set(handles.gauss_mu,'value',handles.coeffs_.l+1);
    switch handles.coeffs.style
        case 'uni'
            set(handles.uniform,'value',1);
            set(handles.gauss,'value',0); 
            set(handles.manual,'value',0); 
        case 'gauss'
            set(handles.uniform,'value',0);
            set(handles.gauss,'value',1); 
            set(handles.manual,'value',0);
            set(handles.gauss_parameters,'visible','on');
        case 'manual'
            set(handles.uniform,'value',0);
            set(handles.gauss,'value',0); 
            set(handles.manual,'value',1);
            set(handles.numbers,'ColumnEditable',true(1,2*handles.coeffs_.l+1));
    end


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes m_coeffs wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = m_coeffs_OutputFcn(hObject, eventdata, handles) 
    if handles.flag
        varargout{1} = handles.coeffs_;
    else
        varargout{1} = handles.coeffs;
    end
    close(handles.figure1);


% --- Executes on button press in manual.
function manual_Callback(hObject, eventdata, handles)
    if get(hObject,'Value')
        set(handles.uniform,'value',0);
        set(handles.gauss,'value',0);   
        set(handles.gauss_parameters,'visible','off');  
        handles.coeffs_.style = 'manual';
        handles.flag = 0;
        guidata(hObject,handles);
    end



% --- Executes on button press in gauss.
function gauss_Callback(hObject, eventdata, handles)
    if get(hObject,'Value')
        set(handles.uniform,'value',0);
        set(handles.manual,'value',0);   
        set(handles.gauss_parameters,'visible','on');     
        handles.coeffs_.style = 'gauss';
        handles.flag = 0;
        guidata(hObject,handles);
    end


% --- Executes on button press in uniform.
function uniform_Callback(hObject, eventdata, handles)
    if get(hObject,'Value')
        set(handles.gauss,'value',0);
        set(handles.manual,'value',0); 
        set(handles.gauss_parameters,'visible','off'); 
        handles.coeffs_.style = 'uni';
        handles.flag = 0;
        guidata(hObject,handles);
    end


% --- Executes on button press in enter.
function enter_Callback(hObject, eventdata, handles)
    if strcmp(handles.coeffs_.style,'manual')
        v = get(handles.numbers,'data');
        handles.coeffs_.vals = v/norm(v);           
    end      
    guidata(hObject, handles);
    uiresume(handles.figure1);


% --- Executes on button press in render.
function render_Callback(hObject, eventdata, handles)
    l = str2num(get(handles.l_value,'string')); %#ok<ST2NM>
    handles.coeffs_.l = l;
    set(handles.numbers,'ColumnName',(-l:1:l));
    switch handles.coeffs_.style
        case 'uni'           
            v = rand(1,2*handles.coeffs_.l+1);
            handles.coeffs_.vals = v/norm(v);
            set(handles.numbers,'data',handles.coeffs_.vals);    
            set(handles.numbers,'ColumnEditable',false(1,2*handles.coeffs_.l+1));            
        case 'gauss'
            mu = get(handles.gauss_mu,'string');
            mu = str2num(mu(get(handles.gauss_mu,'value'))); %#ok<ST2NM>
            sigma = str2double(get(handles.gauss_sigma,'string'));            
            v = exp(-((-l:1:l)-mu).^2/(2*sigma^2));
            handles.coeffs_.vals = v/norm(v);
            set(handles.numbers,'data',handles.coeffs_.vals);    
            set(handles.numbers,'ColumnEditable',false(1,2*handles.coeffs_.l+1));
            handles.coeffs_.style = 'gauss';
        case 'manual'
            set(handles.numbers,'ColumnEditable',true(1,2*handles.coeffs_.l+1));
            set(handles.numbers,'data',zeros(1,2*handles.coeffs_.l+1));
            handles.coeffs_.style = 'manual';
    end
    handles.flag = 1;
    guidata(hObject, handles);


% --- Executes on button press in cancel.
function cancel_Callback(hObject, eventdata, handles)
    handles.flag = 0;
    guidata(hObject,handles);
    uiresume(handles.figure1);



function l_value_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function l_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to l_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in gauss_mu.
function gauss_mu_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function gauss_mu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gauss_mu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function gauss_sigma_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function gauss_sigma_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gauss_sigma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
