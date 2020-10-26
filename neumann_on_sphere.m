function varargout = neumann_on_sphere(varargin)
        % NEUMANN_ON_SPHERE MATLAB code for neumann_on_sphere.fig
        %      NEUMANN_ON_SPHERE, by itself, creates a new NEUMANN_ON_SPHERE or raises the existing
        %      singleton*.
        %
        %      H = NEUMANN_ON_SPHERE returns the handle to a new NEUMANN_ON_SPHERE or the handle to
        %      the existing singleton*.
        %
        %      NEUMANN_ON_SPHERE('CALLBACK',hObject,eventData,handles,...) calls the local
        %      function named CALLBACK in NEUMANN_ON_SPHERE.M with the given input arguments.
        %
        %      NEUMANN_ON_SPHERE('Property','Value',...) creates a new NEUMANN_ON_SPHERE or raises the
        %      existing singleton*.  Starting from the left, property value pairs are
        %      applied to the GUI before neumann_on_sphere_OpeningFcn gets called.  An
        %      unrecognized property name or invalid value makes property application
        %      stop.  All inputs are passed to neumann_on_sphere_OpeningFcn via varargin.
        %
        %      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
        %      instance to run (singleton)".
        %
        % See also: GUIDE, GUIDATA, GUIHANDLES
        
        % Edit the above text to modify the response to help neumann_on_sphere
        
        % Last Modified by GUIDE v2.5 30-Mar-2019 21:39:27
        
        % Begin initialization code - DO NOT EDIT
        gui_Singleton = 1;
        gui_State = struct('gui_Name',       mfilename, ...
                'gui_Singleton',  gui_Singleton, ...
                'gui_OpeningFcn', @neumann_on_sphere_OpeningFcn, ...
                'gui_OutputFcn',  @neumann_on_sphere_OutputFcn, ...
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
        
        %to-do list:
        %merge progress bar into GUI
        %create spherical harmonics table seperatly from the code
        %implement Newton-Raphson to accurate the saddle point location
        %fix bugs
        
        
        % --- Executes just before neumann_on_sphere is made visible.
function neumann_on_sphere_OpeningFcn(hObject, eventgudata, handles, varargin)
        
        %defining defaults
        resolution = 500;
        angular_resolution = 360;
        color_resolution = 100;
        points_size = 40;
        l = 5;
%         l=1;
        m_dist_style = 'uni';
%         coefficients = rand(1,2*l+1);
%         coefficients = [0.4344 0.1537 0.2229 0.3485 0.1706 0.1 0.1441 0.1672 0.1263 0.1685 0.3224];
%         coefficients = [0.4344 0.1537 0.2229 0.3485 0.1706 0.548 0.1441 0.1672 0.1263 0.1685 0.3224];
%         coefficients = [0,1,0];
        coefficients = [0,0,0,1,0,0,0,0,0,0,0];

        
        %creating default structures
        %resolution structure with fields:
        handles.res.resolution = resolution;
        handles.res.angular_resolution = angular_resolution;
        handles.res.color_resolution = color_resolution;
        handles.res.points_size = points_size;
        
        %coefficients structure with fields:
        handles.coeffs.l = l;
        handles.coeffs.style = m_dist_style;
        handles.coeffs.vals = coefficients/norm(coefficients); %vector of numvers. Must be of length 2l+1
        handles.coeffs.hess_chk = false;
        
        %some flags 
        handles.saveable = false;
        handles.show_dom = false;
        handles.show_dom2 = true;
        
        %GUI stuff
        set(handles.dom_table,'data',cell(1,3));
        set(handles.data_table,'data',cell(7,1));
        set(handles.l_text, 'string', sprintf('L-Value: %d', handles.coeffs.l));
        
        
        guidata(hObject, handles);
        uiwait(handles.figure1);
        
        
        % --- Executes on button press in run.
function run_Callback(hObject, eventdata, handles)
        cla;
        tic; %start timing the running
        wb = waitbar(0,'Running...'); %this is the handle for the waitbar which changes along the code, approximatrly with the time that has to pass 
        set(hObject,'string','Running...' );
        set(hObject,'enable','off');
        set(handles.m_dist,'enable','off');
        set(handles.reso,'enable','off');
        

        %gets the coordinate system
        handles.coo = Coordinates(handles.res.resolution, handles.res.angular_resolution, handles.res.color_resolution, handles.res.points_size);
        
        %gets all data about function and neuman lines
        handles.SH = get_SH(handles.coeffs, handles.coo, wb);
        
        %plots everything
        handles.graphs = plot_function(handles.coo,handles.SH,handles.axes,handles.rho_dist,wb);
        
        %wrap up everything
        set(handles.data_table,'data',handles.graphs.data');        
        handles.saveable = true;
        t = toc;
        close(wb);
        set(hObject,'string','Done!');
        set(hObject,'enable','on');
        set(handles.time_text, 'string', sprintf('Elapsed Time: %f s', t));
        set(handles.m_dist,'enable','on');
        set(handles.reso,'enable','on');        
        guidata(hObject,handles);
        
        
        % --- Executes on button press in m_dist.
function m_dist_Callback(hObject, eventdata, handles)
        handles.coeffs = m_coeffs(handles.coeffs);
        set(handles.l_text, 'string', sprintf('L-Value: %d', handles.coeffs.l));
        set(handles.run, 'string', 'Run');
        guidata(hObject,handles);
        
        
        % --- Executes on button press in nodal.
function nodal_Callback(hObject, eventdata, handles)
        if get(hObject,'value')
                handles.graphs.nodal.Visible = 'on';                
        else
                handles.graphs.nodal.Visible = 'off';
        end
        
        
        % --- Executes on button press in extrima.
function extrima_Callback(hObject, eventdata, handles)
        if get(hObject,'value')
                handles.graphs.maxima.Visible = 'on';
                handles.graphs.minima.Visible = 'on';
        else
                handles.graphs.maxima.Visible = 'off';
                handles.graphs.minima.Visible = 'off';
        end
        
        
        % --- Executes on button press in saddle.
function saddle_Callback(hObject, eventdata, handles)
        if get(hObject,'value')
                handles.graphs.saddles.Visible = 'on';
        else
                handles.graphs.saddles.Visible = 'off';
        end
        
        
        % --- Executes on button press in neus.
function neus_Callback(hObject, eventdata, handles)
        if get(hObject,'value')
                for i=1:length(handles.graphs.neumann_lines)
                        handles.graphs.neumann_lines(i).Visible = 'on';
                end
        else
                for i=1:length(handles.graphs.neumann_lines)
                        handles.graphs.neumann_lines(i).Visible = 'off';
                end
        end
        
        
        % --- Executes when entered data in editable cell(s) in data_table.
function data_table_CellEditCallback(hObject, eventdata, handles)
        
        
function reso_Callback(hObject, eventdata, handles)
        handles.res = resolution_function(handles.res);
        guidata(hObject,handles);
        
        
        % --------------------------------------------------------------------
function save_ClickedCallback(hObject, eventdata, handles)
        folder = fullfile(pwd,'saved');
        if ~exist(folder,'dir')
                mkdir(folder);
        end
        if handles.saveable
                SH = handles.SH;
                coo = handles.coo;
                coeffs = handles.coeffs;
                cd(folder);
                file_name = strcat(sprintf('%d',coeffs.l),'_',coeffs.style,datestr(now,'_HH_MM__dd_mm_yy'),'.mat');
                [file,path] = uiputfile('*.mat','Save Simulation',file_name);
                if file~=0
                        save(fullfile(path,file),'SH','coo','coeffs');
                end
                cd ..;
        end
                
        
        % --------------------------------------------------------------------
function load_ClickedCallback(hObject, eventdata, handles)
        folder = fullfile(pwd,'saved');
        if ~exist(folder,'dir')
                mkdir(folder);
        end
        cd(folder);
        [file,path] = uigetfile('*.mat','Load Data');
        cd ..
        if file~=0
                cla;
                wb = 0;
                load(fullfile(path,file),'SH','coo','coeffs');
                
                handles.SH = SH;
                handles.coo = coo;
                handles.coeffs = coeffs;
                
                handles.graphs = plot_function(handles.coo,handles.SH,handles.axes,wb);
                
                set(handles.data_table,'data',handles.graphs.data');
                
                handles.saveable = true;
                guidata(hObject,handles);
        end
        
        
        % --- Outputs from this function are returned to the command line.
function varargout = neumann_on_sphere_OutputFcn(hObject, eventdata, handles)
        
        
        % --- Executes on button press in rho_color.
function rho_color_Callback(hObject, eventdata, handles)
        if get(hObject,'value')             
                for i=1:length(handles.graphs.neu_domains_marks)
                        handles.graphs.neu_domains_marks(i).Visible = 'on';
                        handles.show_dom2 = true;
                        guidata(hObject,handles);
                        set(handles.figure1,'WindowButtonDownFcn',{@axes_ButtonDownFcn, handles});
                end
        else
                for i=1:length(handles.graphs.neu_domains_marks)
                        handles.graphs.neu_domains_marks(i).Visible = 'off';
                        handles.show_dom2 = false;
                        handles.show_dom = false;
                        if ~isempty(handles.graphs.neu_marked)
                                delete(handles.graphs.neu_marked);
                        end
                        guidata(hObject,handles);    
                        set(handles.figure1,'WindowButtonDownFcn',{@axes_ButtonDownFcn, handles});
                end
        end        


% --- Executes on mouse press over axes background.
function axes_ButtonDownFcn(hObject, eventdata, handles)
        if handles.show_dom&&handles.saveable
                mouse = get(hObject,'CurrentPoint');
                if (mouse(1)<=140)&&(mouse(2)<=40)
                        mouse = get(handles.axes,'CurrentPoint');
                        [~, i] = min(pdist2(mouse(1,:), handles.graphs.neu_domains));
                        dom = handles.SH.neumann(i);
                        p = handles.graphs.neu_domains(i,:);
                        p = 1.01*p/norm(p);
                        if ~isempty(handles.graphs.neu_marked)
                                delete(handles.graphs.neu_marked);
                        end
                        handles.graphs.neu_marked = scatter3(p(1),p(2),p(3),70,[0.75, 0, 0.75],'o');
                        handles.graphs.neu_marked.Annotation.LegendInformation.IconDisplayStyle = 'off';
                        set(handles.dom_table, 'data', {dom.area, dom.circumference, dom.rho});
                        guidata(hObject,handles);
                        set(handles.figure1,'WindowButtonDownFcn',{@axes_ButtonDownFcn, handles});
                end
        end
        
        
        % --------------------------------------------------------------------
function rotate_ClickedCallback(hObject, eventdata, handles)
    set(handles.zoom,'state','off');
    set(handles.pan,'state','off');
    set(handles.data_cursor,'state','off');
    set(handles.dom_data,'state','off');
    rotate3d on;
    pan off;
    zoom off;    
    datacursormode off;
    handles.show_dom = false;
    guidata(hObject,handles);
    set(handles.figure1,'WindowButtonDownFcn',{@axes_ButtonDownFcn, handles});


% --------------------------------------------------------------------
function pan_ClickedCallback(hObject, eventdata, handles)    
    set(handles.zoom,'state','off');
    set(handles.rotate,'state','off');
    set(handles.data_cursor,'state','off');
    set(handles.dom_data,'state','off');
    rotate3d off;
    pan on;
    zoom off;  
    datacursormode off;
    handles.show_dom = false;
    guidata(hObject,handles);
    set(handles.figure1,'WindowButtonDownFcn',{@axes_ButtonDownFcn, handles});


% --------------------------------------------------------------------
function zoom_ClickedCallback(hObject, eventdata, handles)
    set(handles.rotate,'state','off');
    set(handles.pan,'state','off');
    set(handles.data_cursor,'state','off');
    set(handles.dom_data,'state','off');
    rotate3d off;
    pan off;
    datacursormode off;
    zoom on;
    handles.show_dom = false;
    guidata(hObject,handles);
    set(handles.figure1,'WindowButtonDownFcn',{@axes_ButtonDownFcn, handles});


% --------------------------------------------------------------------
function data_cursor_ClickedCallback(hObject, eventdata, handles)
        set(handles.rotate,'state','off');
        set(handles.pan,'state','off');
        set(handles.zoom,'state','off');
        set(handles.dom_data,'state','off');
        rotate3d off;
        pan off;
        zoom off;
        datacursormode on;
        handles.show_dom = false;
        guidata(hObject,handles);
        set(handles.figure1,'WindowButtonDownFcn',{@axes_ButtonDownFcn, handles});
        

% --------------------------------------------------------------------
function dom_data_ClickedCallback(hObject, eventdata, handles)
        set(handles.data_cursor,'state','off');
        set(handles.rotate,'state','off');
        set(handles.pan,'state','off');
        set(handles.zoom,'state','off');
        rotate3d off;
        pan off;
        zoom off;
        datacursormode off;
        if handles.show_dom2
                handles.show_dom = true;
                guidata(hObject,handles);
                set(handles.figure1,'WindowButtonDownFcn',{@axes_ButtonDownFcn, handles});
        end
       
        
% --- Executes on button press in stats_dist1.
function stats_dist1_Callback(hObject, eventdata, handles)
        handles.graphs.rho.Visible = 'on';
        handles.graphs.extDeg(1).Visible = 'off';
        handles.graphs.extDeg(2).Visible = 'off';
        axis(handles.rho_dist, handles.graphs.rho_M);
        set(handles.uipanel7, 'title', 'Rho Distribution');


% --- Executes on button press in stats_dist2.
function stats_dist2_Callback(hObject, eventdata, handles)
        handles.graphs.rho.Visible = 'off';
        handles.graphs.extDeg(1).Visible = 'on';
        handles.graphs.extDeg(2).Visible = 'on';
        axis(handles.rho_dist, handles.graphs.extDeg_M);
        set(handles.uipanel7, 'title', 'Extrema Degrees Distribution');
