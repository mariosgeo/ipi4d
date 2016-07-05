function varargout = inversion_parameters(varargin)
%INVERSION_PARAMETERS M-file for inversion_parameters.fig
%      INVERSION_PARAMETERS, by itself, creates a new INVERSION_PARAMETERS or raises the existing
%      singleton*.
%
%      H = INVERSION_PARAMETERS returns the handle to a new INVERSION_PARAMETERS or the handle to
%      the existing singleton*.
%
%      INVERSION_PARAMETERS('Property','Value',...) creates a new INVERSION_PARAMETERS using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to inversion_parameters_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      INVERSION_PARAMETERS('CALLBACK') and INVERSION_PARAMETERS('CALLBACK',hObject,...) call the
%      local function named CALLBACK in INVERSION_PARAMETERS.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help inversion_parameters

% Last Modified by GUIDE v2.5 10-Dec-2010 18:38:40
global mesh
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @inversion_parameters_OpeningFcn, ...
                   'gui_OutputFcn',  @inversion_parameters_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
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


% --- Executes just before inversion_parameters is made visible.
function inversion_parameters_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for inversion_parameters
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes inversion_parameters wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = inversion_parameters_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



% --- Executes on selection change in popupmenu3.
function popupmenu3_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu3


% --- Executes during object creation, after setting all properties.
function popupmenu3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu4.
function popupmenu4_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu4


% --- Executes during object creation, after setting all properties.
function popupmenu4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in log_checkbox.
function log_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to log_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of log_checkbox
global input

input.log_plot=get(gcbo,'Value');
guidata(hObject,handles);

% --- Executes on button press in acb_plot_checkbox.
function acb_plot_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to acb_plot_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of acb_plot_checkbox
global input

input.acb_plot=get(gcbo,'Value');
guidata(hObject,handles);

% --- Executes on button press in electrode_checkbox.
function electrode_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to electrode_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of electrode_checkbox
global input
input.electrode_plot=get(gcbo,'Value');
guidata(hObject,handles);

% --- Executes on button press in resolution_checkbox.
function resolution_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to resolution_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of resolution_checkbox
global input
input.resolution_plot=get(gcbo,'Value');
guidata(hObject,handles);

% --- Executes on button press in parameters_checkbox.
function parameters_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to parameters_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of parameters_checkbox
global input
input.parameter_plot=get(gcbo,'Value');
guidata(hObject,handles);

% --- Executes on button press in fem_checkbox.
function fem_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to fem_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of fem_checkbox
global input

input.fem_plot=get(gcbo,'Value');
guidata(hObject,handles);

% --- Executes on button press in checkbox16.
function checkbox16_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox16
global input

input.pseudo_plot=get(gcbo,'Value');
guidata(hObject,handles);

% --- Executes on button press in checkbox5.
function checkbox5_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox5

global input


if get(gcbo,'Value')==1
    if input.inv_flag==1 || input.inv_flag==2
        input.bgr_res_flag=1;
        
    end
elseif get(gcbo,'Value')==0
    if input.inv_flag==1 || input.inv_flag==2
        input.bgr_res_flag=0;
    end
end

if input.bgr_res_flag==1 || input.bgr_res_flag==2
    set(handles.par_nm,'Enable','On');
elseif input.bgr_res_flag==0
    set(handles.par_nm,'Enable','Off');
end


guidata(hObject,handles);


function num_itr_Callback(hObject, eventdata, handles)
% hObject    handle to num_itr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of num_itr as text
%        str2double(get(hObject,'String')) returns contents of num_itr as a double
global input

input.itn=str2num(get(gcbo,'String'));
if isempty(input.itn)
    msgbox('Error in Iteration Value. Not a number.Intial Value has been automatic entered.','Check Value','error');
    input.itn=5;
    set(gcbo,'String','5');
else
if input.itn>25
    input.itn=25;
    set(gcbo,'String','25');
end
if input.itn<=0 
    input.itn=1;
    set(gcbo,'String','1');
end
set(handles.slider4,'Value',input.itn);
end

guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function num_itr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to num_itr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
global input
set(gcbo,'String','5');
input.itn=get(gcbo,'String');
input.itn=str2num(input.itn);
%set(handles.slider4,'Value',itn);
guidata(hObject,handles);
 

sss=sprintf('%s \n','Select maximum number of iterations. Inversion process will stop','if convergence limits are applied or ','when maximun number of iterations reached.');
set(gcbo,'ToolTipString',sss);

% --- Executes on slider movement.
function slider4_Callback(hObject, eventdata, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global input
input.itn=get(gcbo,'Value');
set(handles.num_itr,'String',num2str(input.itn));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function slider4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
max=25;
min=1;
slider_step(1)=1/(max-min);
slider_step(2)=1/(max-min);

set(gcbo,'sliderstep',slider_step,'max',max,'min',1,'Value',5);
guidata(hObject,handles);


% --- Executes on selection change in inv_flag_popupmenu.
function inv_flag_popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to inv_flag_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns inv_flag_popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from inv_flag_popupmenu
global input
input.par_nm=[];


val=get(gcbo,'Value');

if input.time_lapse_flag==0
    if val==1
        input.inv_flag=1;
        input.bgr_res_flag=0;
        set(handles.checkbox5,'Enable','on');
        set(handles.checkbox5,'Value',0);
    elseif val==2
        input.inv_flag=2;
        input.bgr_res_flag=0;
        set(handles.checkbox5,'Enable','on');
        set(handles.checkbox5,'Value',0);
    elseif val==3
        input.inv_flag=0;
        input.bgr_res_flag=0;
        set(handles.checkbox5,'Enable','on');
        set(handles.checkbox5,'Value',0);        
    elseif val==4
        input.inv_flag=5;
        input.bgr_res_flag=1;
        set(handles.checkbox5,'Value',1);
        set(handles.checkbox5,'Enable','off');
    elseif val==5
        input.inv_flag=6;
        input.bgr_res_flag=1;
        set(handles.checkbox5,'Value',1);
        set(handles.checkbox5,'Enable','off');
    end


elseif input.time_lapse_flag==1
    val=get(gcbo,'Value');
    if val==1
        input.inv_flag=7;
        input.atc_flag=0;
        input.bgr_res_flag=0;
        set(handles.checkbox5,'Value',0);
        set(handles.checkbox5,'Enable','off');  
        
    elseif val==2
        input.inv_flag=7;
        input.atc_flag=1;
        input.bgr_res_flag=1;
        set(handles.checkbox5,'Value',1);
        set(handles.checkbox5,'Enable','on');
    end
end

if (val==1 || val==2) && get(handles.checkbox5,'Value')==1
    input.bgr_res_flag=1;
end

if (input.bgr_res_flag==1 || input.bgr_res_flag==2) 
    set(handles.par_nm,'Enable','On');
    try
        set(handles.pan_nm_filename,'String',filename);
    catch
        set(handles.pan_nm_filename,'String','NO FILE');
    end

elseif input.bgr_res_flag==0
    set(handles.par_nm,'Enable','Off');
    set(handles.pan_nm_filename,'String','NO FILE');   
    clear input.par_nm
end




if input.ip_flag==1 && (input.inv_flag==5 || input.inv_flag==6)
    msgbox('IP and Time-lapse and/or Difference inversion, is not ready in this version. Switching back to defaults','Not ready option','error'); 
    input.inv_flag=1;
    input.bgr_res_param=0;
    clear input.par_nm
end

guidata(hObject,handles);



% --- Executes during object creation, after setting all properties.
function inv_flag_popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to inv_flag_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
global input
set(gcbo,'Value',1);

clear t
if input.time_lapse_flag==0
    set(gcbo,'Value',3);
    t{1}='Occam';
    t{2}='Gauss-Newton';
    t{3}='Levenberg-Marquardt';
    t{4}='Occam Difference';
    t{5}='GN Difference';
    set(gcbo,'String',t);
    input.inv_flag=0;
    sss=sprintf('%s \n','Use Occam or Gauss-Newton for a smoother model. This option is suggegested in cases','of noisy data. In all other cases use Lavenberg-Marquardt. If unsure leave default option.','Notice that in case of difference inversion, a background model is required.');
    set(gcbo,'ToolTipString',sss);
else
    t{1}='4D';
    t{2}='4D-ATC';
    set(gcbo,'String',t);
    input.inv_flag=7;
    input.atc_flag=0;
    set(gcbo,'TooltipString','Use 4D or 4D-ATC. Notice that in case of 4D-ATC, pre-processing file is required.') 
end
   
 
guidata(hObject,handles);

% --- Executes on selection change in jacobian_flag_popupmenu.
function jacobian_flag_popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to jacobian_flag_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns jacobian_flag_popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from jacobian_flag_popupmenu
global input

val=get(gcbo,'Value');
if val==1
    input.jacobian_flag=1;
elseif val==2
    input.jacobian_flag=2;
elseif val==3
    input.jacobian_flag=3;
end

% --- Executes during object creation, after setting all properties.
function jacobian_flag_popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to jacobian_flag_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
global input
input.jacobian_flag=1;
set(gcbo,'Value',1);
guidata(hObject,handles);
sss=sprintf('%s \n','Full Jacobian calculations provide accurate solution but it is time consuming','Quasi newton calculates Jacobian matrix once. Fast but in-accurate solution.');
set(gcbo,'ToolTipString',sss);



% --- Executes on button press in par_nm.
function par_nm_Callback(hObject, eventdata, handles)
% hObject    handle to par_nm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global input

if input.time_lapse_flag==0;
    [input.filename,input.pathname]=uigetfile('*.mod','Select Data File');
else
    [input.filename,input.pathname]=uigetfile('*.mat','Select Data File');
end

input.par_nm=fullfile(input.pathname, input.filename);
set(handles.pan_nm_filename,'String',input.filename);
guidata(hObject,handles);


function lagrn_Callback(hObject, eventdata, handles)
% hObject    handle to lagrn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lagrn as text
%        str2double(get(hObject,'String')) returns contents of lagrn as a double
global input

input.lagrn=str2num(get(gcbo,'String'));
if isempty(input.lagrn)
    msgbox('Error in Lagrangian Value. Not a number.Intial Value has been automatic entered.','Check Value','error');

        input.lagrn=0.1;
        set(gcbo,'String','0.01');
  
    
end


% --- Executes during object creation, after setting all properties.
function lagrn_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lagrn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
global input

        input.lagrn=0.15;
        set(gcbo,'String','0.15');
   
guidata(hObject,handles);
sss=sprintf('%s \n','Set initial Lagrangian value','Large values lead to a smooth model','Small values makes the inversion unstable','Suggested starting value 0.15');
set(gcbo,'ToolTipString',sss);

function acb_min_Callback(hObject, eventdata, handles)
% hObject    handle to acb_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of acb_min as text
%        str2double(get(hObject,'String')) returns contents of acb_min as a double
global input


input.lagrn_min=str2num(get(gcbo,'String'));
if isempty(input.lagrn_min)
    msgbox('Error in ACB Min. Not a number.Intial Value has been automatic entered.','Check Value','error');
    
        input.lagrn_min=0.01;
        set(gcbo,'String','0.01');
   

end

if input.lagrn_max<input.lagrn_min
    msgbox('Maximum value is less than minimun.Intial Values have been automatic entered.','Check Value','error');
    

        input.lagrn_max=1;
        input.lagrn_min=0.01;
        set(gcbo,'String','0.01');
        set(handles.acb_max,'String','1');
    

end  


% --- Executes during object creation, after setting all properties.
function acb_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to acb_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
global input 

        input.lagrn_min=0.01;
        set(gcbo,'String','0.01');

set(gcbo,'Enable','off');

guidata(hObject,handles);
sss=sprintf('%s \n','Lower ACB limit');
set(gcbo,'ToolTipString',sss);


function acb_max_Callback(hObject, eventdata, handles)
% hObject    handle to acb_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of acb_max as text
%        str2double(get(hObject,'String')) returns contents of acb_max as a double
global input 

input.lagrn_max=str2num(get(gcbo,'String'));
if isempty(input.lagrn_max)
    msgbox('Error in ACB Max. Not a number.Intial Value has been automatic entered.','Check Value','error');
    
        input.lagrn_max=1;
        set(gcbo,'String','1');
    
end

if input.lagrn_max<input.lagrn_min
    msgbox('Maximum value is less than minimun.Intial Values have been automatic entered.','Check Value','error');
    
        input.lagrn_max=1;
        input.lagrn_min=0.01;
        set(gcbo,'String','1');
        set(handles.acb_min,'String','0.01');
    

end    



% --- Executes during object creation, after setting all properties.
function acb_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to acb_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
global input

        input.lagrn_max=1;
        set(gcbo,'String','1');

set(gcbo,'Enable','off');
guidata(hObject,handles);
guidata(hObject,handles);
sss=sprintf('%s \n','Upper ACB limit');
set(gcbo,'ToolTipString',sss);

function lagrn_step_Callback(hObject, eventdata, handles)
% hObject    handle to lagrn_step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lagrn_step as text
%        str2double(get(hObject,'String')) returns contents of lagrn_step as a double
global input

input.lagrn_reduction=str2num(get(gcbo,'String'));
if isempty(input.lagrn_lagrn_reduction)
    msgbox('Error in Lagrangian Step. Not a number.Intial Value has been automatic entered.','Check Value','error');
    input.lagrn_reduction=2;
    set(gcbo,'String','2');
end


% --- Executes during object creation, after setting all properties.
function lagrn_step_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lagrn_step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
global input
set(gcbo,'String',2);
input.lagrn_reduction=2;
guidata(hObject,handles);



% --- Executes on button press in res_limit.
function res_limit_Callback(hObject, eventdata, handles)
% hObject    handle to res_limit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of res_limit
global input

input.limit_res=get(gcbo,'Value');

if input.limit_res==1
   set(handles.res_min,'Enable','On');
   set(handles.res_max,'Enable','On');
   set(handles.res_min,'String','0');
   set(handles.res_max,'String','10000000'); 
elseif input.limit_res==0
   set(handles.res_min,'Enable','Off');
   set(handles.res_max,'Enable','Off');
   set(handles.res_min,'String','NoLimit');
   set(handles.res_max,'String','NoLimit');
end


guidata(hObject,handles);

function res_min_Callback(hObject, eventdata, handles)
% hObject    handle to res_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of res_min as text
%        str2double(get(hObject,'String')) returns contents of res_min as a double
global input

input.min_res=str2double(get(gcbo,'String'));
if isempty(input.min_res)
    msgbox('Error in Minimum Value. Not a number.Intial Value has been automatic entered.','Check Value','error');
    input.min_res=0;
    set(gcbo,'String','0');
end
if input.max_res<input.min_res ||input.min_res<0
   msgbox('Maximum Value smaller the minumum.Intial Values have been automatic restored.','Check Value','error'); 
   input.min_res=0;
   set(gcbo,'String','0');
   input.max_res=100000000;
   set(handles.res_max,'String','100000000');
end


% --- Executes during object creation, after setting all properties.
function res_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to res_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
global input
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(gcbo,'Enable','Off');
set(gcbo,'String','NoLimit');
input.min_res=0;
guidata(hObject,handles);
sss=sprintf('%s \n','Lower resistivity limit');
set(gcbo,'ToolTipString',sss);

function res_max_Callback(hObject, eventdata, handles)
% hObject    handle to res_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of res_max as text
%        str2double(get(hObject,'String')) returns contents of res_max as a double
global input

input.max_res=str2double(get(gcbo,'String'));
if isempty(input.max_res)
    msgbox('Error in Minimum Value. Not a number.Intial Value has been automatic entered.','Check Value','error');
    input.max_res=100000000;
    set(gcbo,'String','100000000');
end

if input.max_res<input.min_res || input.min_res<0
   msgbox('Maximum Value smaller the minumum.Intial Values have been automatic restored.','Check Value','error'); 
   input.min_res=0;
   set(handles.res_min,'String','0');
   input.max_res=100000000;
   set(gcbo,'String','100000000');
end

% --- Executes during object creation, after setting all properties.
function res_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to res_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
global input
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(gcbo,'Enable','Off');
set(gcbo,'String','NoLimit');
input.max_res=10000000;
guidata(hObject,handles);
sss=sprintf('%s \n','Lower resistivity limit');
set(gcbo,'ToolTipString',sss);


% --- Executes during object creation, after setting all properties.
function checkbox5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to checkbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
global input
set(gcbo,'Value',0);
input.bgr_res_flag=0;
set(gcbo,'Enable','on');
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function par_nm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to par_nm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

set(gcbo,'Enable','off');
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function acb_plot_checkbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to acb_plot_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
global input

set(gcbo,'Value',0);
input.acb_plot=0;


% --- Executes during object creation, after setting all properties.
function res_limit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to res_limit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
global input

input.limit_res=0;
set(gcbo,'Value',0);
guidata(hObject,handles);
sss=sprintf('%s \n','Limit resistivity values.','Use this choice in cases when extremely low or high resistivity values are found');
set(gcbo,'ToolTipString',sss);

% --- Executes during object creation, after setting all properties.
function log_checkbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to log_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
global input
set(gcbo,'Value',0);
input.log_plot=0;


% --- Executes during object creation, after setting all properties.
function acb_checkbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to acb_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
global input
set(gcbo,'Value',0);
guidata(hObject,handles);

input.acb_flag=0;
sss=sprintf('%s \n','Active constrained balancing (Yi et al, 2003).','Assign Lagrangian between two limits values based on the resolutions oh the parameters','Use this choice in cases of borehole-data');
set(gcbo,'ToolTipString',sss);

% --- Executes on button press in acb_checkbox.
function acb_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to acb_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of acb_checkbox
global input

input.acb_flag=get(gcbo,'Value');
if input.acb_flag==1
    set(handles.acb_max,'Enable','on');
    set(handles.acb_min,'Enable','on');
    
elseif input.acb_flag==0
    set(handles.acb_max,'Enable','off');
    set(handles.acb_min,'Enable','off');
    input.lagrn_min=0.01;
    input.lagrn_max=1;
    set(handles.acb_max,'String','1');
    set(handles.acb_min,'String','0.01');
    
end
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function electrode_checkbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to electrode_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
global input

input.electrode_plot=0;
set(gcbo,'Value',0);


% --- Executes during object creation, after setting all properties.
function resolution_checkbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to resolution_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
global input
input.resolution_plot=0;
set(gcbo,'Value',0);


% --- Executes during object creation, after setting all properties.
function parameters_checkbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to parameters_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
global input
input.parameter_plot=0;
set(gcbo,'Value',0);


% --- Executes during object creation, after setting all properties.
function fem_checkbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fem_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
global input

input.fem_plot=0;
set(gcbo,'Value',0);


% --- Executes during object creation, after setting all properties.
function checkbox16_CreateFcn(hObject, eventdata, handles)
% hObject    handle to checkbox16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
global input

input.pseudo_plot=0;
set(gcbo,'Value',0);


% --- Executes during object creation, after setting all properties.
function pushbutton2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
global procced_flag

proceed_flag=0;
% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global input proceed_flag

%     input.data_type=1;

%     input.elec_array_type=1;           

    input.conv_rate=2;


    input.current=1;


    proceed_flag=1;

    input.loke_colors=0;


% make checks if background file is chosen
if (input.time_lapse_flag==1 && input.inv_flag==7 && input.atc_flag==1 ) ||...
      (input.time_lapse_flag==0 && input.inv_flag>4)
    
    try
        
        if input.time_lapse_flag==0
            importdata(input.par_nm);
            close(gcf)
        else
            load(input.par_nm);
            close(gcf)
        end
        

    
    catch
        msgbox('No ATC file is selected. Select one to proceed');
       
    end
else
    
    
    close(gcf)
end



function edit30_Callback(hObject, eventdata, handles)
% hObject    handle to edit30 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit30 as text
%        str2double(get(hObject,'String')) returns contents of edit30 as a double
global mesh

depth_n=zeros(3,1);
keep=mesh.max_n;
val=str2num(get(gcbo,'String'));
if isempty(val)
    msgbox('Not a valid number');
    set(gcbo,'String',keep)
end



probe_spacing=min(mesh.probe_spacing_x,mesh.probe_spacing_y);
depth_n(1)=0;
depth_n(2)=probe_spacing/2;
m_factor=1.1; % This mean 10% thicker each layer
for i=3:val+1
        tmp_thick=depth_n(i-1)-depth_n(i-2) ;
%         depth_n(i)=depth_n(i-1)+scale*probe_spacing/(2*probe_spacing);
          depth_n(i)=depth_n(i-1)+(m_factor)*tmp_thick;
end



set(handles.uitable2,'Data',depth_n);

% --- Executes during object creation, after setting all properties.
function edit30_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit30 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
global input mesh

set(gcbo,'String',mesh.n_sep);
if length(unique(mesh.num_z_probes))==1
    set(gcbo,'Enable','On');
end



% --- Executes during object creation, after setting all properties.
function uitable2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uitable2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
global mesh


set(gcbo,'Data',mesh.depth_n);
if length(unique(mesh.num_z_probes))==0
    set(gcbo,'Enable','On');
end


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global mesh input big_part
max_n=str2num(get(handles.edit30,'String'));
depth_n=sort(get(handles.uitable2,'data'));
close(gcf)
cla(big_part)
mesh=mesh_gen_3d_new_version(input,mesh,depth_n);
