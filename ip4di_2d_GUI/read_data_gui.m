function varargout = read_data_gui(varargin)
% READ_DATA_GUI MATLAB code for read_data_gui.fig
%      READ_DATA_GUI, by itself, creates a new READ_DATA_GUI or raises the existing
%      singleton*.
%
%      H = READ_DATA_GUI returns the handle to a new READ_DATA_GUI or the handle to
%      the existing singleton*.
%
%      READ_DATA_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in READ_DATA_GUI.M with the given input arguments.
%
%      READ_DATA_GUI('Property','Value',...) creates a new READ_DATA_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before read_data_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to read_data_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help read_data_gui

% Last Modified by GUIDE v2.5 20-Dec-2010 13:20:10

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @read_data_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @read_data_gui_OutputFcn, ...
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


% --- Executes just before read_data_gui is made visible.
function read_data_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to read_data_gui (see VARARGIN)

% Choose default command line output for read_data_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes read_data_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = read_data_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in dc_checkbox.
function dc_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to dc_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of dc_checkbox
global input

input.res2d_flag=0;
input.sip_flag=0;
input.ip_flag=0;
input.dc_flag=1;

set(handles.ip_checkbox,'Value',0);
set(handles.sip_checkbox,'Value',0);
test_file(hObject, eventdata, handles); 

% --- Executes on button press in ip_checkbox.
function ip_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to ip_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ip_checkbox
global input

input.res2d_flag=0;
input.sip_flag=0;
input.ip_flag=1;
input.dc_flag=0;


set(handles.dc_checkbox,'Value',0);
set(handles.sip_checkbox,'Value',0);
test_file(hObject, eventdata, handles); 

% --- Executes on button press in sip_checkbox.
function sip_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to sip_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of sip_checkbox
global input

input.res2d_flag=0;
input.sip_flag=1;
input.ip_flag=0;
input.dc_flag=0;

set(handles.ip_checkbox,'Value',0);
set(handles.dc_checkbox,'Value',0);
test_file(hObject, eventdata, handles); 

% --- Executes on button press in wd_checkbox.
function wd_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to wd_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of wd_checkbox
global input

input.wd_flag=get(gcbo,'Value');

% --- Executes on button press in single_button.
function single_button_Callback(hObject, eventdata, handles)
% hObject    handle to single_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of single_button
global input

val=get(gcbo,'Value');
if val==1; input.time_lapse_flag=0;end
set (handles.radiobutton2,'Value',0);
test_file(hObject, eventdata, handles); 
    

% --- Executes on button press in radiobutton2.
function radiobutton2_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton2
global input

val=get(gcbo,'Value');
if val==1; input.time_lapse_flag=1;end
set (handles.single_button,'Value',0);

test_file(hObject, eventdata, handles);    

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global input mesh

val=get(handles.checkbox5,'Value');
if val==0
    input.data_type=1;
elseif val==1
    input.data_type=2;
end
try
    if input.time_lapse_flag==0
        [input]=read_data(input);
    else
        [input]=read_4d_data(input);   
    end
     close(gcf)
     [mesh]=create_mesh3(input,0);
catch
    msgbox('Something Went Wrong');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global input big_part ip_part
% cla (big_part);
% cla (ip_part);

if input.time_lapse_flag==0
    [filename,pathname]=uigetfile('*.d','Select Data File');
else
    [filename,pathname]=uigetfile('*.dat','Select Data File');    
end

input.mes_in=fullfile(pathname, filename);
set(handles.pushbutton2,'Enable','On');

if input.time_lapse_flag==0
    tmp{1}=(filename);
    set(handles.uitable1,'Data',tmp);
else
    tmp=importdata(input.mes_in);
    set(handles.uitable1,'Data',tmp);    
end

% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
global input

input.res2d_flag=0;
input.sip_flag=0;
input.ip_flag=0;
input.dc_flag=0;

input.time_lapse_flag=0;


% --- Executes during object creation, after setting all properties.
function dc_checkbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dc_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function ip_checkbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ip_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function sip_checkbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sip_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function wd_checkbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to wd_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function single_button_CreateFcn(hObject, eventdata, handles)
% hObject    handle to single_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function radiobutton2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to radiobutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function pushbutton1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function uitable1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uitable1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function pushbutton2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



function test_file(hObject, eventdata, handles)

% first check if DC/IP/SIP is selected
val1=get(handles.dc_checkbox,'Value');
val2=get(handles.ip_checkbox,'Value');
val3=get(handles.sip_checkbox,'Value');


if val1==1 || val2==1 || val3==1
    data_flag=1;
else 
    data_flag=0;
end


% Now check time lapse options
val4=get(handles.single_button,'Value');
val5=get(handles.radiobutton2,'Value');

if val4==1 || val5==1
    time_flag=1;
else
    time_flag=0;
end


if data_flag==1 && time_flag==1
    set(handles.pushbutton1,'Enable','On');
else
    set(handles.pushbutton1,'Enable','Off');
end
