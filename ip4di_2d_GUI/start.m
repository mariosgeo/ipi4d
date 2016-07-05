function varargout = start(varargin)
% START M-file for start.fig
%      START, by itself, creates a new START or raises the existing
%      singleton*.
%
%      H = START returns the handle to a new START or the handle to
%      the existing singleton*.
%
%      START('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in START.M with the given input arguments.
%
%      START('Property','Value',...) creates a new START or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before start_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to start_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help start

% Last Modified by GUIDE v2.5 01-Jan-2012 22:19:46

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @start_OpeningFcn, ...
                   'gui_OutputFcn',  @start_OutputFcn, ...
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


% --- Executes just before start is made visible.
function start_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to start (see VARARGIN)

% Choose default command line output for start
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes start wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = start_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global proceed_flag mesh input tim

if input.time_lapse_flag==1
    h=figure;
    if input.dc_flag==1
        for i=1:input.num_files
            tim(i)=subplot(1,input.num_files,i);
        end
    elseif input.ip_flag==1 || input.sip_flag==1
        for i=1:2*input.num_files
            tim(i)=subplot(2,input.num_files,i);
        end
    end
end
%proceed_flag=1;
if proceed_flag==1
    main(input,mesh)
else
    msgbox('No Inversion Parameters Are Selected. Please use from inversion menu.');
end


% --- Executes during object creation, after setting all properties.
function axes4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
global acb_fig1;
% Hint: place code in OpeningFcn to populate axes4
acb_fig1=gca;
title(acb_fig1,'Lagrangian Distribution');
xlabel(acb_fig1,'Distance (m)');
ylabel(acb_fig1,'Depth (m)');





% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(handles.figure1);

% --- Executes during object creation, after setting all properties.
function pushbutton2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



% --- Executes during object creation, after setting all properties.
function text2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
global textt;
textt=hObject;






% --------------------------------------------------------------------
function Untitled_5_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global mesh
inversion_parameters(mesh);






% --- Executes during object creation, after setting all properties.
function axes8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes8
global big_part

big_part=gca;
cla (big_part);
%title(bigpart,'Inversion Model');

% --- Executes during object creation, after setting all properties.
function axes9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes9
global ip_part
ip_part=gca;
cla (ip_part);

% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
















% --------------------------------------------------------------------
function Untitled_4_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%this dc
global big_part ip_part input mesh
cla (big_part);
cla (ip_part);
[filename,pathname]=uigetfile('*.d','Select Data File');
input.mes_in=fullfile(pathname, filename);



input.res2d_flag=0;
input.sip_flag=0;
input.ip_flag=0;
input.dc_flag=1;
input.time_lapse_flag=0;


[input]=read_data(input);
[mesh]=create_mesh3(input,0);

% --------------------------------------------------------------------
function Untitled_11_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%sip
global big_part ip_part input mesh
cla (big_part);
cla (ip_part);
[filename,pathname]=uigetfile('*.d','Select Data File');
input.mes_in=fullfile(pathname, filename);



input.res2d_flag=0;
input.sip_flag=1;
input.ip_flag=0;
input.dc_flag=0;
input.time_lapse_flag=0;


[input]=read_data(input);
[mesh]=create_mesh3(input,0);


% --------------------------------------------------------------------
function Untitled_12_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%ip
global big_part ip_part input mesh
cla (big_part);
cla (ip_part);
[filename,pathname]=uigetfile('*.d','Select Data File');
input.mes_in=fullfile(pathname, filename);



input.res2d_flag=0;
input.sip_flag=0;
input.ip_flag=1;
input.dc_flag=0;
input.time_lapse_flag=0;


[input]=read_data(input);
[mesh]=create_mesh3(input,0);

% --------------------------------------------------------------------
function Untitled_13_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%time-laspe dc
global big_part ip_part input mesh
cla (big_part);
cla (ip_part);
[filename,pathname]=uigetfile('*.dat','Select Data File');
input.mes_in=fullfile(pathname, filename);



input.res2d_flag=0;
input.sip_flag=0;
input.ip_flag=0;
input.dc_flag=1;
input.time_lapse_flag=1;


[input]=read_4d_data(input);
[mesh]=create_mesh3(input,0);



% --------------------------------------------------------------------
function Untitled_15_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%time-lapse ip
global big_part ip_part input mesh
cla (big_part);
cla (ip_part);
[filename,pathname]=uigetfile('*.dat','Select Data File');
input.mes_in=fullfile(pathname, filename);



input.res2d_flag=0;
input.sip_flag=0;
input.ip_flag=1;
input.dc_flag=0;
input.time_lapse_flag=1;


[input]=read_4d_data(input);
[mesh]=create_mesh3(input,0);





% --------------------------------------------------------------------
function Untitled_16_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%time-lapse sip
global big_part ip_part input mesh
cla (big_part);
cla (ip_part);
[filename,pathname]=uigetfile('*.dat','Select Data File');
input.mes_in=fullfile(pathname, filename);



input.res2d_flag=0;
input.sip_flag=1;
input.ip_flag=0;
input.dc_flag=0;
input.time_lapse_flag=1;


[input]=read_4d_data(input);
[mesh]=create_mesh3(input,0);





% --------------------------------------------------------------------
function Untitled_9_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% res2dinv
global big_part ip_part input mesh
cla (big_part);
cla (ip_part);
[filename,pathname]=uigetfile('*.d','Select Data File');
input.mes_in=fullfile(pathname, filename);



input.res2d_flag=1;
input.sip_flag=0;
input.ip_flag=0;
input.dc_flag=0;
input.time_lapse_flag=0;


[input]=read_data(input);
[mesh]=create_mesh3(input,0);


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_3_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_6_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_2_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_7_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
inversion_results


% --------------------------------------------------------------------
function read_file_dialog_Callback(hObject, eventdata, handles)
% hObject    handle to read_file_dialog (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
read_data_gui


% --------------------------------------------------------------------
function uipushtool3_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
