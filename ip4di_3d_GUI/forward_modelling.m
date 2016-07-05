function varargout = forward_modelling(varargin)
% FORWARD_MODELLING MATLAB code for forward_modelling.fig
%      FORWARD_MODELLING, by itself, creates a new FORWARD_MODELLING or raises the existing
%      singleton*.
%
%      H = FORWARD_MODELLING returns the handle to a new FORWARD_MODELLING or the handle to
%      the existing singleton*.
%
%      FORWARD_MODELLING('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FORWARD_MODELLING.M with the given input arguments.
%
%      FORWARD_MODELLING('Property','Value',...) creates a new FORWARD_MODELLING or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before forward_modelling_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to forward_modelling_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help forward_modelling

% Last Modified by GUIDE v2.5 01-Jan-2012 23:59:28

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @forward_modelling_OpeningFcn, ...
                   'gui_OutputFcn',  @forward_modelling_OutputFcn, ...
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


% --- Executes just before forward_modelling is made visible.
function forward_modelling_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to forward_modelling (see VARARGIN)

% Choose default command line output for forward_modelling
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes forward_modelling wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = forward_modelling_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_2_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global input mesh xi yi zi

clear data
fig=read_data_gui;
uiwait(fig);


set(handles.edit1,'String',mesh.max_n);
set(handles.uitable2,'Data',mesh.depth_n);

data(:,1)=mesh.param_x;
data(:,2)=mesh.param_y;
data(:,3)=mesh.param_z;

data(:,4)=real(mesh.mean_res);
data(:,4)=real(mesh.mean_res);
% for i=1:length(data(:,1))
%    if data(i,3)>=4.8
%        data(i,4)=50;
%    end
% end
for i=1:length(data(:,1))
    if data(i,1)>=30.5
    data(i,4)=500;
    end
end
if input.sip_flag==1 ;data(:,5)=imag(mesh.mean_res);end

set(handles.uitable1,'Data',data);

% create here the grid once
x_unique=unique(mesh.param_x);
y_unique=unique(mesh.param_y);
z_unique=unique(mesh.param_z);

min_x=min(x_unique);
max_x=max(x_unique);

min_y=min(y_unique);
max_y=max(y_unique);

min_z=min(z_unique);
max_z=max(z_unique);

[xi yi zi]=meshgrid(x_unique,y_unique,z_unique);


% --------------------------------------------------------------------
function Untitled_3_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fig=data_2d;
uiwait(fig);
msgbox('Data file created. Now read it...');


% --- Executes when entered data in editable cell(s) in uitable1.
function uitable1_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to uitable1 (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
global mesh xi yi zi real_part imag_part input
data=get(gcbo,'Data');
z_unique=unique(mesh.param_z);
tmp2=data(:,4);

 FF=TriScatteredInterp(mesh.param_x,mesh.param_y,mesh.param_z,tmp2);
vi=FF(xi,yi,zi);

slice(real_part,xi,yi,-zi,vi,[],[],[-z_unique]);
 title(real_part,'Amplitude (Ohm.m)');
 xlabel(real_part,'Distance (m)');
 ylabel(real_part,'Distance (m)');
 zlabel(real_part,'Depth (m)');
 colorbar ('peer',real_part);
 view(real_part,30,10);
 
 if input.sip_flag==1
     tmp2=data(:,5);

 FF=TriScatteredInterp(mesh.param_x,mesh.param_y,mesh.param_z,tmp2);
vi=FF(xi,yi,zi);

slice(imag_part,xi,yi,-zi,vi,[],[],[-z_unique]);
 title(imag_part,'Phase (mrad)');
 xlabel(imag_part,'Distance (m)');
 ylabel(imag_part,'Distance (m)');
 zlabel(imag_part,'Depth (m)');
 colorbar('peer',imag_part);
 view(imag_part,30,10);
 end


% --- Executes during object creation, after setting all properties.
function axes3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes3
global real_part

real_part=gca;


% --- Executes during object creation, after setting all properties.
function axes2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes2
global imag_part

imag_part=gca;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global input mesh
data=get(handles.uitable1,'Data');
if input.sip_flag==0
    mesh.bgr_param=data(:,4);
else
    mesh.bgr_param=complex (data(:,4).*(cos(data(:,5)./1000)) , data(:,4).*(sin(data(:,5)./1000)));    
end
forward_solver(input,mesh);



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
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
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global mesh input big_part xi yi zi
max_n=str2num(get(handles.edit1,'String'));
depth_n=sort(get(handles.uitable2,'data'));


mesh=mesh_gen_3d_new_version(input,mesh,depth_n);
set(handles.edit1,'String',mesh.max_n);
set(handles.uitable2,'Data',mesh.depth_n);

data(:,1)=mesh.param_x;
data(:,2)=mesh.param_y;
data(:,3)=mesh.param_z;

data(:,4)=real(mesh.mean_res);

    

if input.sip_flag==1 ;data(:,5)=imag(mesh.mean_res);end

set(handles.uitable1,'Data',data);

% create here the grid once
x_unique=unique(mesh.param_x);
y_unique=unique(mesh.param_y);
z_unique=unique(mesh.param_z);

min_x=min(x_unique);
max_x=max(x_unique);

min_y=min(y_unique);
max_y=max(y_unique);

min_z=min(z_unique);
max_z=max(z_unique);

[xi yi zi]=meshgrid(x_unique,y_unique,z_unique);
