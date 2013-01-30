function varargout = forward_modelling(varargin)
% FORWARD_MODELLING M-file for forward_modelling.fig
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

% Last Modified by GUIDE v2.5 01-Aug-2011 16:55:05

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
global input mesh
global first_time_call

set(handles.pushbutton2,'Enable','Off');
set(handles.slider1,'Enable','Off');
set(handles.pushbutton3,'Enable','Off');
set(handles.Resolution_button,'Enable','Off');

first_time_call=1;
cla

[filename,pathname]=uigetfile('*.d','Select Data File');

input.mes_in=fullfile(pathname, filename);
input.sip_flag=0;
input.dc_flag=1;
input.time_lapse_flag=0;
input.ip_flag=0;
input.res2d_flag=0;
[input]=read_data(input);
input.data_type=1;
%read_data_small;
[mesh]=create_mesh3(input,0);

set(handles.uitable1,'Data',mesh.depth_n);

if mesh.num_of_bor~=0
    set(handles.pushbutton5,'Enable','off');
else
    set(handles.pushbutton5,'Enable','on');
end
set(handles.edit2,'String',num2str(mesh.max_n));


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global input mesh first_time_call pp_plot ip_plot

   hold on;


faces2=[];
tmp1=1;





% first time call
if first_time_call==1

%     load('skata2.mat');
%     mesh.res_param1=skata2(:,3);
%         mesh.res_param1=tmp_model(:,3);
        
        
         cla
         mesh.node(:,2)=-mesh.node(:,2);
         input.bgr_res_param=mesh.res_param1;

         for i=1:mesh.num_param
             tmp=mesh.faces{i};
             tmp2=mesh.edge(tmp,:);
             tmp3=mesh.node(tmp2,:);
             tmp4=unique(tmp3,'rows');
             
             a=find(tmp4(:,2)==max(tmp4(:,2))); % up
             
             b=find(tmp4(:,2)==min(tmp4(:,2))); % up
             b=b(end:-1:1);
             
             c=find(tmp4(:,1)==min(tmp4(:,1))); % left
             cc=setdiff(c,a);
             c=setdiff(cc,b);
             
             d=find(tmp4(:,1)==max(tmp4(:,1))); % right
             dd=setdiff(d,a);
             d=setdiff(dd,b);
             
             faces=[a;d;b;c];

     patch('faces',faces','vertices',tmp4,'facecolor','flat','FaceVertexCData',(abs(input.bgr_res_param(i)))','edgecolor','k');
    
    hold on
         end
plot(mesh.orig_probe_x,-mesh.orig_probe_z,'o');
            colorbar
            axis equal;% off;
            first_time_call=0;
            

            
            
            
end

xlim([min(mesh.tmp_param(:,3)) max(mesh.tmp_param(:,4))]);
ylim([-max(mesh.tmp_param(:,6)) -min(mesh.tmp_param(:,5))]);

[x y]=ginput(1);
% x=x-left(num_ext);
%bgr_res_param=[1:1:num_param];

if (x<min(mesh.tmp_param(:,3)) || x>max(mesh.tmp_param(:,4))) ||( y<-max(mesh.tmp_param(:,6)) ||  y>-min(mesh.tmp_param(:,5)))
    

    errordlg('Area selected if out of bounds','Selection Error');
else

    % Search which parameter user chose
    tmp_x=abs(mesh.param_x-x);

    [dummy,ind3]=min(tmp_x);
    x_need=mesh.param_x(ind3);
    [ind3]=find(mesh.param_x==x_need);


    tmp_y=abs( -y -mesh.param_y);

    [dummy,ind4]=min(tmp_y);
    y_need=mesh.param_y(ind4);
    [ind4]=find(mesh.param_y==y_need);

    % so we want the parameter withcoordinats (ind3,ind4);
    selected_param=intersect(ind3,ind4);
    % get new value for current parameter    
if input.sip_flag==0
    current=inputdlg('Enter new resistivity value','Forward Model',1,{ num2str(input.bgr_res_param(selected_param))} );
    user_entry = str2double(current); 
elseif input.sip_flag==1
   com_val{1}=num2str(abs(input.bgr_res_param(selected_param)));
   com_val{2}=num2str (  1000*atan2 (   imag(input.bgr_res_param(selected_param)),real(input.bgr_res_param(selected_param))    )   )  ;     
   current=inputdlg({'Enter new amplitude value','Enter new phase value'},'Forward Model',1, com_val );
   mag=str2double(current{1});
   phi=str2double(current{2});
   user_entry=complex (mag*(cos(phi/1000)) , mag*(sin(phi/1000)));
end            
    
     %make some tests in case negative or nan
    if isnan(user_entry)
        errordlg('NaN','Not a Number','modal');
        uicontrol(hObject);
    elseif user_entry<0
        errordlg('Negative value','Negative Value','modal');
    elseif isempty(user_entry)==1
        % do nothing
    else
       input.bgr_res_param(selected_param)=user_entry; 
        
    end
    
    

             tmp=mesh.faces{selected_param};
             tmp2=mesh.edge(tmp,:);
             tmp3=mesh.node(tmp2,:);
             tmp4=unique(tmp3,'rows');
             
             a=find(tmp4(:,2)==max(tmp4(:,2))); % up
             
             b=find(tmp4(:,2)==min(tmp4(:,2))); % up
             b=b(end:-1:1);
             
             c=find(tmp4(:,1)==min(tmp4(:,1))); % left
             cc=setdiff(c,a);
             c=setdiff(cc,b);
             
             d=find(tmp4(:,1)==max(tmp4(:,1))); % right
             dd=setdiff(d,a);
             d=setdiff(dd,b);
             
             faces=[a;d;b;c];

     patch('faces',faces','vertices',tmp4,'facecolor','flat','FaceVertexCData',(abs( input.bgr_res_param(selected_param)))','edgecolor','k');
    plot(mesh.orig_probe_x,-mesh.orig_probe_z,'ko');
    set(handles.pushbutton2,'Enable','On');    
    

end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global input mesh fem


fem=forward_solver(input,mesh)

set(handles.slider1,'Enable','On');
set(handles.slider1,'Min',1);
set(handles.slider1,'Max',input.num_mes);
set(handles.slider1,'SliderStep',[1/(input.num_mes-1) 1/(input.num_mes-1)]);
set(handles.pushbutton3,'Enable','On');
set(handles.Resolution_button,'Enable','On');
set(handles.pushbutton6,'Enable','On');




% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global mesh fem input pp_plot

slider_value = get(hObject,'Value');
set(handles.edit1,'String',num2str(int32(slider_value)));

l=int32(slider_value);


cla




    F = TriScatteredInterp(mesh.param_x,-mesh.param_y,fem.array_jacobian(l,:)');
    qz=F(mesh.xxx,mesh.yyy);
    contourf(mesh.xxx,mesh.yyy,qz);
%     shading interp




colorbar
axis equal;% off; 



hold on
plot(input.ax(l),-input.az(l),'ko','LineWidth',4);

text(input.ax(l),-input.az(l)+mesh.probe_spacing,'A','FontSize',22);

plot(input.bx(l),-input.bz(l),'ko','LineWidth',4);
text(input.bx(l),-input.bz(l)+mesh.probe_spacing,'B','FontSize',22);

plot(input.mx(l),-input.mz(l),'ko','LineWidth',4);
text(input.mx(l),-input.mz(l)+mesh.probe_spacing,'M','FontSize',22);

plot(input.nx(l),-input.nz(l),'ko','LineWidth',4);
text(input.nx(l),-input.nz(l)+mesh.probe_spacing,'N','FontSize',22);




% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


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


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global mesh input fem pp_plot
axes (pp_plot)
for l=1:input.num_mes
    cla
    set(handles.edit1,'String',num2str(l));


        F = TriScatteredInterp(mesh.param_x',-mesh.param_y',fem.array_jacobian(l,:)');
        qz=F(mesh.xxx,mesh.yyy);
        %contourf(xxx,yyy,qz);
        contourf(mesh.xxx,mesh.yyy,qz);
        shading interp


    colorbar
    axis equal;% off; 




hold on
plot(input.ax(l),-input.az(l),'ko','LineWidth',4);

text(input.ax(l),-input.az(l)+mesh.probe_spacing,'A','FontSize',22);

plot(input.bx(l),-input.bz(l),'ko','LineWidth',4);
text(input.bx(l),-input.bz(l)+mesh.probe_spacing,'B','FontSize',22);

plot(input.mx(l),-input.mz(l),'ko','LineWidth',4);
text(input.mx(l),-input.mz(l)+mesh.probe_spacing,'M','FontSize',22);

plot(input.nx(l),-input.nz(l),'ko','LineWidth',4);
text(input.nx(l),-input.nz(l)+mesh.probe_spacing,'N','FontSize',22);
    pause(0.01)
end


% --- Executes on button press in radiobutton3.
function radiobutton3_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton3


% --- Executes on button press in blocks.
function blocks_Callback(hObject, eventdata, handles)
% hObject    handle to blocks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of blocks
if get(gcbo,'Value')==1
    set(handles.contour,'Value',0);
else
    set(handles.contour,'Value',1);
end

% --- Executes on button press in contour.
function contour_Callback(hObject, eventdata, handles)
% hObject    handle to contour (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of contour
if get(gcbo,'Value')==1
    set(handles.blocks,'Value',0);
else
    set(handles.blocks,'Value',1);
end


% --- Executes on button press in Resolution_button.
function Resolution_button_Callback(hObject, eventdata, handles)
% hObject    handle to Resolution_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Dear user since i have allready calculated the jacobian (sensitivity)
% it's now easy for you to calculate the resolution matrix. Since this
% program is meant just for forward modelling, i didn't bother calculatind
% a decent smoothenss matrix. I just use DAMPED LEAST SQUARES. You can
% modify the equation below, so you can insert your smoothness matrix.
% In case someone needs a smoothnrss matrix, please look at the acompaning
% m files below and uncomment below lines, and and change the
% eye(num_param) with ctc (assuming c is the name of the smoothness).

% global num_of_bor
%if num_of_bor==1 || num_of_bor==2
%    smooth_mtx_horizontal;
%    
%elseif num_of_bor==0
%    smooth_mtx_surface;   
%end
global mesh input fem pp_plot

% R = (JTJ +l*I)-1 * JT*J
axes(pp_plot)
resolution_matrix=(fem.array_jacobian.'*fem.array_jacobian + 0.005*eye(mesh.num_param))\(fem.array_jacobian.'*fem.array_jacobian);
% i chose to plot on the diagonal of the matrix
 for i=1:mesh.num_param
    resolution_matrix2(i,1)=resolution_matrix(i,i);
 end
 
 
 % And now plot it...
cla 

        F = TriScatteredInterp(mesh.param_x,-mesh.param_y,resolution_matrix2);
        qz=F(mesh.xxx,mesh.yyy);
        %contourf(xxx,yyy,qz);
        contourf(mesh.xxx,mesh.yyy,qz);
%         shading interp


    colorbar
    axis equal;% off; 
 




% --------------------------------------------------------------------
function Untitled_3_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global input mesh
global first_time_call

set(handles.pushbutton2,'Enable','Off');
set(handles.slider1,'Enable','Off');
set(handles.pushbutton3,'Enable','Off');
set(handles.Resolution_button,'Enable','Off');

first_time_call=1;
cla

[filename,pathname]=uigetfile('*.d','Select Data File');

input.mes_in=fullfile(pathname, filename);
input.sip_flag=1;
input.dc_flag=0;
input.time_lapse_flag=0;
input.ip_flag=0;
input.res2d_flag=0;
[input]=read_data(input);
input.data_type=1;
%read_data_small;
[mesh]=create_mesh3(input,0);


mesh.depth_n=reshape(mesh.depth_n,mesh.max_n,1);
set(handles.uitable1,'Data',mesh.depth_n);

if mesh.num_of_bor~=0
    set(handles.pushbutton5,'Enable','off');
else
    set(handles.pushbutton5,'Enable','on');
end
set(handles.edit2,'String',num2str(mesh.max_n));

% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global mesh input
depth_n=sort(get(handles.uitable1,'data'));
cla(gca)
mesh=create_mesh3(input,depth_n);


% --- Executes during object creation, after setting all properties.
function uitable1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uitable1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double
global mesh 

depth_n=zeros(3,1);
keep=mesh.max_n;
val=int32(str2num(get(gcbo,'String')));
if isempty(val)
    msgbox('Not a valid number');
    set(gcbo,'String',keep)
end




depth_n(1)=0;
depth_n(2)=mesh.probe_spacing/2;
m_factor=1.1; % This mean 10% thicker each layer
for i=3:val
        tmp_thick=depth_n(i-1)-depth_n(i-2) ;
%         depth_n(i)=depth_n(i-1)+scale*probe_spacing/(2*probe_spacing);
          depth_n(i)=depth_n(i-1)+(m_factor)*tmp_thick;
end



set(handles.uitable1,'Data',depth_n);

% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global mesh input

dd_pseudo(mesh,input);

% --- Executes during object creation, after setting all properties.
function pushbutton6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1
global pp_plot

pp_plot=gca;


% --------------------------------------------------------------------
function Untitled_4_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fig=data_2d;
uiwait(fig);
msgbox('Data file created. Now read it...');


% --------------------------------------------------------------------
function Untitled_6_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes during object creation, after setting all properties.
function axes3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes3
global ip_plot

ip_plot=gca;


% --------------------------------------------------------------------
function Untitled_5_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global input first_time_call

first_time_call=1;
[input]=read_data_gui(input);
