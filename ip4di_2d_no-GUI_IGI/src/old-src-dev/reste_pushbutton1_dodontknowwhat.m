
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
%shading interp

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
