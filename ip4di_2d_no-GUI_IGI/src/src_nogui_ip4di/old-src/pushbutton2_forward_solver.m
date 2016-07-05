% --- Executes on button press in pushbutton2.
%function pushbutton2_Callback(hObject, eventdata, handles)
%function pushbutton2_forward_solver(input,mesh,fem)
function pushbutton2_forward_solver
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global input mesh fem

t1=cputime;

fem=forward_solver(input,mesh)

t2=cputime;
disp(sprintf('TIME FOR FORWARD SIMULATION = %f',t2-t1));

%set(handles.slider1,'Enable','On');
%set(handles.slider1,'Min',1);
%set(handles.slider1,'Max',input.num_mes);
%set(handles.slider1,'SliderStep',[1/(input.num_mes-1) 1/(input.num_mes-1)]);
%set(handles.pushbutton3,'Enable','On');
%set(handles.Resolution_button,'Enable','On');
%set(handles.pushbutton6,'Enable','On');
