function varargout = data_2d(varargin)
% data_2d MATLAB code for data_2d.fig
%      data_2d, by itself, creates a new data_2d or raises the existing
%      singleton*.
%
%      H = data_2d returns the handle to a new data_2d or the handle to
%      the existing singleton*.
%
%      data_2d('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in data_2d.M with the given input arguments.
%
%      data_2d('Property','Value',...) creates a new data_2d or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before data_2d_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to data_2d_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help data_2d

% Last Modified by GUIDE v2.5 07-Jan-2011 14:54:09

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @data_2d_OpeningFcn, ...
                   'gui_OutputFcn',  @data_2d_OutputFcn, ...
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


% --- Executes just before data_2d is made visible.
function data_2d_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to data_2d (see VARARGIN)

% Choose default command line output for data_2d
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes data_2d wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = data_2d_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
tmp=str2num(get(gcbo,'String'));
if isempty(tmp)
    msgbox('Error in Electrode spacing. Not a number.Intial Value has been automatic entered.','Check Value','error');

        
        set(gcbo,'String','1');
  
    
end

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



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double
tmp=str2num(get(gcbo,'String'));
if isempty(tmp)
    msgbox('Error in Electrode spacing. Not a number.Intial Value has been automatic entered.','Check Value','error');

        
        set(gcbo,'String','10');
  
    
end

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



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

set(handles.edit3,'String',num2str(get(gcbo,'Value')));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
max=35;
min=1;
slider_step(1)=1/(max-min);
slider_step(2)=1/(max-min);

set(gcbo,'sliderstep',slider_step,'max',max,'min',1,'Value',1);

% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
set(handles.edit4,'String',num2str(get(gcbo,'Value')));
guidata(hObject,handles);



% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
max=5;
min=1;
slider_step(1)=1/(max-min);
slider_step(2)=1/(max-min);

set(gcbo,'sliderstep',slider_step,'max',max,'min',1,'Value',1);



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clear dd

electrode_spacing=str2num(get(handles.edit1,'String'));
number_electrodes=(str2num(get(handles.edit2,'String')));
max_n=str2num(get(handles.edit3,'String'));
max_a=str2num(get(handles.edit4,'String'));

array_type=get(handles.popupmenu1,'Value');

name_file=(get(handles.edit6,'String'));


number_lines=int32(str2num(get(handles.edit7,'String')));
line_spacing=(str2num(get(handles.edit8,'String')));

% first create electrode coordinates
xel=zeros(number_electrodes,1);
yel=zeros(number_electrodes,1);
for i=1:number_electrodes
    xel(i)=(i-1)*electrode_spacing;
    zel(i)=0;
end



pos_elec=[1:1:number_electrodes];

max_pos_elec=nchoosek(pos_elec,2);
tmp_max_mes=size(max_pos_elec);





if array_type==1 %wenner
    tmp11=1;
    for i=1:tmp_max_mes(1)
    a=max_pos_elec(i,1);
    b=max_pos_elec(i,2);
    
    rest_electrodes=setdiff(pos_elec,a);
    rest_electrodes=setdiff(rest_electrodes,b);
    
    rest_electrodes=nchoosek(rest_electrodes,2); % Now this is M and N
    s_size=size(rest_electrodes);
    

    
    for j=1:s_size(1)
        m=rest_electrodes(j,1);
        n=rest_electrodes(j,2); 
            if m-a<=max_a && m>a && n<b && m<n && abs(m-a)==abs(b-n) && abs(m-a)<=max_n && abs(m-a)==abs(m-n) && abs(m-a)==abs(b-n) 

                dd(tmp11,1)=a;
                %dd(tmp11,2)=0;
                dd(tmp11,2)=b;
                %dd(tmp11,4)=0;
                dd(tmp11,3)=rest_electrodes(j,1);
                %dd(tmp11,6)=0;
                dd(tmp11,4)=rest_electrodes(j,2);
                %dd(tmp11,8)=0;
                %dd(tmp11,9)=10;

                tmp11=tmp11+1;

            end
        
    end
    end   
    
    
    
    
    
elseif array_type==2 %sch
 
    tmp11=1;
    for i=1:tmp_max_mes(1)
    a=max_pos_elec(i,1);
    b=max_pos_elec(i,2);
    
    rest_electrodes=setdiff(pos_elec,a);
    rest_electrodes=setdiff(rest_electrodes,b);
    
    rest_electrodes=nchoosek(rest_electrodes,2); % Now this is M and N
    s_size=size(rest_electrodes);
    

    
    for j=1:s_size(1)
        m=rest_electrodes(j,1);
        n=rest_electrodes(j,2); 
            if b-a<=max_a && m>a && n<b && m<n  && abs(m-a)<=max_n && abs(n-b)<=max_n  && abs(m-a)==abs(n-b)

                dd(tmp11,1)=a;
                %dd(tmp11,2)=0;
                dd(tmp11,2)=b;
                %dd(tmp11,4)=0;
                dd(tmp11,3)=rest_electrodes(j,1);
                %dd(tmp11,6)=0;
                dd(tmp11,4)=rest_electrodes(j,2);
                %dd(tmp11,8)=0;
                %dd(tmp11,9)=10;

                tmp11=tmp11+1;

            end
        
    end
    end       
    
    
    
    
    
    
    
    
    
    
    

elseif array_type==3 %dd
    tmp11=1;
    for i=1:tmp_max_mes(1)
    a=max_pos_elec(i,1);
    b=max_pos_elec(i,2);
    
    rest_electrodes=setdiff(pos_elec,a);
    rest_electrodes=setdiff(rest_electrodes,b);
    
    rest_electrodes=nchoosek(rest_electrodes,2); % Now this is M and N
    s_size=size(rest_electrodes);
    

    
    for j=1:s_size(1)
        m=rest_electrodes(j,1);
        n=rest_electrodes(j,2); 
            if b-a<=max_a && abs(min(m,n)-max(a,b))<=max_n && abs(m-n)==abs(b-a) && (min(m,n)>max(a,b))  

                dd(tmp11,1)=a;
                %dd(tmp11,2)=0;
                dd(tmp11,2)=b;
                %dd(tmp11,4)=0;
                dd(tmp11,3)=rest_electrodes(j,1);
                %dd(tmp11,6)=0;
                dd(tmp11,4)=rest_electrodes(j,2);
                %dd(tmp11,8)=0;
                %dd(tmp11,9)=10;

                tmp11=tmp11+1;

            end
        
    end
    end
    
    
    
    
    
elseif array_type==4 %pd
    tmp11=1;
    for i=1:number_electrodes


        a=i;
        %b=max_pos_elec(i,2);

        rest_electrodes=setdiff(pos_elec,a);
        %rest_electrodes=setdiff(rest_electrodes,b);

        rest_electrodes=nchoosek(rest_electrodes,2); % Now this is M and N
        s_size=size(rest_electrodes);



        for j=1:s_size(1)
            m=rest_electrodes(j,1);
            n=rest_electrodes(j,2); 
                if   abs(min(m,n)-a)<=max_n && abs(m-n)<=max_a & (min(m,n)>a)

                    dd(tmp11,1)=a;
                    %dd(tmp11,2)=0;
                    dd(tmp11,2)=0;
                    %dd(tmp11,4)=0;
                    dd(tmp11,3)=rest_electrodes(j,1);
                    %dd(tmp11,6)=0;
                    dd(tmp11,4)=rest_electrodes(j,2);
                    %dd(tmp11,8)=0;
                    %dd(tmp11,9)=10;

                    tmp11=tmp11+1;






                end

        end


    end    
    
  
    
    
    
    
elseif array_type==5 %pp
    tmp11=1;
    for i=1:number_electrodes


        a=i;

        for j=1:number_electrodes

                if   j~=i && abs(j-i)<=max_n

                    dd(tmp11,1)=a;
                    %dd(tmp11,2)=0;
                    dd(tmp11,2)=0;
                    %dd(tmp11,4)=0;
                    dd(tmp11,3)=j;
                    %dd(tmp11,6)=0;
                    dd(tmp11,4)=0;
                    %dd(tmp11,8)=0;
                    %dd(tmp11,9)=10;

                    tmp11=tmp11+1;






                end

        end


    end        
end


% Now add coordinates
clear data        
final=dd;
%Now assign coordinates in each electrode
max_num_mes=size(final);    
max_pos_elec=final; 

tmp11=1;

     for i=1:max_num_mes(1)

         for j=1:number_electrodes

             if max_pos_elec(i,1)==pos_elec(j)
                 data(i,1)=xel(j);
                 data(i,2)=-99;
                 data(i,3)=zel(j);
             end

             if max_pos_elec(i,2)==pos_elec(j)
                 data(i,4)=xel(j);
                 data(i,5)=-99;
                 data(i,6)=zel(j);
             end       

             if max_pos_elec(i,3)==pos_elec(j)
                 data(i,7)=xel(j);
                 data(i,8)=-99;
                 data(i,9)=zel(j);
             end 

             if max_pos_elec(i,4)==pos_elec(j)
                 data(i,10)=xel(j);
                 data(i,11)=-99;
                 data(i,12)=zel(j);
             end        

         end

  

     end
     
     final=[];
     for k=1:number_lines
        yel=double((k-1))*line_spacing;
        data(:,2)=yel;
        data(:,5)=yel;
        data(:,8)=yel;
        data(:,11)=yel;
        final=[final;data];
     end
     
     
     final(:,13)=10;


save(name_file,'final','-ascii');
close (gcf);



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double
tmp=str2num(get(gcbo,'String'));
if isempty(tmp)
    msgbox('Error in Number of paraller lines. Not a number.Intial Value has been automatic entered.','Check Value','error');

        
        set(gcbo,'String','2');
  
    
end

% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double
tmp=str2num(get(gcbo,'String'));
if isempty(tmp)
    msgbox('Error in Number of paraller lines. Not a number.Intial Value has been automatic entered.','Check Value','error');

        
        set(gcbo,'String','1');
  
    
end

% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
