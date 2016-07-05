function varargout = inversion_results(varargin)
% INVERSION_RESULTS M-file for inversion_results.fig
%      INVERSION_RESULTS, by itself, creates a new INVERSION_RESULTS or raises the existing
%      singleton*.
%
%      H = INVERSION_RESULTS returns the handle to a new INVERSION_RESULTS or the handle to
%      the existing singleton*.
%
%      INVERSION_RESULTS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in INVERSION_RESULTS.M with the given input arguments.
%
%      INVERSION_RESULTS('Property','Value',...) creates a new INVERSION_RESULTS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before inversion_results_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to inversion_results_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help inversion_results

% Last Modified by GUIDE v2.5 09-Dec-2010 14:44:03

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @inversion_results_OpeningFcn, ...
                   'gui_OutputFcn',  @inversion_results_OutputFcn, ...
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


% --- Executes just before inversion_results is made visible.
function inversion_results_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to inversion_results (see VARARGIN)
%clear all
% Choose default command line output for inversion_results
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes inversion_results wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = inversion_results_OutputFcn(hObject, eventdata, handles) 
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
global final

[filename,pathname]=uigetfile('*.mat','Select Inversion File');

mes_in=fullfile(pathname, filename);

if ~isempty(mes_in)
 
    final=open(mes_in);
    final=final.final;
    create_text_mesh(handles);
end


% --------------------------------------------------------------------

function uipushtool1_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global final

[filename,pathname]=uigetfile('*.mat','Select Inversion File');

mes_in=fullfile(pathname, filename);


if ~isempty(mes_in)
final=open(mes_in);
final=final.final;



create_text_mesh(handles);
end

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function info_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to info_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

clear handles.axes1
auto_contour(handles)


% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1
global map cmap axes1


map(1,:)=[0 0 128/255];
map(2,:)=[0 0 170/255];
map(3,:)=[0 0 211/255];
map(4,:)=[0 0 255/255];
map(5,:)=[0 128/255 255/255];
map(6,:)=[0 255/255 255/255];
map(7,:)=[0 192/255 128/255];
map(8,:)=[0 255/255 0];
map(9,:)=[0 128/255 0];
map(10,:)=[128/255 192/255 0];
map(11,:)=[255/255 255/255 0];
map(12,:)=[191/255 128/255 0];
map(13,:)=[255/255 128/255 0];
map(14,:)=[255/255 0 0];
map(15,:)=[211/255 0 0];
map(16,:)=[132/255 0 64/255];
map(17,:)=[96/255 0/255 96/255];
map(18,:)=[255/255 255/255 255/255];

cmap=cptcmap('GMT_seis');
cmap=cmap(end:-1:1,:); %MARIOS ADD IN

axes1=gca;



% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


%set(handles.slider1,'max',20,'min',1,'Value',1)
global itr

itr=get(handles.slider1,'Value');
set(handles.text2,'String',num2str(itr));



% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
global itr;
itr=1;

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
global final
% get(handles.slider1,'Max')
% get(handles.slider1,'Min')
% get(handles.slider1,'Value')
if isempty(final)
    max=2;
else
max=final.itn;
end
min=1;
slider_step(1)=1/(max-min);
slider_step(2)=1/(max-min);

set(gcbo,'sliderstep',slider_step,'max',max,'min',1,'Value',1)


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1
log_flag=get(handles.checkbox1,'Value');

% --- Executes during object creation, after setting all properties.
function checkbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(gcbo,'Value',1);



function auto_contour(handles)

global final
global itr XI YI min_x tmp_x max_x tmp_y max_y min_y
global map cmap

cur_m=get(handles.slider3,'Value');






if get(handles.checkbox1,'Value')==1
    
    if final.time_lapse_flag==0
        tmp=final.res_param1(:,itr);
        tmp=log10(tmp);
    else
        tmp=final.d4_res_param1(:,cur_m,itr);
        tmp=log10(tmp); 
    end

    
if final.ip_flag==1
                tmp=log10(final.all_res_model(:,itr));%remmber if is
                tmp1=1000*(final.res_param1(:,itr)-final.all_res_model(:,itr))./(final.res_param1(:,itr));
end    
if final.sip_flag==1 % Plots for DC and IP
        
    if final.time_lapse_flag==0
        tmp=log10(abs(final.res_param1(:,itr)));
        tmp1=1000*atan2(imag(final.res_param1(:,itr)),real(final.res_param1(:,itr)));   
    else
        tmp=log10(abs(final.d4_res_param1(:,cur_m,itr)));
        tmp1=1000*atan2(imag(final.d4_res_param1(:,cur_m,itr)),real(final.d4_res_param1(:,cur_m,itr)));         
    end
        
end
    
    
    
    
    
    
    
    
    elseif get(handles.checkbox1,'Value')==0
    
        
        
        
        
        
        
        if final.time_lapse_flag==0
        tmp=final.res_param1(:,itr);
    else
        tmp=final.d4_res_param1(:,cur_m,itr);
    end
end












if get(handles.popupmenu1,'value')==1
    colors='Jet';
elseif get(handles.popupmenu1,'value')==2
    colors='hsv';
elseif get(handles.popupmenu1,'value')==3
    colors='hot';
elseif get(handles.popupmenu1,'value')==4
    colors='cool';
elseif get(handles.popupmenu1,'value')==5
    colors='gray';
elseif get(handles.popupmenu1,'value')==6
    colors='bone';
elseif get(handles.popupmenu1,'value')==7
    colors='lines';
elseif get(handles.popupmenu1,'value')==8
    colors=map;    
elseif get(handles.popupmenu1,'value')==9
    colors=cmap;   
end
    
if get(handles.popupmenu2,'value')==1
    num_colors='2';
elseif get(handles.popupmenu2,'value')==2
    num_colors='3';
elseif get(handles.popupmenu2,'value')==3
    num_colors='4';
elseif get(handles.popupmenu2,'value')==4
    num_colors='8';
elseif get(handles.popupmenu2,'value')==5
    num_colors='12';
elseif get(handles.popupmenu2,'value')==6
    num_colors='16';
elseif get(handles.popupmenu2,'value')==7
    num_colors='32';
elseif get(handles.popupmenu2,'value')==8
    num_colors='64';
elseif get(handles.popupmenu2,'value')==9
    num_colors='128';
elseif get(handles.popupmenu2,'value')==10
    num_colors='256';   
end

if get(handles.popupmenu1,'value')<=7
    final_color=cellstr([colors,'(',num_colors,')']);
    final_color=char(final_color);
else
    if get(handles.popupmenu1,'value')==8 ;final_color=map(1:17,:);end
    if get(handles.popupmenu1,'value')==9 ;final_color=cmap;end
end



if get(handles.popupmenu3,'value')==1
    stepp=1;
elseif get(handles.popupmenu3,'value')==2
    stepp=2;
elseif get(handles.popupmenu3,'value')==3
    stepp=4;
elseif get(handles.popupmenu3,'value')==4
    stepp=8;
elseif get(handles.popupmenu3,'value')==5
    stepp=12;
elseif get(handles.popupmenu3,'value')==6
    stepp=16;
elseif get(handles.popupmenu3,'value')==7
    stepp=20;
end  





xnew=[min_x:tmp_x/stepp:max_x];
ynew=[0:tmp_y/stepp:max_y];

ynew=-ynew; %Anapoda giati einai arnitika
[xxx,yyy]=meshgrid(xnew,ynew);





    ZI=TriScatteredInterp(final.param_x,-final.param_y,tmp);
    ZII = ZI(xxx,yyy);
   contourf(handles.axes1,xxx,yyy,ZII,17,'EdgeColor','none');
     colorbar('peer',handles.axes1);
     colormap(final_color);
     %shading (handles.axes1,'interp');
     if final.ip_flag==0 ;title(handles.axes1,['Iteration ',num2str(itr),' RMS=> ',num2str(final.RMS(itr))]);end
     if final.ip_flag==1 ;title(handles.axes1,['Iteration ',num2str(itr),' RMS=> ',num2str(final.res_rms(itr))]);end
     axis (handles.axes1,'equal');
     xlabel(handles.axes1,'Distance (m)');
     ylabel(handles.axes1,'Depth (m)');
     xlim(handles.axes1,[min_x max_x]);
     ylim(handles.axes1,[-max_y -min_y]);
    %colorbar(inv_fig);
    drawnow;

    
    
    
if final.ip_flag==1 || final.sip_flag==1    
    ZI=TriScatteredInterp(final.param_x,-final.param_y,tmp1);
    ZII = ZI(xxx,yyy);
   contourf(handles.axes2,xxx,yyy,ZII,17,'EdgeColor','none');
     colorbar('peer',handles.axes2);
     colormap(final_color);
     %shading (handles.axes1,'interp');
     if final.ip_flag==1 ;title(handles.axes2,['Chargebility RMS=> ',num2str(final.ip_rms(itr))]);end
     if final.sip_flag==1; title(handles.axes2,['Phase RMS=> ',num2str(final.RMS(itr))]);end
     axis (handles.axes2,'equal');
     xlabel(handles.axes2,'Distance (m)');
     ylabel(handles.axes2,'Depth (m)');
     xlim(handles.axes2,[min_x max_x]);
     ylim(handles.axes2,[-max_y -min_y]);
    %colorbar(inv_fig);
    drawnow;
    
    
end


    function create_text_mesh(handles)
        

        global final XI YI max_x max_y min_x min_y tmp_x tmp_y
        

        text=cellstr(['FileName :',final.mes_in]);
        text2=cellstr(['Number of Iterations :',num2str(final.itn)]);
        text3=cellstr(['Number of Measurments :',num2str(final.num_mes)]);
        text4=cellstr(['Number of Parameters :',num2str(final.num_param)]);
% %         text5=cellstr(['Number of Boreholes :',num2str(final.num_of_bor)]);
% %         if final.borehole_no_surface_flag==0
% %             text6=cellstr(['+Surface Measurment']);
% %         elseif final.borehole_no_surface_flag==1
% %              text6=cellstr(['No Surface measurments']);
% %         end
         text7=cellstr(['RMS error :',num2str(final.RMS(final.itn))]);

        %set(slide,'Max',final.itn);
        %set(slide,'Min',1);


%         disp_text=[text;text2;text3;text4;text5;text6;text7];
        disp_text=[text;text2;text3;text4;text7];
        set(handles.info_text,'String',disp_text);


        min_x=min(final.param_x);
        max_x=max(final.param_x);

        min_y=min(final.param_y);
        max_y=max(final.param_y);


        xx=union(final.param_x, final.param_x); % take the unique x_line
        yy=union(final.param_y, final.param_y); %take the unique y_line

        tmp_x=xx(2)-xx(1);% find the x step
        tmp_y=yy(2)-yy(1);% find the y_step
        yy=-yy; %Now creative y as negative



        [XI YI]=meshgrid(xx,yy);% and now create the gridata   
        
        
        max_itr=final.itn;

        
        slider_step(1)=1/(max_itr-1);
        slider_step(2)=1/(max_itr-1);

        set(handles.slider1,'sliderstep',slider_step,'max',max_itr,'min',1,'Value',1);
        
        
        
        
        
        % get(handles.slider1,'Value')
if final.time_lapse_flag==1
    max_files=final.num_files;
    slider_step(1)=1/(max_files-1);
    slider_step(2)=1/(max_files-1);

    set(handles.slider3,'sliderstep',slider_step,'max',max_files,'min',1,'Value',1)
    
else
    max_files=1;
end
    
    

        
        
        
        
        
        
        if final.acb_flag==1
            set(handles.pushbutton2,'Visible','on');
        else
         set(handles.pushbutton2,'Visible','off');   
        end
        
        if final.inv_flag==5 || final.inv_flag==6  || final.bgr_res_flag==1
            set(handles.pushbutton6,'Visible','on');
            set(handles.pushbutton7,'Visible','on');
            set(handles.pushbutton9,'Visible','on');
        else
           set(handles.pushbutton6,'Visible','off');
           set(handles.pushbutton7,'Visible','off');
           set(handles.pushbutton9,'Visible','off');
        end
        
        


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu1 contents as cell array
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

set(gcbo,'Value',1);


% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2


% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(gcbo,'Value',8);


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
set(gcbo,'Value',6);


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global final
global itr XI YI min_x tmp_x max_x tmp_y max_y min_y

global map cmap


% if get(handles.checkbox1,'Value')==1
%     
% tmp=final.res_param2(:,itr+1);
% tmp=log10(tmp);
% elseif get(handles.checkbox1,'Value')==0
%    tmp=final.res_param2(:,itr+1);
% end


if get(handles.popupmenu1,'value')==1
    colors='Jet';
elseif get(handles.popupmenu1,'value')==2
    colors='hsv';
elseif get(handles.popupmenu1,'value')==3
    colors='hot';
elseif get(handles.popupmenu1,'value')==4
    colors='cool';
elseif get(handles.popupmenu1,'value')==5
    colors='gray';
elseif get(handles.popupmenu1,'value')==6
    colors='bone';
elseif get(handles.popupmenu1,'value')==7
    colors='lines';
elseif get(handles.popupmenu1,'value')==8
    colors=map;    
elseif get(handles.popupmenu1,'value')==9
    colors=cmap;   
end
    
if get(handles.popupmenu2,'value')==1
    num_colors='2';
elseif get(handles.popupmenu2,'value')==2
    num_colors='3';
elseif get(handles.popupmenu2,'value')==3
    num_colors='4';
elseif get(handles.popupmenu2,'value')==4
    num_colors='8';
elseif get(handles.popupmenu2,'value')==5
    num_colors='12';
elseif get(handles.popupmenu2,'value')==6
    num_colors='16';
elseif get(handles.popupmenu2,'value')==7
    num_colors='32';
elseif get(handles.popupmenu2,'value')==8
    num_colors='64';
elseif get(handles.popupmenu2,'value')==9
    num_colors='128';
elseif get(handles.popupmenu2,'value')==10
    num_colors='256';   
end

if get(handles.popupmenu1,'value')<=7
    final_color=cellstr([colors,'(',num_colors,')']);
    final_color=char(final_color);
else
    if get(handles.popupmenu1,'value')==8 ;final_color=map(1:17,:);end
    if get(handles.popupmenu1,'value')==9 ;final_color=cmap;end
end

if get(handles.popupmenu3,'value')==1
    stepp=1;
elseif get(handles.popupmenu3,'value')==2
    stepp=2;
elseif get(handles.popupmenu3,'value')==3
    stepp=4;
elseif get(handles.popupmenu3,'value')==4
    stepp=8;
elseif get(handles.popupmenu3,'value')==5
    stepp=12;
elseif get(handles.popupmenu3,'value')==6
    stepp=16;
elseif get(handles.popupmenu3,'value')==7
    stepp=20;
end  





xnew=[min_x:tmp_x/stepp:max_x];
ynew=[0:tmp_y/stepp:max_y];

ynew=-ynew; %Anapoda giati einai arnitika
[xxx,yyy]=meshgrid(xnew,ynew);
    ZI=TriScatteredInterp(final.param_x,-final.param_y,final.ACB(:,itr));
    ZII = ZI(xxx,yyy);

   contourf(handles.axes1,xxx,yyy,ZII,17,'EdgeColor','none');
     colorbar('peer',handles.axes1);
     colormap(final_color);
     %shading (handles.axes1,'interp');
     title('Lagrangian Distribution');
     axis (handles.axes1,'equal');
     xlabel(handles.axes1,'Distance (m)');
     ylabel(handles.axes1,'Depth (m)');
     xlim(handles.axes1,[min_x max_x]);
     ylim(handles.axes1,[-max_y -min_y]);
    %colorbar(inv_fig);
    drawnow;


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global final
global itr XI YI min_x tmp_x max_x tmp_y max_y min_y

global map cmap


% if get(handles.checkbox1,'Value')==1
%     
% tmp=final.res_param2(:,itr+1);
% tmp=log10(tmp);
% elseif get(handles.checkbox1,'Value')==0
%    tmp=final.res_param2(:,itr+1);
% end


if get(handles.popupmenu1,'value')==1
    colors='Jet';
elseif get(handles.popupmenu1,'value')==2
    colors='hsv';
elseif get(handles.popupmenu1,'value')==3
    colors='hot';
elseif get(handles.popupmenu1,'value')==4
    colors='cool';
elseif get(handles.popupmenu1,'value')==5
    colors='gray';
elseif get(handles.popupmenu1,'value')==6
    colors='bone';
elseif get(handles.popupmenu1,'value')==7
    colors='lines';
elseif get(handles.popupmenu1,'value')==8
    colors=map;    
elseif get(handles.popupmenu1,'value')==9
    colors=cmap;   
end
    
if get(handles.popupmenu2,'value')==1
    num_colors='2';
elseif get(handles.popupmenu2,'value')==2
    num_colors='3';
elseif get(handles.popupmenu2,'value')==3
    num_colors='4';
elseif get(handles.popupmenu2,'value')==4
    num_colors='8';
elseif get(handles.popupmenu2,'value')==5
    num_colors='12';
elseif get(handles.popupmenu2,'value')==6
    num_colors='16';
elseif get(handles.popupmenu2,'value')==7
    num_colors='32';
elseif get(handles.popupmenu2,'value')==8
    num_colors='64';
elseif get(handles.popupmenu2,'value')==9
    num_colors='128';
elseif get(handles.popupmenu2,'value')==10
    num_colors='256';   
end

if get(handles.popupmenu1,'value')<=7
    final_color=cellstr([colors,'(',num_colors,')']);
    final_color=char(final_color);
else
    if get(handles.popupmenu1,'value')==8 ;final_color=map(1:17,:);end
    if get(handles.popupmenu1,'value')==9 ;final_color=cmap;end
end

if get(handles.popupmenu3,'value')==1
    stepp=1;
elseif get(handles.popupmenu3,'value')==2
    stepp=2;
elseif get(handles.popupmenu3,'value')==3
    stepp=4;
elseif get(handles.popupmenu3,'value')==4
    stepp=8;
elseif get(handles.popupmenu3,'value')==5
    stepp=12;
elseif get(handles.popupmenu3,'value')==6
    stepp=16;
elseif get(handles.popupmenu3,'value')==7
    stepp=20;
end  





xnew=[min_x:tmp_x/stepp:max_x];
ynew=[0:tmp_y/stepp:max_y];

ynew=-ynew; %Anapoda giati einai arnitika
[xxx,yyy]=meshgrid(xnew,ynew);

    ZI=TriScatteredInterp(final.param_x,-final.param_y,real(final.resolution_1)');
    ZII = ZI(xxx,yyy);
   contourf(handles.axes1,xxx,yyy,ZII,17,'EdgeColor','none');
     colorbar('peer',handles.axes1);
     colormap(final_color);
     %shading (handles.axes1,'interp');
     title('Model Resolution');
     axis (handles.axes1,'equal');
     xlabel(handles.axes1,'Distance (m)');
     ylabel(handles.axes1,'Depth (m)');
     xlim(handles.axes1,[min_x max_x]);
     ylim(handles.axes1,[-max_y -min_y]);
    %colorbar(inv_fig);
    drawnow;




% --- Executes on button press in pushbutton4.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global final
global itr XI YI min_x tmp_x max_x tmp_y max_y min_y
global map cmap



% if get(handles.checkbox1,'Value')==1
%     
% tmp=final.res_param2(:,itr+1);
% tmp=log10(tmp);
% elseif get(handles.checkbox1,'Value')==0
%    tmp=final.res_param2(:,itr+1);
% end


if get(handles.popupmenu1,'value')==1
    colors='Jet';
elseif get(handles.popupmenu1,'value')==2
    colors='hsv';
elseif get(handles.popupmenu1,'value')==3
    colors='hot';
elseif get(handles.popupmenu1,'value')==4
    colors='cool';
elseif get(handles.popupmenu1,'value')==5
    colors='gray';
elseif get(handles.popupmenu1,'value')==6
    colors='bone';
elseif get(handles.popupmenu1,'value')==7
    colors='lines';
elseif get(handles.popupmenu1,'value')==8
    colors=map;    
elseif get(handles.popupmenu1,'value')==9
    colors=cmap;   
end
    
if get(handles.popupmenu2,'value')==1
    num_colors='2';
elseif get(handles.popupmenu2,'value')==2
    num_colors='3';
elseif get(handles.popupmenu2,'value')==3
    num_colors='4';
elseif get(handles.popupmenu2,'value')==4
    num_colors='8';
elseif get(handles.popupmenu2,'value')==5
    num_colors='12';
elseif get(handles.popupmenu2,'value')==6
    num_colors='16';
elseif get(handles.popupmenu2,'value')==7
    num_colors='32';
elseif get(handles.popupmenu2,'value')==8
    num_colors='64';
elseif get(handles.popupmenu2,'value')==9
    num_colors='128';
elseif get(handles.popupmenu2,'value')==10
    num_colors='256';   
end

if get(handles.popupmenu1,'value')<=7
    final_color=cellstr([colors,'(',num_colors,')']);
    final_color=char(final_color);
else
    if get(handles.popupmenu1,'value')==8 ;final_color=map(1:17,:);end
    if get(handles.popupmenu1,'value')==9 ;final_color=cmap;end
end

if get(handles.popupmenu3,'value')==1
    stepp=1;
elseif get(handles.popupmenu3,'value')==2
    stepp=2;
elseif get(handles.popupmenu3,'value')==3
    stepp=4;
elseif get(handles.popupmenu3,'value')==4
    stepp=8;
elseif get(handles.popupmenu3,'value')==5
    stepp=12;
elseif get(handles.popupmenu3,'value')==6
    stepp=16;
elseif get(handles.popupmenu3,'value')==7
    stepp=20;
end  





xnew=[min_x:tmp_x/stepp:max_x];
ynew=[0:tmp_y/stepp:max_y];

ynew=-ynew; %Anapoda giati einai arnitika
[xxx,yyy]=meshgrid(xnew,ynew);


tmp=(final.res_param1(:,final.itn))./final.bgr_param; % Only for the last iteration

    ZI=TriScatteredInterp(final.param_x,-final.param_y,tmp);
    ZII = ZI(xxx,yyy);
   contourf(handles.axes1,xxx,yyy,ZII,17,'EdgeColor','none');
     colorbar('peer',handles.axes1);
     colormap(final_color);
     %shading (handles.axes1,'interp');
     title('Ratio');
     axis (handles.axes1,'equal');
     xlabel(handles.axes1,'Distance (m)');
     ylabel(handles.axes1,'Depth (m)');
     xlim(handles.axes1,[min_x max_x]);
     ylim(handles.axes1,[-max_y -min_y]);
    %colorbar(inv_fig);
    drawnow;


% --- Executes on button press in pushbutton4.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global final
global itr XI YI min_x tmp_x max_x tmp_y max_y min_y
global map cmap



% if get(handles.checkbox1,'Value')==1
%     
% tmp=final.res_param2(:,itr+1);
% tmp=log10(tmp);
% elseif get(handles.checkbox1,'Value')==0
%    tmp=final.res_param2(:,itr+1);
% end


if get(handles.popupmenu1,'value')==1
    colors='Jet';
elseif get(handles.popupmenu1,'value')==2
    colors='hsv';
elseif get(handles.popupmenu1,'value')==3
    colors='hot';
elseif get(handles.popupmenu1,'value')==4
    colors='cool';
elseif get(handles.popupmenu1,'value')==5
    colors='gray';
elseif get(handles.popupmenu1,'value')==6
    colors='bone';
elseif get(handles.popupmenu1,'value')==7
    colors='lines';
elseif get(handles.popupmenu1,'value')==8
    colors=map;    
elseif get(handles.popupmenu1,'value')==9
    colors=cmap;   
end
    
if get(handles.popupmenu2,'value')==1
    num_colors='2';
elseif get(handles.popupmenu2,'value')==2
    num_colors='3';
elseif get(handles.popupmenu2,'value')==3
    num_colors='4';
elseif get(handles.popupmenu2,'value')==4
    num_colors='8';
elseif get(handles.popupmenu2,'value')==5
    num_colors='12';
elseif get(handles.popupmenu2,'value')==6
    num_colors='16';
elseif get(handles.popupmenu2,'value')==7
    num_colors='32';
elseif get(handles.popupmenu2,'value')==8
    num_colors='64';
elseif get(handles.popupmenu2,'value')==9
    num_colors='128';
elseif get(handles.popupmenu2,'value')==10
    num_colors='256';   
end

if get(handles.popupmenu1,'value')<=7
    final_color=cellstr([colors,'(',num_colors,')']);
    final_color=char(final_color);
else
    if get(handles.popupmenu1,'value')==8 ;final_color=map(1:17,:);end
    if get(handles.popupmenu1,'value')==9 ;final_color=cmap;end
end

if get(handles.popupmenu3,'value')==1
    stepp=1;
elseif get(handles.popupmenu3,'value')==2
    stepp=2;
elseif get(handles.popupmenu3,'value')==3
    stepp=4;
elseif get(handles.popupmenu3,'value')==4
    stepp=8;
elseif get(handles.popupmenu3,'value')==5
    stepp=12;
elseif get(handles.popupmenu3,'value')==6
    stepp=16;
elseif get(handles.popupmenu3,'value')==7
    stepp=20;
end  





xnew=[min_x:tmp_x/stepp:max_x];
ynew=[0:tmp_y/stepp:max_y];

ynew=-ynew; %Anapoda giati einai arnitika
[xxx,yyy]=meshgrid(xnew,ynew);


tmp=100*((- (final.bgr_param) + final.res_param2( :,(final.itn) ) )./ (final.bgr_param )); % Only for the last iteration

    ZI=TriScatteredInterp(final.param_x,-final.param_y,tmp);
    ZII = ZI(xxx,yyy);
   contourf(handles.axes1,xxx,yyy,ZII,17,'EdgeColor','none');
     colorbar('peer',handles.axes1);
     colormap(final_color);
     %shading (handles.axes1,'interp');
     title('% Difference');
     axis (handles.axes1,'equal');
     xlabel(handles.axes1,'Distance (m)');
     ylabel(handles.axes1,'Depth (m)');
     xlim(handles.axes1,[min_x max_x]);
     ylim(handles.axes1,[-max_y -min_y]);
    %colorbar(inv_fig);
    drawnow;



% --- Executes during object creation, after setting all properties.
function pushbutton7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --------------------------------------------------------------------
function Untitled_4_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_5_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%gcf=handles.axes1;
set(gcf, 'PaperPositionMode', 'auto')   % Use screen size
print -djpeg  -r300 tst.jpg

% --------------------------------------------------------------------
function Untitled_6_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global final itr


if final.time_lapse_flag==0 
    if final.dc_flag==1
        tmp=[final.param_x -final.param_y final.res_param1(:,final.itn) log10(final.res_param1(:,final.itn))];
    elseif final.ip_flag==1
        tmp=[final.param_x -final.param_y final.all_res_model(:,final.itn) log10(final.all_res_model(:,final.itn)) final.chargeb];
    elseif final.sip_flag==1
        tmp=[final.param_x -final.param_y abs(final.res_param1(:,final.itn))...
            abs(final.res_param1(:,final.itn) ) ...
            1000*atan2(  imag(final.res_param1(:,final.itn)),  real(final.res_param1(:,final.itn))  )   ];
    end
else
    if final.dc_flag==1
        tmp=[final.param_x -final.param_y final.d4_res_param1(:,:,final.itn) log10(final.d4_res_param1(:,:,final.itn))];
    elseif final.ip_flag==1
%         tmp=[final.param_x;-final_param_y;final.all_res_model(:,itn);log10(final.all_res_model(:,itn)),final.chargeb];
    elseif final.sip_flag==1
        tmp=[final.param_x -final.param_y abs(final.d4_res_param(:,:,final.itn)) log10(abs((final.d4_res_param(:,:,final.itn)))) 1000*atan(imag(final.d4_res_param(:,:,final.itn)),real(final.d4_res_param(:,:,final.itn)))];
    end
    
end

save('results.txt','tmp','-ascii');




% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
colormapeditor;


% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global final
global itr XI YI min_x tmp_x max_x tmp_y max_y min_y
global map cmap



% if get(handles.checkbox1,'Value')==1
%     
% tmp=final.res_param2(:,itr+1);
% tmp=log10(tmp);
% elseif get(handles.checkbox1,'Value')==0
%    tmp=final.res_param2(:,itr+1);
% end


if get(handles.popupmenu1,'value')==1
    colors='Jet';
elseif get(handles.popupmenu1,'value')==2
    colors='hsv';
elseif get(handles.popupmenu1,'value')==3
    colors='hot';
elseif get(handles.popupmenu1,'value')==4
    colors='cool';
elseif get(handles.popupmenu1,'value')==5
    colors='gray';
elseif get(handles.popupmenu1,'value')==6
    colors='bone';
elseif get(handles.popupmenu1,'value')==7
    colors='lines';
elseif get(handles.popupmenu1,'value')==8
    colors=map;    
elseif get(handles.popupmenu1,'value')==9
    colors=cmap;   
end
    
if get(handles.popupmenu2,'value')==1
    num_colors='2';
elseif get(handles.popupmenu2,'value')==2
    num_colors='3';
elseif get(handles.popupmenu2,'value')==3
    num_colors='4';
elseif get(handles.popupmenu2,'value')==4
    num_colors='8';
elseif get(handles.popupmenu2,'value')==5
    num_colors='12';
elseif get(handles.popupmenu2,'value')==6
    num_colors='16';
elseif get(handles.popupmenu2,'value')==7
    num_colors='32';
elseif get(handles.popupmenu2,'value')==8
    num_colors='64';
elseif get(handles.popupmenu2,'value')==9
    num_colors='128';
elseif get(handles.popupmenu2,'value')==10
    num_colors='256';   
end

if get(handles.popupmenu1,'value')<=7
    final_color=cellstr([colors,'(',num_colors,')']);
    final_color=char(final_color);
else
    if get(handles.popupmenu1,'value')==8 ;final_color=map(1:17,:);end
    if get(handles.popupmenu1,'value')==9 ;final_color=cmap;end
end

if get(handles.popupmenu3,'value')==1
    stepp=1;
elseif get(handles.popupmenu3,'value')==2
    stepp=2;
elseif get(handles.popupmenu3,'value')==3
    stepp=4;
elseif get(handles.popupmenu3,'value')==4
    stepp=8;
elseif get(handles.popupmenu3,'value')==5
    stepp=12;
elseif get(handles.popupmenu3,'value')==6
    stepp=16;
elseif get(handles.popupmenu3,'value')==7
    stepp=20;
end  





xnew=[min_x:tmp_x/stepp:max_x];
ynew=[0:tmp_y/stepp:max_y];

ynew=-ynew; %Anapoda giati einai arnitika
[xxx,yyy]=meshgrid(xnew,ynew);




    ZI=TriScatteredInterp(final.param_x,-final.param_y,final.bgr_param);
    ZII = ZI(xxx,yyy);
   contourf(handles.axes1,xxx,yyy,ZII,17,'EdgeColor','none');
     colorbar('peer',handles.axes1);
     colormap(final_color);
     %shading (handles.axes1,'interp');
     title('Background');
     axis (handles.axes1,'equal');
     xlabel(handles.axes1,'Distance (m)');
     ylabel(handles.axes1,'Depth (m)');
     xlim(handles.axes1,[min_x max_x]);
     ylim(handles.axes1,[-max_y -min_y]);
    %colorbar(inv_fig);
    drawnow;


% --- Executes during object creation, after setting all properties.
function pushbutton9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on slider movement.
function slider3_Callback(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global cur_m

cur_m=get(handles.slider3,'Value');
set(handles.text5,'String',num2str(cur_m));


% --- Executes during object creation, after setting all properties.
function slider3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


global final
% get(handles.slider1,'Max')
% get(handles.slider1,'Min')
% get(handles.slider1,'Value')
if isempty(final)
    max_files=1;
else
    if final.time_lapse_flag==1 
        max_files=final.num_files;
    else
        max_files=1;
    end
end
min_files=1;
slider_step(1)=1/(max_files-min_files);
slider_step(2)=1/(max_files-min_files);
try
set(gcbo,'sliderstep',slider_step,'max',max_files,'min',1,'Value',1)
end

% --- Executes during object creation, after setting all properties.
function text3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function text5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
