function varargout = preprocessing(varargin)
% PREPROCESSING MATLAB code for preprocessing.fig
%      PREPROCESSING, by itself, creates a new PREPROCESSING or raises the existing
%      singleton*.
%
%      H = PREPROCESSING returns the handle to a new PREPROCESSING or the handle to
%      the existing singleton*.
%
%      PREPROCESSING('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PREPROCESSING.M with the given input arguments.
%
%      PREPROCESSING('Property','Value',...) creates a new PREPROCESSING or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before preprocessing_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to preprocessing_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help preprocessing

% Last Modified by GUIDE v2.5 10-Dec-2010 13:36:44

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @preprocessing_OpeningFcn, ...
                   'gui_OutputFcn',  @preprocessing_OutputFcn, ...
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


% --- Executes just before preprocessing is made visible.
function preprocessing_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to preprocessing (see VARARGIN)

% Choose default command line output for preprocessing
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes preprocessing wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = preprocessing_OutputFcn(hObject, eventdata, handles) 
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

lagran_min=str2num(get(gcbo,'String'));

if isempty(lagran_min)
    msgbox('Error in Lagrangian Value. Not a number.Intial Value has been automatic entered.','Check Value','error');
    set(gcbo,'String','0.01');
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
lagran_max=str2num(get(gcbo,'String'));

if isempty(lagran_max)
    msgbox('Error in Lagrangian Value. Not a number.Intial Value has been automatic entered.','Check Value','error');
    set(gcbo,'String','0.1');
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


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global act
try
    
    
    save(get(handles.edit3,'String'),'act');
    close (gcf)
catch
msgbox('Something went wrong');    
    
    
end
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
global mes_in
[filename,pathname]=uigetfile('*.inv','Select Data File');
mes_in=fullfile(pathname,filename);


 


% --- Executes on button press in pushbutton.
function pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global mes_in
global act


lagrn_min=str2num(get(handles.edit1,'String'));
lagrn_max=str2num(get(handles.edit2,'String'));


if lagrn_max<lagrn_min
    msgbox('Error. Minimun Lagranginal value is larger than Maximum. Values restored');
    lagrn_min=0.001;
    lagrn_max=0.1;
    set(handles.edit2,'String','0.1');
    set(handles.edit1,'String','0.01');
end

cmap=cptcmap('GMT_seis');
cmap=cmap(end:-1:1,:); %MARIOS ADD IN


tmp=importdata(mes_in);
num_files=length(tmp);


    [misc,tin,misc,misc,misc,tmp15,misc,nitr,misc,tmp16_lag,misc,tpar,misc,tmes]=...
                                    textread(tmp{1},'%s %d %s %s %s %f %s %d %s %f %s %d %s %d',1);

  num_param=tpar;                              
                                
try
    for i=1:num_files
        tmp2(:,:,i)=importdata(tmp{i},' ',1);
    end
catch
    disp('this version does not support files with different number of measurments on each file. Future wil...');
    return
end
    tmp=tmp2(:,1,1);
    x=tmp.data(1:num_param,1,1);
    y=tmp.data(1:num_param,2,1);
    
    num_param=length(x);
    
    res_model(:,1)=tmp.data(1:num_param,3,1);
    
    sip_flag=0;
    try %if in first file, there is 10th column, then check for IP or SIP
    imag_model(:,1)=tmp.data(1:num_param,4,1);
    sip_flag=1;
    end
    
    
    for i=1:num_files
        tmp=tmp2(:,1,i);
        res_model(:,i)=tmp.data(1:num_param,3);
       if sip_flag==1; imag_model(:,i)=tmp.data(1:num_param,4);end
    end
        
    

[xi,yi]=meshgrid(unique(x),unique(y));
l_final=[];
% compeare two files per time
% figure
for i=1:num_files-1
    
    % First search for resistivity (or amplitude)
    tact(:,1)=abs(log10(res_model(:,i+1))-log10(res_model(:,i)))+eps;
%     tact(:,1)=res_model(:,i+1)./res_model(:,i);
    SP_max=max(tact(:,1));
    SP_min=min(tact(:,1));
    L1=zeros(num_param,1);
    
    for j=1:num_param
        L1(j)=log10(lagrn_min) +  ( (log10(lagrn_max) - log10(lagrn_min) ) / (log10(SP_max) - log10(SP_min)) )*(log10(tact(j,1)) - log10(SP_min));    
        L1(j)=10^(L1(j));
    end

    L1=lagrn_max-L1+lagrn_min;
    
    subplot(num_files,sip_flag+1,i);
 tact(:,1)=res_model(:,i+1)-res_model(:,i);
%     FF=TriScatteredInterp(x,y,L1);
    FF=TriScatteredInterp(x,y,tact(:,1));
    vi=FF(xi,yi);

    contourf(xi,yi,vi,17,'EdgeColor','None');
    colormap(cmap);
    colorbar
    axis equal
    xlim([min(x) max(x)]);
    ylim([min(y) max(y)]);
    caxis([-100 100]);
    % If we have phase
    if sip_flag==1
       tact(:,2)=abs(imag_model(:,i+1)-imag_model(:,i));
       SP_max=max(tact(:,2));
       SP_min=min(tact(:,2));
       L2=zeros(num_param,1);
        for j=1:num_param
            L2(j)=log10(lagrn_min) +  ( (log10(lagrn_max) - log10(lagrn_min) ) / (log10(SP_max) - log10(SP_min)) )*(log10(tact(j,2)) - log10(SP_min));    
            L2(j)=10^(L2(j));
        end       
       L2=lagrn_max-L2+lagrn_min;
       subplot(sip_flag+1,num_files,i+num_files);
       FF=TriScatteredInterp(x,y,L2);
       vi=FF(xi,yi);

       contourf(xi,yi,vi,17,'EdgeColor','None');
       colormap(cmap);
       colorbar 
%        axis equal
       xlim([min(x) max(x)]);
       ylim([min(y) max(y)]);
    end
    
    % Keep final 
    LL=zeros(num_param,1);
    if sip_flag==1
        for j=1:length(L1);
            LL(j)=min(L1(j),L2(j));
        end
    end
    
    
    if sip_flag==1
        l_final=[l_final;LL];
    else
        l_final=[l_final;L1];
    end
    
end
    
 act=zeros(num_files*num_param,   num_files*num_param);
 for i=1:length(l_final)
    act(i,i)=l_final(i);
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
