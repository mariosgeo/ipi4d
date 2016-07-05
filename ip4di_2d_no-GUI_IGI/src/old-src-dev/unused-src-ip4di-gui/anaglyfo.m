function varargout = anaglyfo(varargin)
% ANAGLYFO M-file for anaglyfo.fig
%      ANAGLYFO, by itself, creates a new ANAGLYFO or raises the existing
%      singleton*.
%
%      H = ANAGLYFO returns the handle to a new ANAGLYFO or the handle to
%      the existing singleton*.
%
%      ANAGLYFO('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ANAGLYFO.M with the given input arguments.
%
%      ANAGLYFO('Property','Value',...) creates a new ANAGLYFO or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before anaglyfo_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to anaglyfo_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help anaglyfo

% Last Modified by GUIDE v2.5 20-Oct-2009 13:16:13

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @anaglyfo_OpeningFcn, ...
                   'gui_OutputFcn',  @anaglyfo_OutputFcn, ...
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


% --- Executes just before anaglyfo is made visible.
function anaglyfo_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to anaglyfo (see VARARGIN)

% Choose default command line output for anaglyfo
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes anaglyfo wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = anaglyfo_OutputFcn(hObject, eventdata, handles) 
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


anaglyfo_data=get(handles.uitable1,'Data');
save('anaglyfo_data.mat','anaglyfo_data');

close(gcf) 

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function uitable1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uitable1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

load('orig_probe_x');
data(:,1)=orig_probe_x;
data(:,2)=0;    
columnname={'Electrode','Height'};
columneditable =  [false true];
figure
set(gcbo,'Data',data,'ColumnName', columnname,'ColumnEditable', columneditable);





% --- Executes when selected cell(s) is changed in uitable1.
function uitable1_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to uitable1 (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)
global add_x_points proceed_flag

selection = eventdata.Indices;




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

load('orig_probe_x');
data=get(gcbo,'Data');
data(:,1)=orig_probe_x;

selection = eventdata.Indices;
% if selection(1,1)==1
%     msgbox('Do not modify distances','Error','Error');
%     set(gcbo,'Data',data);    
% end

lala=get(gcbo,'Data');
lala2=lala(selection(1),2);
if isnan(lala2)
  msgbox('Only Numbers are valid','Error','Error'); 
  data(selection(1),2)=0;
  set(gcbo,'Data',data);  
end
