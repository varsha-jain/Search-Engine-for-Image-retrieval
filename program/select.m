function varargout = select(varargin)
% SELECT MATLAB code for select.fig
%      SELECT, by itself, creates a new SELECT or raises the existing
%      singleton*.
%
%      H = SELECT returns the handle to a new SELECT or the handle to
%      the existing singleton*.
%
%      SELECT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SELECT.M with the given input arguments.
%
%      SELECT('Property','Value',...) creates a new SELECT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before select_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to select_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help select

% Last Modified by GUIDE v2.5 25-Jan-2014 20:09:21

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @select_OpeningFcn, ...
                   'gui_OutputFcn',  @select_OutputFcn, ...
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
end
% End initialization code - DO NOT EDIT


% --- Executes just before select is made visible.
function select_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to select (see VARARGIN)

% Choose default command line output for select
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
end



% UIWAIT makes select wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = select_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;end






% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename,filepath]=uigetfile({'*.jpg;*.tif;*.png;*.gif','All Image Files';...
          '*.*','All Files' },'mytitle',...
          'D:\varsha\database\');
      if ~ischar(filename) % on cancel press you display a message of error with errordlg
    errordlg('Error!','No file selected'); % displays an error message by means of errordlg function
   return;
      end
 image=imread(strcat(filepath,filename)); 
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 srcFiles = dir('D:\varsha\database\*.jpg');
imgs=cell(length(srcFiles),1);
gr=cell(length(srcFiles),1);
hi=cell(length(srcFiles),1);
a=imread(strcat(filepath,filename)); 
R1=a(:,:,1);
G1=a(:,:,2);
B1=a(:,:,3);

bins=0:1:255;

hred1=hist(R1(:),bins);
hgreen1=hist(G1(:),bins);
hblue1=hist(B1(:),bins);

%thresh=2.5830e+03;
thresh=30;
count=1;

for i = 1 : length(srcFiles)%hist of all the images in db
     filename = strcat('D:\varsha\database\',srcFiles(i).name);
    imgs{i}=imread(filename);
R=imgs{i}(:,:,1);
G=imgs{i}(:,:,2);
B=imgs{i}(:,:,3);

e{i}=sqrt(sum((hred1-hist(R(:),bins)).^2+(hgreen1-hist(G(:),bins)).^2+(hblue1-hist(B(:),bins)).^2));
if e{i}<=thresh
temp{count}=imgs{i};
count=count+1;
end

%hred{i}=hist(R(:),bins);
%hgreen{i}=hist(G(:),bins);
%hblue{i}=hist(B(:),bins);
end
 threshnew=1000;
if count==2
for i = 1 : length(srcFiles)      
      filename = strcat('D:\varsha\database\',srcFiles(i).name);
    imgs{i}=imread(filename);
R=imgs{i}(:,:,1);
G=imgs{i}(:,:,2);
B=imgs{i}(:,:,3);

e{i}=sqrt(sum((hred1-hist(R(:),bins)).^2+(hgreen1-hist(G(:),bins)).^2+(hblue1-hist(B(:),bins)).^2));
if e{i}<=threshnew
temp{count}=imgs{i};
count=count+1;
end

end
else
    return;
end

%for i = 1 : length(srcFiles)
%e{i}=sqrt(sum(hred1-hred{i}).^2+(hgreen1-hgreen{i}).^2+(hblue1-hblue{i}).^2);
%if e{i}<=thresh
%t{count}=imgs{i};
%count=count+1;
%end
%end

    

x=cell(length(srcFiles),1);

x=sortrows(e');
dcount=(count-1)./2;
vgcount=count-1;
dfcount=floor(dcount);
figure(1)
for i=1:vgcount
subplot(15,16,i);
h=imshow(temp{i});
end
vcount=dfcount+1;
totalretrieved=count;
relevantretrieved=57;
totalimagesindb=1000;
recall=(100*relevantretrieved)./113;
precision=(100*relevantretrieved)./count;end


 
 
 
 
 
 
