function varargout = GUIimp(varargin)
% GUIIMP MATLAB code for GUIimp.fig
%      GUIIMP, by itself, creates a new GUIIMP or raises the existing
%      singleton*.
%
%      H = GUIIMP returns the handle to a new GUIIMP or the handle to
%      the existing singleton*.
%
%      GUIIMP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUIIMP.M with the given input arguments.
%
%      GUIIMP('Property','Value',...) creates a new GUIIMP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUIimp_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUIimp_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUIimp

% Last Modified by GUIDE v2.5 08-Mar-2014 23:16:42

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUIimp_OpeningFcn, ...
                   'gui_OutputFcn',  @GUIimp_OutputFcn, ...
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


% --- Executes just before GUIimp is made visible.
function GUIimp_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUIimp (see VARARGIN)

% Choose default command line output for GUIimp
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUIimp wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUIimp_OutputFcn(hObject, eventdata, handles) 
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
[filename,filepath]=uigetfile({'*.jpg;*.tif;*.png;*.gif','All Image Files';...
          '*.*','All Files' },'mytitle',...
          'D:\varsha\images\');
      if ~ischar(filename) % on cancel press you display a message of error with errordlg
    errordlg('Error!','No file selected'); % displays an error message by means of errordlg function
   return;
      end
 handles.ima=imread(strcat(filepath,filename));
 %handles.f=filename;
 %axes(handles.axes1);
 imshow(handles.ima);
 guidata(hObject,handles);

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value=get(get(handles.select,'SelectedObject'),'String');
%val=get(get(handle.select,'SelecctedObject'),'String')
switch value
    case 'RGB'
       
            
        srcFiles = dir('D:\varsha\images\*.jpg');
imgs=cell(length(srcFiles),1);
gr=cell(length(srcFiles),1);
hi=cell(length(srcFiles),1);
%a=imread(strcat(filepath,filename)); 
subplot(2,2,1);
imshow(handles.ima);
R1=handles.ima(:,:,1);
G1=handles.ima(:,:,2);
B1=handles.ima(:,:,3);
bins=0:1:255;
hred1=hist(R1(:),bins);
hgreen1=hist(G1(:),bins);
hblue1=hist(B1(:),bins);
for i = 1 : length(srcFiles)%hist of all the images in db
     filename = strcat('D:\varsha\images\',srcFiles(i).name);
    imgs{i}=imread(filename);
R=imgs{i}(:,:,1);
G=imgs{i}(:,:,2);
B=imgs{i}(:,:,3);
e{i}=sqrt(sum((hred1-hist(R(:),bins)).^2+(hgreen1-hist(G(:),bins)).^2+(hblue1-hist(B(:),bins)).^2));
end
figure(1)
[x,index]=sortrows(e');
gf=num2cell(index);
cz=1;
for i=1 :50
    ds{cz}=imgs{gf{i}};
    
    subplot(7,8,i);
h=imshow(ds{cz});
cz=cz+1;
end
    
  
    case 'HSV'
        srcFiles = dir('D:\varsha\images\*.jpg');
imgs=cell(length(srcFiles),1);
gr=cell(length(srcFiles),1);
hi=cell(length(srcFiles),1);
%a=imread('D:\varsha\database\800.jpg');
f=rgb2hsv(handles.ima);
%figure(1);
%imshow(f);
H1=f(:,:,1);
S1=f(:,:,2);
V1=f(:,:,3);
bins=0:1:49;
hh1=hist(H1(:),bins);
hs1=hist(S1(:),bins);
hv1=hist(V1(:),bins);
x=cell(length(srcFiles),1);
hsvim=cell(length(srcFiles),1);
for i = 1 : length(srcFiles)%hist of all the images in db
     filename = strcat('D:\varsha\images\',srcFiles(i).name);
    imgs{i}=imread(filename);
    hsvim{i}=rgb2hsv(imgs{i});
RH=hsvim{i}(:,:,1);
GS=hsvim{i}(:,:,2);
BV=hsvim{i}(:,:,3);
e{i}=sqrt(sum((hh1-hist(RH(:),bins)).^2+(hs1-hist(GS(:),bins)).^2+(hv1-hist(BV(:),bins)).^2));
end
figure(2)
[x,index]=sortrows(e');
gf=num2cell(index);
cz=1;
for i=1 : 50
    ds{cz}=imgs{gf{i}};
    cz=cz+1;
    subplot(7,8,i);
h=imshow(ds{i});
end
    case 'Kmeans'
    %function colorMoments = colorMoments(image)
% input: image to be analyzed and extract 2 first moments from each R,G,B
% output: 1x6 vector containing the 2 first color momenst from each R,G,B
% channel

srcFiles = dir('D:\varsha\images\*.jpg');
imgs=cell(length(srcFiles),1);
cd=cell(length(srcFiles),1);
%a=imread('D:\imagevary\imagevary\705.jpg')
[rowsA colsA numberOfColorChannelsA] = size(handles.ima);

x=rgb2gray(handles.ima);
glcm1=graycomatrix(x);
prop1=graycoprops(glcm1);
cell2=struct2cell(prop1);
f1(1,1)=cell2{1,1};
f1(1,2)=cell2{2,1};
f1(1,3)=cell2{3,1};
f1(1,4)=cell2{4,1};

%image=imread('C:\Users\tanvi\database\0.jpg');
% extract color channels
R = double(handles.ima(:, :, 1));
G = double(handles.ima(:, :, 2));
B = double(handles.ima(:, :, 3));

% compute 2 first color moments from each channel
meanRed = mean( R(:) );
stdRed  = std( R(:) );
meanGreen = mean( G(:) );
stdGreen  = std( G(:) );
meanBlue = mean( B(:) );
stdBlue  = std( B(:) );

% construct output vector
colorMom = zeros(1, 6);
colorMom(1, :) = [meanRed stdRed meanGreen stdGreen meanBlue stdBlue];

count=1;
count1=1;
count2=1;
for i=1 : length(srcFiles)
    filename = strcat('D:\varsha\images\',srcFiles(i).name);
   imgs{i}=imread(filename);
   cd{i}=imgs{i};
   [rowsB, colsB ,numberOfColorChannelsB] = size(imgs{i});
% See if lateral sizes match.
if rowsB ~= rowsA || colsA ~= colsB
% Size of B does not match A, so resize B to match A's size.
imgs{i} = imresize(imgs{i}, [rowsA colsA]);
end

R = double(imgs{i}(:, :, 1));
G = double(imgs{i}(:, :, 2));
B = double(imgs{i}(:, :, 3));

% compute 2 first color moments from each channel
meanR(i) = mean( R(:) );
stdR(i) = std( R(:) );
meanG(i) = mean( G(:) );
stdG(i)  = std( G(:) );
meanB(i) = mean( B(:) );
stdB(i)  = std( B(:) );
if meanR(i)>meanG(i)&meanR(i)>meanB(i)
    clusterred{count}=imgs{i};
    count=count+1;
else if meanG(i)>meanR(i)&meanG(i)>meanB(i)
    clustergreen{count1}=imgs{i};
    count1=count1+1;
    
    else if meanB(i)>meanR(i)&meanB(i)>meanG(i)
    clusterblue{count2}=imgs{i};
    count2=count2+1;
        end
    end
end


% construct output vector
%colorMoments(i) = zeros(i, 6);
colorMoments(i, :) = [meanR(i) stdR(i) meanG(i) stdG(i) meanB(i) stdB(i)];

end


if meanRed>meanGreen&meanRed>meanBlue
    imgred=handles.ima;
    t=1;
else if meanGreen>meanRed&meanGreen>meanBlue
    imggreen=handles.ima;
    t=2;
    else if meanBlue>meanGreen&meanBlue>meanRed
    imgblue=handles.ima;
    t=3;
        end
    end
end


% clear workspace
%clear('R', 'G', 'B', 'meanR', 'stdR', 'meanG', 'stdG', 'meanB', 'stdB');

%end

%figure(1);
%for i=1 : length(clusterred)
 %   subplot(30,20,i);
%imshow(clusterred{1,i});
%end
if t==2
for i=1 : length(clustergreen)
%i=1;
   % filename = strcat('C:\Documents and Settings\user\Desktop\database\',srcFiles(i).name);
   %imgs{i}=imread(filename);
   gray{i}=rgb2gray(clustergreen{i});
   glcma{i}=graycomatrix(gray{i});
   prop{i}=graycoprops(glcma{i});
   c1{i}=struct2cell(prop{i});
   %c1{i}=(c{i})';
   
   greenfea(i,1)=c1{1,i}{1,1};
   greenfea(i,2)=c1{1,i}{2,1};
   greenfea(i,3)=c1{1,i}{3,1};
   greenfea(i,4)=c1{1,i}{4,1};
end

%further clustering
   count3=1;
   count4=1;
   
   clustergreen1(count3)=clustergreen(1);
   green1feature(count3,:)=[greenfea(1,1) greenfea(1,2) greenfea(1,3) greenfea(1,4)];
   count3=count3+1;
   
   clustergreen2(count4)=clustergreen(2);
   green2feature(count4,:)=[greenfea(2,1) greenfea(2,2) greenfea(2,3) greenfea(2,4)];
   count4=count4+1;
   
for i=3 : length(clustergreen)
    
   e1=sqrt((greenfea(i,1)-greenfea(1,1)).^2+(greenfea(i,2)-greenfea(1,2)).^2+(greenfea(i,3)-greenfea(1,3)).^2+(greenfea(i,4)-greenfea(1,4)).^2);
   e2=sqrt((greenfea(i,1)-greenfea(2,1)).^2+(greenfea(i,2)-greenfea(2,2)).^2+(greenfea(i,3)-greenfea(2,3)).^2+(greenfea(i,4)-greenfea(2,4)).^2);
    
   if e1<e2
       clustergreen1(count3)=clustergreen(i);
       dist1(count3)=e1;
       green1feature(count3,:)=[greenfea(i,1) greenfea(i,2) greenfea(i,3) greenfea(i,4)];
       count3=count3+1;
      
   else 
       clustergreen2(count4)=clustergreen(i);
       dist2(count4)=e2;
       green2feature(count4,:)=[greenfea(i,1) greenfea(i,2) greenfea(i,3) greenfea(i,4)];
       count4=count4+1;
           
   end
end


  
else if t==1
for i=1 : length(clusterred)
%i=1;
   % filename = strcat('C:\Documents and Settings\user\Desktop\database\',srcFiles(i).name);
   %imgs{i}=imread(filename);
   gray{i}=rgb2gray(clusterred{i});
   glcma{i}=graycomatrix(gray{i});
   prop{i}=graycoprops(glcma{i});
   c1{i}=struct2cell(prop{i});
   %c1{i}=(c{i})';
   
   redfea(i,1)=c1{1,i}{1,1};
   redfea(i,2)=c1{1,i}{2,1};
   redfea(i,3)=c1{1,i}{3,1};
   redfea(i,4)=c1{1,i}{4,1}; 
end

%further clustering
   count5=1;
   count6=1;
   
   clusterred1(count5)=clusterred(1);
   red1feature(count5,:)=[redfea(1,1) redfea(1,2) redfea(1,3) redfea(1,4)];
   count5=count5+1;
   
   clusterred2(count6)=clusterred(2);
   red2feature(count6,:)=[redfea(2,1) redfea(2,2) redfea(2,3) redfea(2,4)];
   count6=count6+1;
   
for i=3 : length(clusterred)
    
   e3=sqrt((redfea(i,1)-redfea(1,1)).^2+(redfea(i,2)-redfea(1,2)).^2+(redfea(i,3)-redfea(1,3)).^2+(redfea(i,4)-redfea(1,4)).^2);
   e4=sqrt((redfea(i,1)-redfea(2,1)).^2+(redfea(i,2)-redfea(2,2)).^2+(redfea(i,3)-redfea(2,3)).^2+(redfea(i,4)-redfea(2,4)).^2);
    
   if e3<e4
      clusterred1(count5)=clusterred(i);
      red1feature(count5,:)=[redfea(i,1) redfea(i,2) redfea(i,3) redfea(i,4)];
      dist3(count5)=e3;
      count5=count5+1;
      
   else 
       clusterred2(count6)=clusterred(i);
       red2feature(count6,:)=[redfea(i,1) redfea(i,2) redfea(i,3) redfea(i,4)];
       dist4(count6)=e4;
       count6=count6+1;
           
   end
end




else if t==3
for i=1 : length(clusterblue)
%i=1;
   % filename = strcat('C:\Documents and Settings\user\Desktop\database\',srcFiles(i).name);
   %imgs{i}=imread(filename);
   gray{i}=rgb2gray(clusterblue{i});
   glcma{i}=graycomatrix(gray{i});
   prop{i}=graycoprops(glcma{i});
   c1{i}=struct2cell(prop{i});
   %c1{i}=(c{i})';
   
   bluefea(i,1)=c1{1,i}{1,1};
   bluefea(i,2)=c1{1,i}{2,1};
   bluefea(i,3)=c1{1,i}{3,1};
   bluefea(i,4)=c1{1,i}{4,1}; 
   
end

%further clustering
   count7=1;
   count8=1;
   
   clusterblue1(count7)=clusterblue(1);
   blue1feature(count7,:)=[bluefea(1,1) bluefea(1,2) bluefea(1,3) bluefea(1,4)];
   count7=count7+1;
   
   clusterblue2(count8)=clusterblue(2);
   blue2feature(count8,:)=[bluefea(2,1) bluefea(2,2) bluefea(2,3) bluefea(2,4)];
   count8=count8+1;
   
   
  for i=3 : length(clusterblue)
    
   e5=sqrt((bluefea(i,1)-bluefea(1,1)).^2+(bluefea(i,2)-bluefea(1,2)).^2+(bluefea(i,3)-bluefea(1,3)).^2+(bluefea(i,4)-bluefea(1,4)).^2);
   e6=sqrt((bluefea(i,1)-bluefea(2,1)).^2+(bluefea(i,2)-bluefea(2,2)).^2+(bluefea(i,3)-bluefea(2,3)).^2+(bluefea(i,4)-bluefea(2,4)).^2);
   
   if e5<e6
       clusterblue1(count7)=clusterblue(i);
       blue1feature(count7,:)=[bluefea(i,1) bluefea(i,2) bluefea(i,3) bluefea(i,4)];
       dist5(count7)=e5;
       count7=count7+1;
      
   else 
       clusterblue2(count8)=clusterblue(i);
       blue2feature(count8,:)=[bluefea(i,1) bluefea(i,2) bluefea(i,3) bluefea(i,4)];
       dist6(count8)=e6;
       count8=count8+1;
           
   end
end 
    end
    end
end

% Query img

if t==2
e7=sqrt((f1(1,1)-greenfea(1,1)).^2+(f1(1,2)-greenfea(1,2)).^2+(f1(1,3)-greenfea(1,3)).^2+(f1(1,4)-greenfea(1,4)).^2);
e8=sqrt((f1(1,1)-greenfea(2,1)).^2+(f1(1,2)-greenfea(2,2)).^2+(f1(1,3)-greenfea(2,3)).^2+(f1(1,4)-greenfea(2,4)).^2);

dist1=dist1';
dist2=dist2';

if e7<e8
    for i=1 : length(clustergreen1)
    finaldist(i)=sqrt((dist1(i)-e7).^2);
    end
else
    for i=1 : length(clustergreen2)
    finaldist(i)=sqrt((dist2(i)-e8).^2);
    end

end
end

if t==1
 e9=sqrt((f1(1,1)-redfea(1,1)).^2+(f1(1,2)-redfea(1,2)).^2+(f1(1,3)-redfea(1,3)).^2+(f1(1,4)-redfea(1,4)).^2);
e10=sqrt((f1(1,1)-redfea(2,1)).^2+(f1(1,2)-redfea(2,2)).^2+(f1(1,3)-redfea(2,3)).^2+(f1(1,4)-redfea(2,4)).^2);

dist3=dist3';
dist4=dist4';

if e9<e10
    for i=1 : length(clusterred1)
    finaldist(i)=sqrt((dist3(i)-e9).^2);
    end
else
    for i=1 : length(clusterred2)
    finaldist(i)=sqrt((dist4(i)-e10).^2);
    end

end
end

if t==3
e11=sqrt((f1(1,1)-bluefea(1,1)).^2+(f1(1,2)-bluefea(1,2)).^2+(f1(1,3)-bluefea(1,3)).^2+(f1(1,4)-bluefea(1,4)).^2);
e12=sqrt((f1(1,1)-bluefea(2,1)).^2+(f1(1,2)-bluefea(2,2)).^2+(f1(1,3)-bluefea(2,3)).^2+(f1(1,4)-bluefea(2,4)).^2);

dist5=dist5';
dist6=dist6';

if e11<e12
    for i=1 : length(clusterblue1)
    finaldist(i)=sqrt((dist5(i)-e11).^2);
    end
else
    for i=1 : length(clusterblue2)
    finaldist(i)=sqrt((dist6(i)-e12).^2);
    end

end
end

finaldist=finaldist';
%num=finaldist(index);
[num index]=sort(finaldist);

%green
if t==2 
    if e7<e8
        
    clustergreen1=clustergreen1';    
    figure(1);
    
for i=1 :50
    subplot(7,8,i);
%     imshow(clustergreen2{i});
    imshow(clustergreen1{index(i,1)});
end

    else if e7>e8
            
    clustergreen2=clustergreen2';    
    figure(1);
    
for i=1 :50
    subplot(7,8,i);
    imshow(clustergreen2{index(i,1)});
end
            
    end
    end
end

%red
if t==1 
    if e9<e10
        
    clusterred1=clusterred1';    
    figure(1);
    
for i=1 :50
    subplot(7,8,i);
    imshow(clusterred1{index(i,1)});
end

    else if e9>e10
            
    clusterred2=clusterred2';    
    figure(1);
    
for i=1 :50
    subplot(7,8,i);
    imshow(clusterred2{index(i,1)});
end
            
    end
    end
end

%blue
if t==3 
    if e11<e12
        
    clusterblue1=clusterblue1';    
    figure(1);
    
for i=1 :50
    subplot(7,8,i);
%     imshow(clustergreen2{i});
    imshow(clusterblue1{index(i,1)});
end

    else if e11>e12
            
    clusterblue2=clusterblue2';    
    figure(1);
    
for i=1 :50
    subplot(7,8,i);
    imshow(clusterblue2{index(i,1)});
end
            
    end
    end
end





%     [m n]=size(array);
% cen1=sum(array(i,1))/m;
% cen2=sum(array(i,2))/m;
% cen3=sum(array(i,3))/m;
% cen4=sum(array(i,4))/m;
% 
%     centroid(1,1)=cen1;
%     centroid(1,2)=cen2;
%     centroid(1,3)=cen3;
%     centroid(1,4)=cen4;
    

    
    case 'Kfcg'
        srcFiles = dir('D:\varsha\images\*.jpg');
imgs=cell(length(srcFiles),1);
gr=cell(length(srcFiles),1);
hi=cell(length(srcFiles),1);
%a=imread('D:\varsha\database\120.jpg');
R1=handles.ima(:,:,1);
G1=handles.ima(:,:,2);
B1=handles.ima(:,:,3);
bins=0:1:255;
hred1=hist(R1(:),bins);
hgreen1=hist(G1(:),bins);
hblue1=hist(B1(:),bins);
x=cell(length(srcFiles),1);
for i = 1 : length(srcFiles)%hist of all the images in db
     filename = strcat('D:\varsha\images\',srcFiles(i).name);
    imgs{i}=imread(filename);
R=imgs{i}(:,:,1);
G=imgs{i}(:,:,2);
B=imgs{i}(:,:,3);
e{i}=sqrt(sum((hred1-hist(R(:),bins)).^2+(hgreen1-hist(G(:),bins)).^2+(hblue1-hist(B(:),bins)).^2));
end
%figure(1)
[x,index]=sortrows(e');
gf=num2cell(index);
cz=1;
for i=1 :200
    ds{cz}=imgs{gf{i}};
    
    %subplot(7,8,i);
%h=imshow(ds{cz});
cz=cz+1;
end

%Aq=imread('D:\varsha\database\0.jpg');
xzq=rgb2gray(handles.ima);
%xzq=imresize(xzq,[128 128]);
myFun2= @(block_struct1) graycomatrix(block_struct1.data);
I_GLC= blockproc(xzq,[9 9],myFun2);
Ba=im2col(I_GLC,[9 9],'sliding');
[ma na]=size(I_GLC);
for i=1:ma
ta=Ba(:,i);
ra=reshape(ta,9,9);
va(i)=graycoprops(ta);  %creating a structure of properties for each subimage
end

%to convert struct to vector format

 for iq=1:ma
ca{iq}=struct2cell(va(iq));
ua{iq}=reshape((ca{iq})',1,4);% u is d feature vector to be used in kfcg
 uq(iq,1)=ua{1,iq}{1,1};
  uq(iq,2)=ua{1,iq}{1,3};

   uq(iq,3)=ua{1,iq}{1,4};

 end


vs=va';
for ig=1:ma
arrayf(ig,1)=vs(ig).Contrast;
% arrayf(ig,2)=vs(i).Correlation;
arrayf(ig,2)=vs(ig).Energy;
arrayf(ig,3)=vs(ig).Homogeneity;
end
[szrq szcq]=size(arrayf);
cenq(1,1)=sum(arrayf(:,1)/szrq);
% cenq(1,2)=sum(arrayf(:,2)/szrq);

cenq(1,2)=sum(arrayf(:,2)/szrq);
cenq(1,3)=sum(arrayf(:,3)/szrq);
codevectorsq=cenq;

countq1=1;
countq2=1;


for i=1:ma
if uq(i,1)<cenq(1,1)
clusters1{countq1,1}=uq(i,:);
countq1=countq1+1;
else
clusters2{countq2,1}=uq(i,:);
countq2=countq2+1;
end
end

[sizeclusters1m,sizeclusters1n]=size(clusters1);
for i=1:sizeclusters1n

arrayclusters1(i,1)=clusters1{1,i}(1,1);
%arrayclusters1(i,2)=clusters1{1,i}{1,2}; 
arrayclusters1(i,2)=clusters1{1,i}(1,2);
arrayclusters1(i,3)=clusters1{1,i}(1,3);
end


[numberq1,spaceq1]=size(arrayclusters1);
centroidclusters1(1,1)=sum(arrayclusters1(:,1))/numberq1;
centroidclusters1(1,2)=sum(arrayclusters1(:,2))/numberq1;

centroidclusters1(1,3)=sum(arrayclusters1(:,3))/numberq1;

[sizeclusters2m,sizeclusters2n]=size(clusters2);
for i=1:sizeclusters2n

arrayclusters2(i,1)=clusters2{1,i}(1,1);
%arraycluster2(i,2)=cluster2{1,i}{1,2}; 
arrayclusters2(i,2)=clusters2{1,i}(1,2);
arrayclusters2(i,3)=clusters2{1,i}(1,3);
end

[numberq2,spaceq2]=size(arrayclusters2);
centroidclusters2(1,1)=sum(arrayclusters2(:,1))/numberq2;
centroidclusters2(1,2)=sum(arrayclusters2(:,2))/numberq2;
centroidclusters2(1,3)=sum(arrayclusters2(:,3))/numberq2;

%for all images


qw1=1;
qw2=1;
for w = 1 : 200%hist of all the images in db
%for w = 1 : length(srcFiles)
     %filename = strcat('D:\varsha\database\',srcFiles(w).name);
    imgs{w}=ds{w};
%A=imread();
%A=imread('C:\Users\tanvi\database\80.jpg');

xz=rgb2gray(imgs{w});

myFun1= @(block_struct) graycomatrix(block_struct.data);
I_GLCM= blockproc(xz,[9 9],myFun1);
B=im2col(I_GLCM,[9 9],'sliding');
[m n]=size(I_GLCM);
for i=1:m
t=B(:,i);
r=reshape(t,9,9);
v(i)=graycoprops(t);  %creating a structure of properties for each subimage
end

%to convert struct to vector format

 for i=1:m
c{i}=struct2cell(v(i));
u{i,1}=reshape((c{i})',1,4); % u is d feature vector to be used in kfcg
prop{w,1}{i}=u{i,1};
 end
 
 v1=v';
for i=1:m
array{w,1}(i,1)=v1(i).Contrast;
% array{w,1}(i,2)=v1(i).Correlation;
array{w,1}(i,2)=v1(i).Energy;
array{w,1}(i,3)=v1(i).Homogeneity;
end


[number,space]=size(array{w,1});
centroid=sum(array{w,1})/number;

%kfcg algo 1st iteration

%codevector=u{1};
codevector=centroid;

counta1=1;
counta2=1;

[sizearraym,sizearrayn]=size(array{w,1});
for i=1:sizearraym
    if array{w,1}(i,1)<centroid(1,1)
        for j=1:3
            clustera1{qw1,1}(counta1,j)=array{w,1}(i,j);
        end
        counta1=counta1+1;
    else
        for j=1:3
            clustera2{qw2,1}(counta2,j)=array{w,1}(i,j);
        end
        counta2=counta2+1;
    end


end
%calculating centroid of cluster1


[sizeclust1,sizeclust2]=size(clustera1{qw1,1});
centroidclust1(w,1)=sum(clustera1{qw1,1}(:,1))/sizeclust1;
centroidclust1(w,2)=sum(clustera1{qw1,1}(:,2))/sizeclust1;

centroidclust1(w,3)=sum(clustera1{qw1,1}(:,3))/sizeclust1;

%calculating centroid of cluster 2
[sizeclustm2,sizeclustn2]=size(clustera2{qw2,1});
centroidclust2(w,1)=sum(clustera2{qw2,1}(:,1))/sizeclustm2;
centroidclust2(w,2)=sum(clustera2{qw2,1}(:,2))/sizeclustm2;

centroidclust2(w,3)=sum(clustera2{qw2,1}(:,3))/sizeclustm2;


qw1=qw1+1;
qw2=qw2+1;



end

for ir=1:200
    e1(ir,1)=sqrt((centroidclust1(ir,1)-centroidclusters1(1,1)).^2+(centroidclust1(ir,2)-centroidclusters1(1,2)).^2+(centroidclust1(ir,3)-centroidclusters1(1,3)).^2);
        e2(ir,1)=sqrt((centroidclust2(ir,1)-centroidclusters2(1,1)).^2+(centroidclust2(ir,2)-centroidclusters2(1,2)).^2+(centroidclust2(ir,3)-centroidclusters2(1,3)).^2);
eq(ir,1)=e1(ir,1)+e2(ir,1);
        
    
    
end
[num index]=sortrows(eq);figure(1)
for i=1:50
    subplot(7,8,i);
    imshow(imgs{index(i,1)});
end
   
%cvb=cell2mat(ds);

    case 'GLCM'
        srcFiles = dir('D:\varsha\images\*.jpg');
imgs=cell(length(srcFiles),1);
%gr=cell(length(srcFiles),1);
%hi=cell(length(srcFiles),1);
%codebook=cell(1000,1);
for w = 1 : length(srcFiles)%hist of all the images in db
%for w = 1 : length(srcFiles)
     filename = strcat('D:\varsha\images\',srcFiles(w).name);
    imgs{w}=imread(filename);
A=imread(filename);
%A=imread('C:\Users\tanvi\database\80.jpg');

xz=rgb2gray(A);
gray=graycomatrix(xz);
gp{w,1}=graycoprops(gray);
fea{w,1}=struct2cell(gp{w,1});
% f(w,1)=fea{w,1}{1,1};
% f(w,2)=fea{w,1}{2,1};
% f(w,3)=fea{w,1}{3,1};
% f(w,4)=fea{w,1}{4,1};


f(w,1)=gp{w,1}.Contrast;
f(w,2)=gp{w,1}.Correlation;
f(w,3)=gp{w,1}.Energy;
f(w,4)=gp{w,1}.Homogeneity;

end


%B=imread('D:\imagevary\imagevary\\800.jpg');

x=rgb2gray(handles.ima);
gray1=graycomatrix(x);
gp1=graycoprops(gray1);
fea1=struct2cell(gp1);
f1(1,1)=fea1{1,1};
f1(1,2)=fea1{2,1};
f1(1,3)=fea1{3,1};
f1(1,4)=fea1{4,1};


for w=1:length(srcFiles)
   
    e{w}=sqrt((f(w,1)-f1(1,1)).^2+(f(w,2)-f1(1,2)).^2+(f(w,3)-f1(1,3)).^2+(f(w,4)-f1(1,4)).^2);
end

e1=e';
[num,index]=sortrows(e1);
figure(1)
for i=1:50
    subplot(7,9,i);
    h=imshow(imgs{index(i,1)});
    
end

    case 'Block'
        srcFiles = dir('D:\varsha\images\*.jpg');
imgs=cell(length(srcFiles),1);
%gr=cell(length(srcFiles),1);
%hi=cell(length(srcFiles),1);
%codebook=cell(1000,1);
for w = 1 : 1000%hist of all the images in db
%for w = 1 : length(srcFiles)
     filename = strcat('D:\varsha\images\',srcFiles(w).name);
    imgs{w}=imread(filename);
A=imread(filename);
%A=imread('C:\Users\tanvi\database\80.jpg');

xz=rgb2gray(A);

myFun1= @(block_struct) graycomatrix(block_struct.data);
I_GLCM= blockproc(xz,[9 9],myFun1);
B=im2col(I_GLCM,[9 9],'sliding');
[m n]=size(I_GLCM);
for i=1:m
t=B(:,i);
r=reshape(t,9,9);
v(i)=graycoprops(t);  %creating a structure of properties for each subimage
end

%to convert struct to vector format

 for i=1:m
c{i}=struct2cell(v(i));
u{i}=reshape((c{i})',1,4); % u is d feature vector to be used in kfcg
prop{w,1}{i}=u{i};
 end
 
 v1=v';
for i=1:m
array{w,1}(i,1)=v1(i).Contrast;
%array(i,2)=v1(i).Correlation;
array{w,1}(i,2)=v1(i).Energy;
array{w,1}(i,3)=v1(i).Homogeneity;
end
end

%Aq=imread('C:\Users\tanvi\database\800.jpg');
xzq=rgb2gray(handles.ima);
%xzq=imresize(xzq,[128 128]);
myFun2= @(block_struct1) graycomatrix(block_struct1.data);
I_GLC= blockproc(xzq,[9 9],myFun2);
Ba=im2col(I_GLC,[9 9],'sliding');
[ma na]=size(I_GLC);
for i=1:ma
ta=Ba(:,i);
ra=reshape(ta,9,9);
va(i)=graycoprops(ta);  %creating a structure of properties for each subimage
end

%to convert struct to vector format

 for iq=1:ma
ca{iq}=struct2cell(va(iq));
ua{iq}=reshape((ca{iq})',1,4); % u is d feature vector to be used in kfcg
 end

vs=va';
for ig=1:ma
arrayf(ig,1)=vs(ig).Contrast;
%array(i,2)=v1(i).Correlation;
arrayf(ig,2)=vs(ig).Energy;
arrayf(ig,3)=vs(ig).Homogeneity;
end
[szrq szcq]=size(arrayf);
cenq(1,1)=sum(arrayf(:,1)/szrq);
cenq(1,2)=sum(arrayf(:,2)/szrq);

cenq(1,3)=sum(arrayf(:,3)/szrq);

for w=1:1000
    [szr szc]=size(array{w,1});
    cen(1,1)=sum(array{w,1}(:,1)/szr);
    cen(1,2)=sum(array{w,1}(:,2)/szr);
    cen(1,3)=sum(array{w,1}(:,3)/szr);
    
    e(w,1)=sqrt((cen(1,1)-cenq(1,1)).^2+(cen(1,2)-cenq(1,2)).^2+(cen(1,3)-cenq(1,3)).^2);

    %cen=sum(arrayf{w,1}/szr);
end

[num index]=sortrows(e);
figure(1);
for i=1:50
    subplot(7,8,i);
    imshow(imgs{index(i,1)});
end

end


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1
index_selected = get(hObject,'Value');
list = get(hObject,'String');
item_selected = list{index_selected}
handles.item=item_selected
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
