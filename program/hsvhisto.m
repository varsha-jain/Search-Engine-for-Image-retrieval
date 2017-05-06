srcFiles = dir('D:\varsha\database\*.jpg');
imgs=cell(length(srcFiles),1);
gr=cell(length(srcFiles),1);
hi=cell(length(srcFiles),1);
a=imread('D:\varsha\database\800.jpg');
f=rgb2hsv(a);
figure(1);
imshow(f);
H1=f(:,:,3);
S1=f(:,:,2);
V1=f(:,:,3);
bins=0:1:49;
hh1=hist(H1(:),bins);
hs1=hist(S1(:),bins);
hv1=hist(V1(:),bins);

%thresh=2.5830e+03;
thresh=100;
count=1;
x=cell(length(srcFiles),1);
hsvim=cell(length(srcFiles),1);
for i = 1 : length(srcFiles)%hist of all the images in db
     filename = strcat('D:\varsha\database\',srcFiles(i).name);
    imgs{i}=imread(filename);
    hsvim{i}=rgb2hsv(imgs{i});
RH=hsvim{i}(:,:,3);
GS=hsvim{i}(:,:,2);
BV=hsvim{i}(:,:,3);
e{i}=sqrt(sum((hh1-hist(RH(:),bins)).^2+(hs1-hist(GS(:),bins)).^2+(hv1-hist(BV(:),bins)).^2));
end
    
%hred{i}=hist(R(:),bins);
%hgreen{i}=hist(G(:),bins);
%hblue{i}=hist(B(:),bins);
%end


% threshnew=100;
%if count==2
%for i = 1 : length(srcFiles)      
%      filename = strcat('D:\varsha\database\',srcFiles(i).name);
%    imgs{i}=imread(filename);
%hsvim{i}=rgb2hsv(imgs{i});
%RH=hsvim{i}(:,:,5);
%GS=hsvim{i}(:,:,3);
%BV=hsvim{i}(:,:,5);

%e{i}=sqrt(sum((hh1-hist(RH(:),bins)).^2+(hs1-hist(GS(:),bins)).^2+(hv1-hist(BV(:),bins)).^2));

%if e{i}<=threshnew
%temp{count}=imgs{i};
%count=count+1;
%end
%end
%else
  %  break;
%end
figure(2)
[x,index]=sortrows(e');
gf=num2cell(index);
cz=1;
for i=1 : 40
    ds{cz}=imgs{gf{i}};
    cz=cz+1;
    subplot(7,8,i);
h=imshow(ds{i});
end
%dcount=(count-1)./2;
%vgcount=count-1;
%dfcount=floor(dcount);
%{figure(3)
%for i=1:vgcount
%subplot(15,16,i);
%h=imshow(temp{i});
%end
%vcount=dfcount+1;
%totalretrieved=count;
%relevantretrieved=57;
%totalimagesindb=1000;
%recall=(100*relevantretrieved)./113;
%precision=(100*relevantretrieved)./count;
