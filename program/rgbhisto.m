srcFiles = dir('D:\varsha\database\*.jpg');
imgs=cell(length(srcFiles),1);
gr=cell(length(srcFiles),1);
hi=cell(length(srcFiles),1);
a=imread('D:\varsha\database\800.jpg');
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
x=cell(length(srcFiles),1);
for i = 1 : length(srcFiles)%hist of all the images in db
     filename = strcat('D:\varsha\database\',srcFiles(i).name);
    imgs{i}=imread(filename);
R=imgs{i}(:,:,1);
G=imgs{i}(:,:,2);
B=imgs{i}(:,:,3);

e{i}=sqrt(sum((hred1-hist(R(:),bins)).^2+(hgreen1-hist(G(:),bins)).^2+(hblue1-hist(B(:),bins)).^2));

%if e{i}<=thresh
%temp{count}=imgs{i};
%count=count+1;
%end
%hred{i}=hist(R(:),bins);
%hgreen{i}=hist(G(:),bins);
%hblue{i}=hist(B(:),bins);
end
%x=sortrows(e');

 %threshnew=1000;
%if count==2
%for i = 1 : length(srcFiles)      
%      filename = strcat('D:\varsha\database\',srcFiles(i).name);
 %   imgs{i}=imread(filename);
%R=imgs{i}(:,:,1);
%G=imgs{i}(:,:,2);
%B=imgs{i}(:,:,3);

%e{i}=sqrt(sum((hred1-hist(R(:),bins)).^2+(hgreen1-hist(G(:),bins)).^2+(hblue1-hist(B(:),bins)).^2));
%if e{i}<=threshnew
%temp{count}=imgs{i};
%count=count+1;
%end

%end
%else
 %   break;
%end
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
%dcount=(count-1)./2;
%vgcount=count-1;
%dfcount=floor(dcount);
%figure(1)
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
