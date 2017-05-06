A=imread('C:\Documents and Settings\AS\My Documents\Downloads\z1.jpg');
A=rgb2gray(A);
myFun1= @(block_struct) graycomatrix(block_struct.data);
I_GLCM= blockproc(A,[3 3],myFun1);
B=im2col(I_GLCM,[3 3],'sliding');
[m n]=size(I_GLCM);
for i=1:m
t=B(:,i);
r=reshape(t,3,3);
v(i)=graycoprops(t);  %creating a structure of properties for each subimage
end

%to convert struct to vector format

 for i=1:m
c{i}=struct2cell(v(i));
u{i}=reshape((c{i})',1,4); % u is d feature vector to be used in kfcg
end

%to convert structure to array for computing centroid
v1=v';
for i=1:m
array(i,1)=v1(i).Contrast;
array(i,2)=v1(i).Correlation;
array(i,3)=v1(i).Energy;
array(i,4)=v1(i).Homogeneity;
end

%calculate centroid=  sum_of(vectors)/number_of(vectors)

 [number,space]=size(array);
centroid=sum(array)/number;

%kfcg algo 1st iteration

%codevector=u{1};
codevector=centroid;
count1=1;
count2=1;
%for i=2:m
%if u{i}{1}<u{1}{1}

for i=1:m
if u{i}{1}<centroid(1,1)
cluster1{count1}=u{i};
count1=count1+1;
else
cluster2{count2}=u{i};
count2=count2+1;
end
end


 

