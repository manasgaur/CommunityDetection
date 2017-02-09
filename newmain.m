
disp('Denser the heat map, better is the result. This is because we are looking for overlapping community');
disp(' In this program we compare Cosine Similarity, Jaccard Similarity and Girvan Newman Model for overlapping community detection in fMRI data');
clear all;
clc;
load time_series_112.mat
mat=time_series(1).val;
mat=mat';
[r c]=size(mat);
T=140;
interval_size=2;%                                                   SUBJECT TO CHANGE
nIntervals=T-interval_size+1;%c/interval_size;
data=zeros(nIntervals,r,interval_size);
nintervals=nIntervals;
%data(1,:,:)=temp;

lb=1; ub=interval_size;
for i=1:nintervals
temp=mat(:,lb:ub);
data(i,:,:)=temp;
lb=lb+1;
ub=ub+1;
if ub>140
break;
end
end

minmean=zeros(interval_size,1);
minstdev=zeros(interval_size,1);
corrmat=zeros(interval_size,r,r);
thresh=zeros(interval_size,1);

for i=1:interval_size
temp=squeeze(data(:,:,i));
corrmat(i,:,:)=corr(temp,'type','pearson');
end

for i=1:interval_size
temp=squeeze(corrmat(i,:,:));
minmean(i,1)=min(mean(temp));
minstdev(i,1)=min(std(temp));
thresh(i,1)=minmean(i,1)-(2*minstdev(i,1));
end

 filtercormat=zeros(interval_size,r,r);
for i=1:interval_size
for j=1:r
for k=1:r
if corrmat(i,j,k)> thresh(i,1)
filtercormat(i,j,k)=corrmat(i,j,k);
end
end
end
end

% checking the data with outlier and without thresholding
svd_cormat=svd(squeeze(corrmat(1,:,:)));
svd_filtercormat=svd(squeeze(filtercormat(1,:,:)));
figure, hold on
loglog(1:112,svd_cormat,'-r*');
loglog(1:112,svd_filtercormat,'-.bs');
hold off
xlabel 'number of nodes';
ylabel 'singular values';
title 'comparing the data with and without thresholding, layer 1'

svd_cormat=svd(squeeze(corrmat(2,:,:)));
svd_filtercormat=svd(squeeze(filtercormat(2,:,:)));
figure, hold on
loglog(1:112,svd_cormat,'-r*');
loglog(1:112,svd_filtercormat,'-.bs');
hold off
xlabel 'number of nodes';
ylabel 'singular values';
title 'comparing the data with and without thresholding, layer 2'

% checking the energy the graph with outliers and without outliers
[~,eoutliers]=eigs(squeeze(corrmat(1,:,:)));
[~,ewoutliers]=eigs(squeeze(filtercormat(1,:,:)));
Energy_outliers=sum(abs(real(diag(eoutliers))));
Energy_without_outliers=sum(abs(real(diag(ewoutliers))));
display('energy of graph with thresholding');
display(Energy_outliers);
display('energy of graph without thresholding');
display(Energy_without_outliers);

% performing the cosine clustering 
clustermat=zeros(interval_size,r,r);
for i=1:interval_size
S0=zeros(r,r);S1=ones(r,r);
while norm((S1-S0),'fro')>0.005
S0=S1;
for j=1:r
U=squeeze(filtercormat(1,j,:));
for k=1:r
if j==k
continue;
else
V=squeeze(filtercormat(1,k,:));
S1(j,k)=dot(U,V)/(norm(U)*norm(V));
end
end
end
end
clustermat(i,:,:)=S1;
end
%getting the cosine similarity matrix computing the community based on
%sparse graph
cosine_graph=zeros(r,r);
for k=1:interval_size
    gph=squeeze(clustermat(k,:,:));
for i=1:r
    for j=1:r
        if gph(i,j)<=0
            cosine_graph(i,j)=0;
        else
            cosine_graph(i,j)=1;
        end
    end
end
figure, hold on
spy(cosine_graph);
xlabel 'nodes';
ylabel 'nodes';
title(sprintf('cosine similarity %d',k));
hold off;
end
G=squeeze(clustermat(1,:,:));
[community,cutvalue]=ncuts(G);
if community > cutvalue
    disp('number of communities detected by cosine similarity for layer 1')
    disp(community); disp('normalized cut values');disp(cutvalue);
else
    disp('number of communities detected by cosine similarity for layer 1')
    disp(cutvalue); disp('normalized cut values');disp(community);
end
G=squeeze(clustermat(2,:,:));
[community,cutvalue]=ncuts(G);
if community > cutvalue
    disp('number of communities detected by cosine similarity for layer 2')
    disp(community); disp('normalized cut values');disp(cutvalue);
else
    disp('number of communities detected by cosine similarity for layer 2')
    disp(cutvalue); disp('normalized cut values');disp(community);
end

% identifying the most similar node in 1 neighborhood
clusterval=zeros(interval_size,r,interval_size);
for i=1:interval_size
temp=squeeze(clustermat(i,:,:));
for j=1:size(temp,1)
[M,I]=max(temp(j,:));
clusterval(i,j,1)=I;clusterval(i,j,2)=M;
end
end

cosinecluster=[squeeze(clusterval(1,:,:)) squeeze(clusterval(2,:,:))];
% Optional plot
%figure, hold on
%scatter(1:112,cosinecluster(:,1),'d');
%scatter(1:112,cosinecluster(:,3),'filled','red');
%xlabel 'nodes';
%ylabel 'nodes';
%title 'identifying 1-nearest neighbor based on cosine distance and comparing it in different layer';
%hold off;


%clustering based on euclidean distance
euclid_distance=zeros(interval_size,r,r);
for itr=1:interval_size
tmp=squeeze(filtercormat(itr,:,:));
mtx=zeros(r,r);
S0=zeros(r,r);S1=ones(r,r);
while norm((S1-S0),'fro')> 0.005
    S0=S1;
for j=1:r
for k=1:r
if j==k
mtx(j,k)=0;
else
mtx(j,k)=pdist2(tmp(j,:),tmp(k,:));
end
end
end
for i=1:r
    rangei=max(mtx(i,:));
    for j=1:r
        S1(i,j)=mtx(i,j)/rangei;
    end
end
end
euclid_distance(itr,:,:)=S1;
end
elayer=squeeze(euclid_distance(1,:,:));
[cutvalue,community]=ncuts(elayer);
disp('number of communities detected by euclidean distance for layer 1')
disp(community); disp('normalized cut values');disp(cutvalue);
elayer=squeeze(euclid_distance(2,:,:));
[cutvalue,community]=ncuts(elayer);
disp('number of communities detected by euclidean distance for layer 2')
disp(community); disp('normalized cut values');disp(cutvalue);

euclid_cluster=zeros(interval_size,r,interval_size);
for i=1:interval_size
tmp=squeeze(euclid_distance(i,:,:));
for j=1:size(tmp,1)
[M,I]=min(tmp(j,:));
euclid_cluster(i,j,1)=I;euclid_cluster(i,j,2)=M;
end
end
euclid_data=[squeeze(euclid_cluster(1,:,1))' squeeze(euclid_cluster(2,:,1))'];
% plotting the euclid data
figure, hold on
scatter(1:112,euclid_data(:,1),'d');
scatter(1:112,euclid_data(:,2),'filled','red');
xlabel 'nodes';
ylabel 'nodes';
title 'modularity maximization based on euclidean distance';
hold off;

%getting the euclidean similarity matrix computing the community based on
%sparse graph
euclid_graph=zeros(r,r);
for k=1:interval_size
    gph=squeeze(euclid_distance(k,:,:));
for i=1:r
    med=median(gph(i,:));
    for j=1:r
        if gph(i,j)>0
            euclid_graph(i,j)=0;
        else
            euclid_graph(i,j)=1;
        end
    end
end
figure, hold on
spy(euclid_graph);
xlabel 'nodes';
ylabel 'nodes';
title(sprintf('euclidean distance graph%d',k));
hold off;
end

% calculating the jaccard similarity
jaccardmat=zeros(interval_size,r,r);
for i=1:interval_size
    temp=squeeze(filtercormat(i,:,:));
    % making the matrix into adjacency matrix
    adj_temp=zeros(r,r);
    for m=1:r
        for n=1:r
            if temp(m,n)>=0
                adj_temp(m,n)=1;
            end
        end
    end
S0=zeros(r,r);S1=ones(r,r);
while norm((S1-S0),'fro')>0.005
    S0=S1;
for j=1:r
for k=1:r
if j==k
S1(j,k)=1;
else
S1(j,k)=jaccardsimilarity(adj_temp(j,:),adj_temp(k,:));
end
end
end
end
jaccardmat(i,:,:)=S1;
end

elayer=squeeze(euclid_distance(1,:,:));
[cutvalue,community]=ncuts(elayer);
disp('number of communities detected by jaccard similarity for layer 1')
disp(community); disp('normalized cut values');disp(cutvalue);
elayer=squeeze(euclid_distance(2,:,:));
[cutvalue,evalue]=ncuts(elayer);
disp('number of communities detected by jaccard similarity for layer 2')
disp(community); disp('normalized cut values');disp(cutvalue);

% after calculating jaccard similarity, calculate 1-neighbor neighborhood
for i=1:r
    jmat=squeeze(jaccardmat(1,:,:));
[M,I]=max(jmat(i,:));
jaccard_clust(i,1)=I;jaccard_clust(i,2)=M;
end
total_jaccardmat=[jaccard_clust];
for i=1:r
jmat=squeeze(jaccardmat(2,:,:));
[M,I]=max(jmat(i,:));
jaccard_clust(i,1)=I;jaccard_clust(i,2)=M;
end
total_jaccardmat=[total_jaccardmat jaccard_clust];

%plotting the jaccard similarity clusters
figure, hold on
scatter(1:112,total_jaccardmat(:,3),'filled','red');
scatter(1:112,total_jaccardmat(:,1),'d');
xlabel 'nodes';
ylabel 'nodes';
title 'jaccard similarity based community formation'; 
hold off

% efficiency of different metric used so far
% 1. cosine similarity
efficiency=zeros(r,1);
for i=1:r
if cosinecluster(i,1)==cosinecluster(i,3)
efficiency(i,1)=1;
end
end
Percent_cosine_efficiency=(sum(efficiency)/r)*100;
%2. euclidean with 1 neighbour
 efficiency=zeros(r,1);
for i=1:r
if euclid_data(i,1)==euclid_data(i,2)
efficiency(i,1)=1;
end
end

Percent_euclid_efficiency=(sum(efficiency)/r)*100;
% 3. jaccard Similarity
efficiency=zeros(r,1);
for i=1:r
if total_jaccardmat(i,3)==total_jaccardmat(i,1)
efficiency(i,1)=1;
end
end
Percent_jaccard_efficiency=(sum(efficiency)/r)*100;

Percent_jaccard_bin_efficiency=(sum(efficiency)/r)*100;
display(Percent_cosine_efficiency);
display(Percent_euclid_efficiency);
display(Percent_jaccard_efficiency);

%modularity optimization using newman girvan model
GNM=zeros(r,r);
m=0.5*sum(sum(temp));
for i=1:r
K=sum(temp(i,:));
for j=1:r
K2=sum(temp(j,:));
GNM(i,j)=(K*K2)/m;
end
end

% modularity maximization in brain networks layer1
maxitr=1000;
kd=eye(r,r);
A=squeeze(filtercormat(1,:,:));
gamma=0.01;
layer1=zeros(maxitr,2);
i=1;
 while i< maxitr && gamma~=1
mat=(A-gamma*GNM)*kd;
Q=sum(sum(mat));
if Q<0
break;
else
layer1(i,1)=gamma;layer1(i,2)=trace(mat);
gamma=gamma+0.01;
QT(i,1)=Q;
i=i+1;
end
 end
 
% normalized cuts
[cutvalue,community]=ncuts(mat);
disp('number of communities detected by girvan newman model and kronecker delta function for layer 1')
disp(layer1(i,2)); disp('normalized cut values');disp(cutvalue);
% modularity minimization in brain networks layer2
A=squeeze(filtercormat(2,:,:));
gamma=0.01;
layer2=zeros(maxitr,2);
i=1;
 while i< maxitr && gamma~=1
mat=(A-gamma*GNM)*kd;
Q=sum(sum(mat));
if Q<0
break;
else
layer2(i,1)=gamma;layer2(i,2)=trace(mat);
gamma=gamma+0.01;
QT(i,2)=Q;
i=i+1;
end
 end
[cutvalue,community]=ncuts(mat);
disp('number of communities detected by grivan newman model and kronecker delta function for layer 2')
disp(layer2(i,2)); disp('normalized cut values');disp(cutvalue);
%finalmodularity=[modularity1' modularity2'];
hold off
figure, hold on
scatter(layer1(:,1),layer1(:,2),'s');
scatter(layer2(:,1),layer2(:,2),'filled','red');
xlabel 'gamma value';
ylabel 'number of communities';
title 'modularity optimization using girvan newman model and kronecker delta function';
hold off

%calculating interlayer graph similarity and generating one similarity
%matrix and comparing it that of wavelet correlation matrix
% matrix has been formulated from blondel et al paper
 A=squeeze(filtercormat(1,:,:));
B=squeeze(filtercormat(2,:,:));
m=size(A,1); n=size(B,1);
S=zeros(n,m); S_new=ones(n,m); % initialize S:

while norm(S_new-S,'fro')>0.001
  S=S_new;
  % do an iteration twice
  S_new=(B*S*transpose(A)+transpose(B)*S*A)/norm(B*S*transpose(A)+transpose(B)*S*A,'fro');
  S_new=(B*S_new*transpose(A)+transpose(B)*S_new*A)/norm(B*S_new*transpose(A)+transpose(B)*S_new*A,'fro');
end

S=S_new;
[cutvalue,community]=ncuts(S);
disp('number of communities detected by Blondel Similarity measure of two graphs for layer 1')
disp(community); disp('normalized cut values');disp(cutvalue);
disp('Blondel Similarity measure doesnt work for fMRI data'); 
% clustering based on hub and authorities
% first formulate adjacency matrix for defining hub and authorities
authority=zeros(r,interval_size);
hub=zeros(r,interval_size);
for k=1:interval_size
U=ones(r,1);
temp=squeeze(filtercormat(k,:,:));
adj_temp=zeros(r,r);
for i=1:r
for j=1:r
if temp(i,j)<0 || i==j
adj_temp(i,j)=0;
else
adj_temp(i,j)=1;
end
end
end
X=adj_temp;
T_old=0; T_new=1;
while (abs(T_new-T_old))>0.0005
T_old=T_new;
V=X'*U;
U=X*V;
sumV=sum(V/(norm(V,'fro')));
sumU=sum(U/(norm(U,'fro')));
T_new=sumV+sumU;
end
authority(:,k)=V;
hub(:,k)=U;
end

interlayer_adj=zeros(r,r);
for i=1:r
for j=1:r
if S(i,j)<=0
interlayer_adj(i,j)=0;
else
interlayer_adj(i,j)=1;
end
end
end
figure,hold on
spy(sparse(interlayer_adj));
xlabel 'nodes';
ylabel 'nodes';
title 'community detection interlayer similarity based on hubs and authority'
hold off

%for k=1:interval_size
%i=rank(temp);
%i=10;
%while i~=1
%xold=zeros(r,1);
%xnew=ones(r,1);
%while norm(xnew-xold,'fro')>0.01
%xold=xnew;
%tmp=temp*xnew;
%tmp=tmp/norm(tmp,'fro');
%xnew=tmp;
%end
%temp=temp-((xnew'*temp*xnew)*xnew*xnew');
%i=i-1;
%end
%evdmat(k,:,:)=temp;
%end
%for k=1:interval_size
%mtx=zeros(r,r);
%X=squeeze(evdmat(k,:,:));
%for i=1:r
%for j=1:r
%if X(i,j)<=0
%mtx(i,j)=0;
%else
%mtx(i,j)=1;
%end
%end
%end
%spy(sparse(mtx));
%end
%evdlayer=squeeze(evdmat(1,:,:));
%[cutvalue,community]=ncuts(evdlayer);
%disp('number of communities detected by evd decomposition for layer 1')
%disp(community); disp('normalized cut values');disp(cutvalue);
%evdlayer=squeeze(evdmat(2,:,:));
%[cutvalue,community]=ncuts(evdlayer);
%disp('number of communities detected by evd decomposition for layer 2')
%disp(community); disp('normalized cut values');disp(cutvalue);
%% temp end 
