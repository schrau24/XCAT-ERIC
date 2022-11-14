%Motion detection: step 1

%Li Feng, NYU, 12/18/2017

function [Coil,Res_Signal_Post]=MC_Step1(ZIP,n1);

ntviews=size(ZIP,2);
close all
k=1;N=n1;
clear Nav

% Do PCA along each coil element and to get the first two PCs from each
% coil, as described in the paper
for ii=1:size(ZIP,3)
    SI=ZIP(:,end-N+1:end,ii)';
    covariance=cov(SI);
    [PC, V] = eig(covariance);
    V = diag(V);
    [junk, rindices] = sort(-1*V);
    V = V(rindices);
    PC = PC(:,rindices);
    SI = (PC' * SI')';
    for jj=1:3
        tmp=smooth(SI(:,jj),6,'lowess');
        tmp=tmp-min(tmp(:));
        tmp=tmp./max(tmp(:));
        Nav(:,k)=tmp;k=k+1;
    end
end

% coil clustering
% code obtained from Tao Zhang, stanford
thresh = 0.97;
[Res_Signal, cluster] = CoilClustering(Nav, thresh);
tmp=find(cluster~=0);
while length(tmp)<=2
    thresh=thresh-0.01;
    [Res_Signal, cluster] = CoilClustering(Nav, thresh);
    tmp=find(cluster~=0);
end
thresh
Res_Signal=Res_Signal-min(Res_Signal(:));
Res_Signal_Post=Res_Signal./max(Res_Signal(:));

% find the "good" coil elements used for estimation of all motion later
cluster=abs(reshape(cluster,[3,size(Nav,2)/3]));
cluster=sum(cluster,1);
Coil=find(cluster>0);

% if mean(Res_Signal_Post)<0.5
%     Res_Signal_Post=-Res_Signal_Post;
%     Res_Signal_Post=Res_Signal_Post-min(Res_Signal_Post(:));
%     Res_Signal_Post=Res_Signal_Post./max(Res_Signal_Post(:));
% end

%figure,imagesc(abs(ZIP(:,end-n1+1:end,1))),axis image, axis off, colormap(gray),title('Respiratory Motion')
%hold on
%plot(-Res_Signal_Post(:)*100+220,'r')
