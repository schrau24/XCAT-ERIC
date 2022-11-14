%Motion detection: step 2

%Li Feng, NYU, 12/18/2017

function [SI,corrm,Res_Signal,ZIP1] = MC_Step2(ZIP,Coil,n1,Res_Signal_Post);

%Do PCA along all the selected coils concatated together
ntviews = size(ZIP,2);
ZIP1 = ZIP(:,:,Coil);
SI = permute(ZIP1,[1,3,2]); %Concatate the selected coils
SI = abs(reshape(SI,[size(SI,1)*size(SI,2),ntviews])');
covariance = cov(SI);
[PC, V] = eig(covariance);
V = diag(V);
[junk, rindices] = sort(-1*V);
V = V(rindices);
PC = PC(:,rindices);
SI = (PC' * SI')';
SI = SI(:,1:5);%Consider the first 5 principal component only

% Do some smoothing
clear SI_STD
for ii = 1:size(SI,2)
    SI(:,ii) = smooth(SI(:,ii),6,'lowess');
end

% check the correlation with the motion signal detected in the first step
% calculate the correlation for the late phase only
for ch = 1:size(SI,2)
    corrm(ch) = xcov(Res_Signal_Post,SI(end-n1+1:end,ch),0,'coef');
end
corrm_abs = abs(corrm);
Coil_Index_PCA = find(corrm_abs == max(corrm_abs))
Res_Signal = SI(:,Coil_Index_PCA);

%flip the signal if the corrlation coefficient is negative
if corrm(Coil_Index_PCA)~=corrm_abs(Coil_Index_PCA);
    Res_Signal = -Res_Signal;
end

% % if the first couple of points on Res_Signal are absurdly high, set to the
% % mean
% ind = find(Res_Signal>mean(Res_Signal(end-n1+1:end))+3*std(Res_Signal(end-n1+1:end)));
% if ~isempty(find(ind<10))
%     Res_Signal(ind(ind<10)) = mean(Res_Signal(end-n1+1:end));
%     
% else
%     % if the first couple of points on Res_Signal are absurdly low, set to the
%     % mean
%     ind = find(Res_Signal<mean(Res_Signal(end-n1+1:end))-3*std(Res_Signal(end-n1+1:end)));
%     Res_Signal(ind(ind<10)) = mean(Res_Signal(end-n1+1:end));
% end

Res_Signal = Res_Signal-min(Res_Signal(:));
Res_Signal = Res_Signal./max(Res_Signal(:));