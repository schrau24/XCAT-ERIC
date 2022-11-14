function [Res_Signal2]=MC_Step3_ems(Res_Signal, n2)
% Motion detection: step 3

% find the pre- and post-contrast peaks, subtract the mean of peaks for each
[~, Peak_Index]=findpeaks(double(Res_Signal),'MinPeakProminence',0.05);
figure,plot(Res_Signal);
hold on;
plot(Peak_Index,Res_Signal(Peak_Index),'ro');

if length(Peak_Index) > 1
    
    % find peaks before injection
    pkDiffsOutliers = Peak_Index(isoutlier(abs(diff(Res_Signal(Peak_Index)))));
    preInd = Peak_Index(Peak_Index<=n2);
    
    if isempty(preInd)
        preInd = Peak_Index(1:2);
    end
    if length(preInd) == 1 % add the next peak as well
        preInd(2) = Peak_Index(find(Peak_Index>preInd,1,'first'));
    end
    
    preRes = Res_Signal(1:preInd(end)-1) - mean(Res_Signal(1:preInd(end)-1));
    preRes = preRes/max(abs(preRes(:)));

    if ~isempty(pkDiffsOutliers)
        postInd = Peak_Index(find(Peak_Index>pkDiffsOutliers(end),1,'first'));
    else
        postInd = Peak_Index(find(Peak_Index>max(preInd)));
    end
    postRes = Res_Signal(postInd:end);
    p = polyfit(1:length(postRes),postRes,2);
    postRes = postRes - polyval(p,1:length(postRes))';
    postRes = postRes/max(abs(postRes(:)));
    
    % during injection, fit
    midInjection = Res_Signal(preInd(end):postInd-1);
    p = polyfit(1:length(midInjection),midInjection,1);
    midInjection = midInjection - polyval(p,1:length(midInjection))';
    midInjection = midInjection/max(abs(midInjection(:)));
    
    Res_Signal2 = cat(1,preRes,midInjection,postRes) + 1;
    
end

% figure; plot(Res_Signal2)