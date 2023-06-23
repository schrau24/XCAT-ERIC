function sg_signal_final = extractSG(app, kspace, isCartesian, isRadial)

nSlcs = numel(unique(app.kz_samples(1,:)));
% find center k-space indices
if isCartesian
    kyy = squeeze(app.ky_samples(1,:)) - floor(app.matrixAP/2) - 1;
    kzz = squeeze(app.kz_samples(1,:)) - floor(app.matrixFH/2) - 1;
    centerLineIdx = find(kyy==0 & kzz==0);
    SG_ksp = abs(kspace(:,centerLineIdx,:));
    ZIP = SG_ksp;
else
    centerLineIdx = 1:nSlcs:size(kspace,2);
    % sort into slices
    for sl = 1:nSlcs
       ksp(:,:,sl,:) = permute(kspace(:,sl:nSlcs:end,:),[1 2 4 3]); 
    end
    if isRadial
        ZIP = squeeze(ksp(round(size(ksp,1)/2),:,:,:));
    else % is spiral!
        ZIP = squeeze(ksp(1,:,:,:));
    end
    ZIP = permute(ZIP,[2,1,3]);
end

% Respiratory motion detection
ZIP = abs(ZIP);
% get ratio of ZIP that is post injection
postInjection = 1-(app.injectionTime.Value+25)/(app.scanTime.Value*60);
preInjection = (app.injectionTime.Value)/(app.scanTime.Value*60);
if postInjection < 0
    postInjection = 1;
    preInjection = 0;
end
n1 = ceil(size(ZIP,2)*postInjection);
n2 = floor(size(ZIP,2)*preInjection);

% STEP 1: find the coil elements with good representation of respiratory motion
% from the late enhancement spokes
[Coil,Res_Signal_Post] = MC_Step1(ZIP,n1);

% STEP 2: Estimate motion signal using PCA from the concatated coil elements
% Those coil elements were selected in the first step
[SI,corrm,Res_Signal,ZIP1] = MC_Step2(ZIP,Coil,n1,Res_Signal_Post);

% Step 3: You noticed that the signal is not flat, due to the contrast
% injection. So, now let's estimate the envelop of the signal and substract it
if preInjection > 0
    Res_Signal = MC_Step3_ems(Res_Signal,n2);
end
close all

% check to flip resp signal if expiration is on bottom
if sum(Res_Signal < 0.5) > sum(Res_Signal > 0.5)
    Res_Signal = -Res_Signal + 1;
end
figure(1);clf;
plot(Res_Signal(:)*100+220,'r','LineWidth',2)
xlabel('k_0 readout index')
ylabel('resp signal (a.u.)')
title('Respiratory Motion')
box off

sg_signal = zeros(size(app.ky_samples,2),1);
for ii=1:size(Res_Signal,1)
    if ii==1
        sg_signal(1:centerLineIdx(ii+1)-1)=Res_Signal(ii,1);
    elseif ii < size(Res_Signal,1)
        sg_signal(centerLineIdx(ii):centerLineIdx(ii+1)-1)=Res_Signal(ii,1);
    else
        sg_signal(centerLineIdx(ii):end)=Res_Signal(ii,1);
    end
end
sg_signal_final = sg_signal;


% % compare with true RespWaveTimes
% Times = 0:TR/1000:TR/1000*(length(sg_signal)-1);
% figure(5); clf;
% plot(Times,sg_signal,'k','LineWidth',2); hold on;
% h = arrayfun(@(a)xline(a,'-.r','LineWidth',2),RespWaveTimes/1000);
% xlim([-1 max(RespWaveTimes)/1000+1])
% legend('estimated resp wave', 'ground truth resp wave times')
% xlabel('scan time (s)'); ylabel('resp signal (a.u.)')