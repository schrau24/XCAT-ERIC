function [transforms, tra_table] = respRegister(IMG, app)

%% draw ROI, perform image registration to get translations

% 3D registration of only the pancreas
recon = squeeze(IMG(:,:,:,1));
recon = recon/max(recon(:));
temp = squeeze(abs(recon(:,:,round(size(IMG,3)/2))));
figure(1); clf; imagesc(temp); colormap gray, axis equal off
h = drawrectangle(gca);
pause;
maskSz = round(h.Position);
close(figure(1))

% take only the central 50% of slices
magnitude2 = abs(IMG(maskSz(2):(maskSz(2)+maskSz(4)), maskSz(1):(maskSz(1)+maskSz(3)),...
    round(0.1*size(IMG,3)):round(0.9*size(IMG,3)),:));

clear transforms registered
for phase = 2:size(IMG,4)
    fixed = magnitude2(:,:,:,1);
    moving = magnitude2(:,:,:,phase);
    [T, RegisteredImage] = registerImages3D(moving,fixed);
    transforms(:,phase-1) = squeeze(T.T(4,1:3))';    % only store the translational components
    registered(:,:,:,phase-1) = RegisteredImage;
end

% % to view the results
% a = repmat(rot90(flip(flip(fixed,1),3),-1),[1 1 1 size(IMG,4)]);
% b = cat(4,rot90(flip(flip(fixed,1),3),-1),rot90(flip(flip(registered,1),3),-1));
% View4D(cat(2,rot90(flip(flip(magnitude2,1),3),-1),a,b),1,'axisnames',...
%     {'','','moving, end-expiration, registered'}, 'FigureName','Registration Result',...
%     'FramePanelTitle','Resp Frames')

transforms = round(transforms,2);
tra_table = table((2:size(transforms,2)+1)',app.anteriorposteriorRes.Value*squeeze(transforms(1,:))',...
    app.leftrightRes.Value*squeeze(transforms(2,:))',...
    app.footheadRes.Value*squeeze(transforms(3,:))',...
    'VariableNames',{'phase','AP','RL','FH'});
% report the results
tra_table