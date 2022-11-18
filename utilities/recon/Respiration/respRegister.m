function [transforms, tra_table] = respRegister(IMG, app)

%% draw ROI, perform image registration to get translations

% perform registration on a single sagittal slice chosen by user with
% user defined ROI
mag = IMG;
tmp_mag = mean(mag,4);

% use imtool3d, close figure when done and mask is saved
clear tool 
tool = imtool3D_3planes(tmp_mag);
tool.setAspectRatio([app.leftrightRes.Value app.leftrightRes.Value app.footheadRes.Value]) % set voxel size
msgbox(sprintf('1) select sagittal slice \n2) draw rectangle over area to register \n3) right-click and select poly2mask \n4) close figure'));
waitfor(tool.getHandles.fig);
h = tool.getTool;
mask = single(h.getMask(1));
%
% find the sagittal slice that was chosen
[x,y,z]=ind2sub(size(tmp_mag),find(mask));
sl = unique(x);

% mask IMG, take only that slice, and perform registration
mag_mask = squeeze(mag(sl,min(y):max(y), min(z):max(z),:));

clear transforms registered
for phase = 2:size(IMG,4)
    fixed = mag_mask(:,:,1);
    moving = mag_mask(:,:,phase);
    [T, RegisteredImage] = registerSagittalImages(moving,fixed);
    transforms(:,phase-1) = squeeze(T(3,1:2))';    % only store the translational components
    registered(:,:,phase-1) = RegisteredImage;
end
tra_table = table((1:size(transforms,2))',squeeze(transforms(1,:))',squeeze(transforms(2,:))',...
    'VariableNames',{'phase','FH','AP'});
tra_table