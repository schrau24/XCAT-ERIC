function [outputImg] = nonCartesian_CS(app, kspace, trajSorted, kspaceSorted, TVlambda)
warning('off','all')
%% Coil sensitivity calculation
disp('Calculating coil sensitivity..')
% build k matrix from all currSpokes
x = app.kx_samples*app.matrixRL;
y = app.ky_samples*app.matrixAP;
z = (app.kz_samples-app.matrixFH/2);
kloc = cat(3,y,x,z); kloc = permute(kloc, [3 1 2]);

% % create 3d nufft object and sample
% % gpuDevice(1);  % reset gpu
% obj = nufft_3d(kloc,[app.matrixAP app.matrixRL app.matrixFH],'radial',1,'gpu',0);
% 
% %reconstruction (inverse transform)
% maxit = 1; % 0 or 1 for gridding, higher values for conjugate gradient
% ref = obj.iNUFT(kspace,maxit); % plain reconstruction
% 
% [recon,b1] = adapt_array_3d(ref,size(kspace,3),1);
% % b2 = b1;
% % for ii = 1:size(b2,4)
% %     tmp = b2(:,:,:,ii);
% %     tmp(find(abs(recon(:))<0.1*max(abs(recon(:))))) = 0;
% %     b2(:,:,:,ii) = tmp;
% % end
% % b1 = b2/max(abs(b2(:)));
% % clear tmp

[b1, recon] = calculateb1Maps(app, kspace);

figure(101);
tmp = abs(b1(:,:,size(b1,3)/2,:)); tmp = tmp/max(tmp(:));
montage(tmp,'Size',[2 4]);
title('calculated coil sensitivities');
drawnow;
clear tmp b2
%% Initial guessing and TV weight determination

% Set default values for the reconstruction parameters
recoparams = struct;
recoparams.inner_iterations    = 5;
recoparams.outer_iterations    = 3;

% build k matrix
x = trajSorted(:,:,1,:)*app.matrixRL;
y = trajSorted(:,:,2,:)*app.matrixAP;
z = (trajSorted(:,:,3,:)-app.matrixFH/2);
kloc = cat(3,y,x,z); kloc = permute(kloc, [3 1 2 4]);

% % initialEstimate
% param = [];
% param.E = MCNUFFT_true3D(double(kloc),double(b1));
% param.y = double(permute(kspaceSorted, [5 1 2 3 4]));
% initialEstimate = param.E'*permute(kspaceSorted,[5 1 2 3 4]);
% initialEstimate = initialEstimate/max(abs(initialEstimate(:)));

[nx, ny, nSl, nCh] = size(b1);
kspaceSorted_slices = ifft(kspaceSorted,[],3);
shiftdim=3;
shiftdata = permute(kspaceSorted_slices,[shiftdim 1:shiftdim-1 shiftdim+1:length(size(kspaceSorted_slices))]);
pixelshift = ceil(size(shiftdata,1)/2);
s_img_size = size(shiftdata);
shiftdata = [shiftdata(pixelshift+1:end,:);shiftdata(1:pixelshift,:)];
shiftdata = reshape(shiftdata,s_img_size);
kspaceSorted_slices = ipermute(shiftdata,[shiftdim 1:shiftdim-1 shiftdim+1:length(size(shiftdata))]);

% create parpool
delete(gcp('nocreate'))
parpool();
% Calculate gridding solutions as initial guess and to properly scale
% the TV weight
initialEstimate = zeros([nx, ny, nSl, size(kloc,4)]);

unqSlcs = length(unique(squeeze(kloc(3,1,:,1))));
kloc2 = kloc(:,:,1:unqSlcs:end,:); kloc2(3,:,:) = 0;
parfevalOnAll(@warning,0,'off','all');
parfor sl = 1:nSl
    tempparam = [];
    tempparam.y = double(permute(kspaceSorted_slices(:,:,sl,:,:), [3 1 2 4 5]));
    tempparam.SG = 1;
    tempparam.E = nufft_3d_slices(kloc2,double(b1(:,:,sl,:)));
    
    initialEstimate(:,:,sl,:) = tempparam.E'*tempparam.y;
    tempparam = [];
end
delete(gcp('nocreate'))
initialEstimate = initialEstimate/max(abs(initialEstimate(:)));

% Load into local variables to avoid communication overhead
innerIterations = recoparams.inner_iterations;
outerIterations = recoparams.outer_iterations;

if app.CompressedsensingreconCheckBox.Value
    %% CS reconstruction
    disp('GRASP 4D reconstruction...');
    
    % create parpool
    delete(gcp('nocreate'))
    parpool();
    parfevalOnAll(@warning,0,'off','all');
    parfor sl = 1:nSl
        fprintf('Calculating slice %d...\n', sl);
        param = [];
        
        % permute the soft-weighting, to update!
        param.SG = 1;%double(permute(SoftWeight,[1,2,5,4,3]));
        
        param.E = nufft_3d_slices(kloc2,double(b1(:,:,sl,:)));
        param.y = double(permute(kspaceSorted_slices(:,:,sl,:,:), [3 1 2 4 5]));
        param.TV = TV_Temp3D;
        param.TVWeight = TVlambda;
        param.TVWeightRes = 0;
        param.L1Weight = 0;
        param.nite = innerIterations;
        param.display = 0;
        
        recon_cs = initialEstimate(:,:,sl,:);
        recon_cs(isnan(recon_cs)) = 0;
        
        for n = 1:outerIterations
            recon_cs = CSL1NlCg_4DRadial_SG(recon_cs,param);
        end
        data(:,:,sl,:) = abs(recon_cs);
        param = [];
    end
    outputImg = data/max(abs(data(:)));
    
    delete(gcp('nocreate'))
else
    outputImg = abs(initialEstimate);
end
