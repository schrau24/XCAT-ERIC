function [outputImg] = nonCartesian_CS(app, kspace, trajSorted, kspaceSorted, TVlambda)
warning('off','all')
%% Coil sensitivity calculation
disp('Calculating coil sensitivity..')

[b1, recon] = calculateb1Maps(app, kspace);

figure(101);
tmp = abs(b1(:,:,round(size(b1,3)/2),:)); tmp = tmp/max(tmp(:));
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
x = trajSorted(:,:,1,:);
y = trajSorted(:,:,2,:);
z = trajSorted(:,:,3,:);
kloc = cat(3,y,x); kloc = permute(kloc, [3 1 2 4]);

nSlcs = length(unique(squeeze(z(1,:,1,1))));
nRO = size(kspaceSorted,2)/nSlcs;
kk = zeros(size(kspaceSorted,1),nRO,size(b1,4),size(kspaceSorted,4),app.matrixFH);
for sl = 1:nSlcs
   kk(:,:,:,:,sl) = kspaceSorted(:,sl:nSlcs:end,:,:);
end
kspaceSorted2 = permute(kk,[1 2 5 3 4]); clear kk;
[~,kz_ord] = sort(app.kz_samples(1,1:nSlcs));
kspaceSorted2 = kspaceSorted2(:,:,kz_ord,:,:);
% pad to correct size
kspaceSorted2 = padarray(kspaceSorted2,[0 0 app.matrixFH-nSlcs 0 0],'pre');

[nx, ny, nSl, nCh] = size(b1);
kspaceSorted_slices = ifft(kspaceSorted2,[],3);
shiftdim=3;
shiftdata = permute(kspaceSorted_slices,[shiftdim 1:shiftdim-1 shiftdim+1:length(size(kspaceSorted_slices))]);
pixelshift = ceil(size(shiftdata,1)/2);
s_img_size = size(shiftdata);
shiftdata = [shiftdata(pixelshift+1:end,:);shiftdata(1:pixelshift,:)];
shiftdata = reshape(shiftdata,s_img_size);
kspaceSorted_slices = ipermute(shiftdata,[shiftdim 1:shiftdim-1 shiftdim+1:length(size(shiftdata))]);
clear kspaceSorted2

% Calculate gridding solutions as initial guess and to properly scale
% the TV weight
initialEstimate = zeros([nx, ny, nSl, size(kloc,4)]);
kloc = kloc(:,:,1:nSlcs:end,:);
X = squeeze(kloc(1,:,:,:));
Y = squeeze(kloc(2,:,:,:));
Traj_nt = X+1i*Y;
DensityComp_nt = sqrt(X.^2+Y.^2);
kdata_nt = kspaceSorted_slices.*repmat(sqrt(permute(DensityComp_nt,[1 2 4,5,3])),[1,1,app.matrixFH,size(b1,4),1]);
% create parpool
delete(gcp('nocreate'))
parpool();
parfevalOnAll(@warning,0,'off','all');
parfor sl = 1:nSl
    tempparam = [];
    tempparam.y = double(kdata_nt(:,:,sl,:,:));
    tempparam.SG = 1;
    tempparam.E = MCNUFFT3D(double(Traj_nt),double(DensityComp_nt),double(b1(:,:,sl,:)));
    
    initialEstimate(:,:,sl,:) = tempparam.E'*tempparam.y;
    tempparam = [];
end
delete(gcp('nocreate'))

if app.CompressedsensingreconCheckBox.Value
    %% CS reconstruction
    % Load into local variables to avoid communication overhead
    globalTVWeight = max(abs(initialEstimate(:)))*TVlambda;
    innerIterations = recoparams.inner_iterations;
    outerIterations = recoparams.outer_iterations;
    disp('GRASP 4D reconstruction...');
    
    % create parpool
    delete(gcp('nocreate'))
    parpool();
    parfevalOnAll(@warning,0,'off','all');
    parfor sl = 1:nSl
        fprintf('Calculating slice %d...\n', sl);
        param = [];
        
        % permute the soft-weighting
        param.SG = 1;%double(permute(SoftWeight,[1,2,5,4,3]));
        
        % Limit data for MCNUFFT3D operator to one slice for enabling
        % multithreaded calculation of slices.
        param.y = double(kdata_nt(:,:,sl,:,:));
        
        % Prepare the operator for the forward operation
        param.E = MCNUFFT3D(double(Traj_nt),double(DensityComp_nt),double(b1(:,:,sl,:)));
        param.TV = TV_Temp3D;
        param.TVWeight = globalTVWeight;
        param.TVWeightRes = 0;
        param.L1Weight = 0;
        param.nite = innerIterations;
        param.display = 1;
        param.slice = sl;
        
        recon_cs = initialEstimate(:,:,sl,:);
        
        for n = 1:outerIterations
            recon_cs = CSL1NlCg_4DRadial_SG(recon_cs,param);
        end
        data(:,:,sl,:) = abs(recon_cs);
        param = [];
    end
    outputImg = data/max(abs(data(:)));
    
    delete(gcp('nocreate'))
else
    outputImg = abs(initialEstimate/max(abs(initialEstimate(:))));
end
