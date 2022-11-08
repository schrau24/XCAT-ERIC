function [kspace] = sampleKSpace(app, isCartesian, SNR)
warning('off','all')

FOV = size(app.IMG_CP,1:3);
if isCartesian
    IMGs = app.IMG_CP;
else
    % pad matrix to sample on a square matrix in RL/AP
    IMGs = padarray(app.IMG_CP,[(max(FOV(1:2))-FOV(1))/2 ...
        (max(FOV(1:2))-FOV(2))/2 0 0]);
    FOV = size(IMGs,1:3);
end

% Calculate simulated coils in 3D
nCh = 8;
Coils = Simulate_Coils(FOV,nCh);

% also load typical slice profile
load([app.appPath 'utilities/sampling/SliceProfile.mat']);
SP = SP';
% interp slice profile to number of slices
if FOV(3) ~= 35
    SP = resample(SP,FOV(3),35);
    % ensure it is symmetric
    SP(round(FOV(3)/2+1:end)) = SP(floor(FOV(3)/2):-1:1);
end

% using timing and TR, we need to sort our readouts to determine what
% contrast phase they belong to
nFE = size(app.kx_samples,2);
TR = app.TR_sim.Value;
readoutTiming = 0:TR:(nFE-1)*TR;     % in ms
nPhases = size(IMGs,4);
sortedROs = sortReadOuts(nPhases, app.timing*1000, nFE, readoutTiming);

%% create and save simulated k-space
kspace = complex(zeros([size(app.kx_samples,1:2),nCh]));
h = waitbar(0,'sampling k-space');
% gpuDevice(1);   % reset gpu
for currPhase = 1:nPhases
    
    if mod(currPhase,2) == 1
        waitbar(currPhase/nPhases,h)
    end
    currSpokes = find(sortedROs==currPhase);
    
    if ~isempty(currSpokes)
        % Grab the current RespPhase and RRPhase image
        IMG = IMGs(:,:,:,currPhase);
        
        % two types of sampling:
        % 1. Cartesian/Pseudo-spiral, simply use ky,kz points
        % 2. Radial/Spiral, use nufft for kx,ky,kz points
        if isCartesian
            % loop over coils and extract kspace lines
            for coil = 1:nCh
                tmp = Coils(:,:,:,coil).*(IMG).*(repmat(permute(SP, [2 3 1]),[FOV(1) FOV(2) 1]));
                % 3D FFT
                tmp_KSP = fftshift(fftn(tmp,FOV));
                for RO = 1:length(currSpokes)
                    kk(:,RO) = squeeze(tmp_KSP(app.ky_samples(1,currSpokes(RO)),:,app.kz_samples(1,currSpokes(RO))));
                end
                kspace(:,currSpokes,coil) = kk;
            end
            clear IMG kk
        else    % nonCartesian sampling using NUFFT and all coils at once
            tmp = Coils.*IMG.*(repmat(permute(SP, [2 3 1]),[FOV(1) FOV(2) 1]));
            
            % build k matrix from all currSpokes
            x = app.kx_samples(:,currSpokes)*app.matrixRL;
            y = app.ky_samples(:,currSpokes)*app.matrixAP;
            z = (app.kz_samples(:,currSpokes)-app.matrixFH/2);
            kloc = cat(3,y,x,z); kloc = permute(kloc, [3 1 2]);
            
            %             % create 3d nufft object and sample
            %             obj = nufft_3d(kloc,[app.matrixAP app.matrixRL app.matrixFH],'radial',1);
            %             kk = obj.fNUFT(tmp);
            %             kspace(:,currSpokes,:) = reshape(kk,[size(x,1),size(x,2),nCh]);
            
            tmp2 = fft(ifftshift(tmp,3),[],3);   % take slices into k-space
            sl_array = unique(z);
            kk = zeros(size(x,1),size(x,2),nCh);
            for sl = 1:length(sl_array)
                currSl = sl_array(sl);
                ind = find(kloc(3,1,:) == currSl);
                kloc2 = kloc(:,:,ind); kloc2(3,:,:) = 0;
                obj = nufft_3d(kloc2,[app.matrixAP app.matrixRL 1],'radial',1,'gpu',0);
                kk(:,ind,:) = reshape(obj.fNUFT(squeeze(tmp2(:,:,currSl+abs(min(sl_array))+1,:))), ...
                    size(x,1), length(ind), size(tmp2,4));
            end
            kspace(:,currSpokes,:) = kk;
            
            clear IMG kk
        end
    end
end
close(h)

% the final step is to add noise as a percentage of max kspace signal
% which means Pnoise/Psignal = 0.1 or Psignal/Pnoise = 10
% then SNRdb = 10*log10(Psignal/Pnoise)
if SNR < 100
    % add noise, IMGs is normalized to 1 so noise = 1/SNR
    NOISE = (1/SNR)*randn(size(kspace)).*exp(1i*2*pi*(rand(size(kspace))));
    kspace = kspace + NOISE;
end