function [kspaceDCESorted, trajDCESorted] = sortDCEKspace(kspace, app, isCartesian)

% re-sort for DCE, use timings from trajectory calculation to be consistent
nROperFrame = app.nROperFrame;
nFE = size(kspace,2);

% WRITE LABEL
for ii=1:app.nt
    if ii==1
        idx = 1:(nROperFrame);
    elseif ii < app.nt
        idx = (nROperFrame*(ii-1)+1):(nROperFrame*ii);
    else
        idx = (nROperFrame*(ii-1)+1):nFE;
    end
    DCEphases(idx) = ii;
end

if isCartesian
    labels = ones([app.matrixAP app.matrixFH app.nt],'uint8');
    
    kspace_sorted = zeros([app.matrixRL size(kspace,3) prod(size(labels))],'single');
    kspace_averages = zeros([app.matrixRL size(kspace,3) prod(size(labels))],'uint8');
    kspace = single(kspace);
    
    for ro = 1:nFE
        % grab necessary phase
        DCEPhase   = DCEphases(ro);
        
        % the ky and kz locations
        currky = app.ky_samples(1,ro); currkz = app.kz_samples(1,ro);
        
        % check to see if this readout line already has data, if yes average,
        % if no simply fill in kspace
        labelInd = sub2ind(size(labels),currky,currkz,DCEPhase);
        k = kspace_sorted(:,:,labelInd);
        currKLine = squeeze(kspace(:,ro,:));
        kspace_sorted(:,:,labelInd) = kspace_sorted(:,:,labelInd) + currKLine;
        kspace_averages(:,:,labelInd) = kspace_averages(:,:,labelInd) + 1;
        
    end
    
    % Normalize by number of averages
    kspace_sorted = kspace_sorted./single(kspace_averages);
    kspace_sorted(isnan(kspace_sorted)) = complex(0);     % correct for NaN because of division by zero in case of missing k-lines
    kspace_sorted(isinf(kspace_sorted)) = complex(0);
    
    kspace_sorted = reshape(kspace_sorted,[app.matrixRL size(kspace,3) size(labels)]);
    % replace our old kspace with the new one
    kspace_sorted = permute(kspace_sorted,[1 3 4 2 5 6]); clear kspace_averages
    
    % dummy variable for Cartesian
    trajDCESorted = [];
    
    % calculate undersampling
    mask = kspace_sorted ~= 0;
    mask_dims = size(mask);
    m = mask(:,:,:,1,:);
    R = numel(m)/sum(m(:)) /4*pi;
    tmpImg = squeeze(squeeze(mask(round(mask_dims(1)/2),:,:,1,:)));
    figure(200);
    montage(permute(tmpImg,[1 2 4 3]))
    clear pic
    
else
    
    % to make all data size the same with no zeros,
    % find the min number of readouts over all resp frames
    min_lines = numel(find(DCEphases==1));
    for phase = 2:max(DCEphases)
        if numel(find(DCEphases==phase)) < min_lines
            min_lines = numel(find(DCEphases==phase));
        end
    end
    
    kspace_sorted = zeros(size(kspace,1),min_lines, size(kspace,3), max(DCEphases));
    traj_sorted = zeros(size(kspace,1),min_lines, 3, max(DCEphases));
    
    for phase = 1:max(DCEphases)
        kk = zeros(size(kspace,1),min_lines,3);
        
        % check to see which readouts actually fall into our phase
        C = find(DCEphases==phase);
        kspace_sorted(:,:,:,phase) = kspace(:,C(1:min_lines),:);
        traj_sorted(:,:,:,phase) = cat(3,app.kx_samples(:,C(1:min_lines)),...
            app.ky_samples(:,C(1:min_lines)), app.kz_samples(:,C(1:min_lines)));
    end
    trajDCESorted = traj_sorted;
    
    NyquistSampling = pi/2*app.matrixRL;
                    R = NyquistSampling/(size(kspace_sorted,2)/app.matrixFH);
end
fprintf('R by mask: %.2f\n', R);

kspaceDCESorted = kspace_sorted;