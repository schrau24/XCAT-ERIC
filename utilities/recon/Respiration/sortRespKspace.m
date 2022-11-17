function [kspace_sorted, traj_sorted] = sortRespKspace(kspace,RESPPhases,app,isCartesian)
if isCartesian
    labels = ones([app.matrixAP app.matrixFH max(RESPPhases)],'uint8');
    
    kspace_sorted = zeros([app.matrixRL size(kspace,3) prod(size(labels))],'single');
    kspace_averages = zeros([app.matrixRL size(kspace,3) prod(size(labels))],'uint8');
    kspace = single(kspace);
    %     h = waitbar(0,'sorting data...');
    nFE = size(app.ky_samples,2);
    for ro = 1:nFE
        % grab necessaryresp phase
        currPhase   = RESPPhases(ro);
        
        % the ky and kz locations
        currky = app.ky_samples(1,ro); currkz = app.kz_samples(1,ro);
        
        % check to see if this readout line already has data, if yes average,
        % if no simply fill in kspace
        labelInd = sub2ind(size(labels),currky,currkz,currPhase);
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
    
    mask = kspace_sorted ~= 0;
    mask_dims = size(mask);
    m = mask(:,:,:,1,:);
    R = numel(m)/sum(m(:)) /4*pi;
    tmpImg = squeeze(squeeze(mask(round(mask_dims(1)/2),:,:,1,:)));
    fig=figure(200);
    montage(permute(tmpImg,[1 2 4 3]),'Size',[1 max(RESPPhases)])
    title('Masks for respiratory recon')
    clear pic
    fprintf('R by mask: %.2f\n', R);
    
    % not needed for Cartesian
    traj_sorted = [];
else    % sort for nonCartesian
    
    % to make all data size the same with no zeros,
    % find the min number of readouts over all resp frames, and that is a
    % multiple of nSlcs
    nSlcs = length(unique(squeeze(app.kz_samples(1,:))));
    min_lines = numel(find(RESPPhases==1));
    for phase = 2:max(RESPPhases)
        if numel(find(RESPPhases==phase)) < min_lines
            min_lines = numel(find(RESPPhases==phase));
        end
    end
    min_lines = floor(min_lines/nSlcs)*nSlcs;
    
    kspace_sorted = zeros(size(kspace,1),min_lines, size(kspace,3), max(RESPPhases));
    traj_sorted = zeros(size(kspace,1),min_lines, 3, max(RESPPhases));
    
    for phase = 1:max(RESPPhases)        
        % check to see which readouts actually fall into our phase
        C = find(RESPPhases==phase);
        kspace_sorted(:,:,:,phase) = kspace(:,C(1:min_lines),:);
        traj_sorted(:,:,:,phase) = cat(3,app.kx_samples(:,C(1:min_lines)),...
            app.ky_samples(:,C(1:min_lines)), app.kz_samples(:,C(1:min_lines)));
    end
    
end