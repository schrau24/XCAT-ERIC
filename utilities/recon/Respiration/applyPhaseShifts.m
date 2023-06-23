function kspace_corrected = applyPhaseShifts(kspace, RESPPhases, transforms, app)
% the trajectory
KYY = app.ky_samples;
KYY = 0.5*KYY/max(abs(KYY(:)));    % normalize

KZZ = app.kz_samples;
KZZ = 0.5*KZZ/max(abs(KZZ(:)));    % normalize

KSP = kspace;

% correct AP/RL/FH motion,
AP_offsets = -transforms(2,:);
FH_offsets = -transforms(1,:);
for phase = 1:size(transforms,2)
    
    idx = find(RESPPhases==phase+1);
    disp(['resp phase=' num2str(phase+1) ', length=' num2str(numel(idx)) ...
        ', AP=' num2str(round(AP_offsets(phase),1)) ...
        ', FH=' num2str(round(FH_offsets(phase),1))])

    % phase in kspace
    KSP(:,idx) = KSP(:,idx) .* exp(2*pi*1i*(...
        AP_offsets(phase) * KYY(:,idx) + ...
        FH_offsets(phase) * KZZ(:,idx) ));
end
kspace_corrected = KSP;