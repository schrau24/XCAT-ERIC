function kspace_corrected = applyPhaseShifts(kspace, RESPPhases, transforms, app)
% the trajectory
KXX = app.kx_samples;
KXX = 0.5*KXX/max(abs(KXX(:)));    % normalize

KYY = app.ky_samples;
KYY = 0.5*KYY/max(abs(KYY(:)));    % normalize

KZZ = app.kz_samples;
KZZ = 0.5*KZZ/max(abs(KZZ(:)));    % normalize

KSP = kspace;

% correct AP/RL/FH motion, RL offset can be set to 0
RL_offsets = -transforms(2,:)*0;
AP_offsets = -transforms(1,:);
FH_offsets = -transforms(3,:);
for phase = 1:size(transforms,2)
    
    idx = find(RESPPhases==phase+1);
    disp(['resp phase=' num2str(phase+1) ', length=' num2str(numel(idx)) ...
        ', AP=' num2str(round(AP_offsets(phase),1)) ...
        ', RL=' num2str(round(RL_offsets(phase),1)) ...
        ', FH=' num2str(round(FH_offsets(phase),1))])

    % phase in kspace
    KSP(:,idx) = KSP(:,idx) .* exp(2*pi*1i*(...
        RL_offsets(phase) * KXX(:,idx) + ...
        AP_offsets(phase) * KYY(:,idx) + ...
        FH_offsets(phase) * KZZ(:,idx) ));
end
kspace_corrected = KSP;