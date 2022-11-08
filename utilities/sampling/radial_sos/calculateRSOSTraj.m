function [KXX, KYY, KZZ, nArms, nSlcs] = calculateRSOSTraj(kx, nArms, origSlcs)

% other parameters are determined from popup inputdlg
prompt = {'Golden angle (degrees):','aligned (0) or rotating (1):',...
    'Half scan factor in FH (0.51 to 1):'};
dlgtitle = 'Stack-of-stars parameters';
dims = [1 50];
definput = {'111.246','0','1'};
answer = inputdlg(prompt,dlgtitle,dims,definput);

%--Golden angle increment--%
GA = str2double(answer{1});
dangleinc = GA*pi/180;

isRotating = str2double(answer{2});
Zh  = str2double(answer{3});                    % Half scan in Z
nSlcs = ceil(origSlcs*Zh);                      % update slices
if mod(nSlcs,2) == 1                            % make even number of slices
    nSlcs = nSlcs+1;
end

% nArms should be a multiple of nSlcs
nArms = nSlcs*ceil(nArms/nSlcs);

%--Loop through spiral arms--%
KXX = zeros(length(kx),nArms);
KYY = zeros(length(kx),nArms);

count = 0;
for arm=1:nArms
    %         aa(arm)=count;
    rot = count*dangleinc;
    KXX(:,arm) = kx*sin(rot);
    KYY(:,arm) = kx*cos(rot);
    if isRotating
        count = count+1;
    else
        if mod(arm,nSlcs) == 0
            count = count+1;
        end
    end
end

% make z locations too, using reversed central ordering
steps = origSlcs-1:-1:1;
steps(2:2:end) = -steps(2:2:end);
loc_order = ones(origSlcs,1);
for s = 2:length(steps)+1
    loc_order(s) = loc_order(s-1) + steps(s-1);
end
kz = flip(loc_order,1)';
% slices to remove for half-scan
kz(kz<(origSlcs+1-nSlcs)) = [];
kz = repmat(kz,[1 size(KYY,2)/nSlcs]);
KZZ = repmat(kz,[length(kx) 1]);
