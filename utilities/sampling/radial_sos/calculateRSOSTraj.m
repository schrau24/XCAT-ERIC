function [KXX, KYY, KZZ, nArms, nSlcs] = calculateRSOSTraj(kx, nArms, origSlcs)

% other parameters are determined from popup inputdlg
prompt = {'Golden angle (degrees):','aligned (0) or rotating (1):'};
dlgtitle = 'Stack-of-spirals parameters';
dims = [1 50];
definput = {'111.246','0'};
if origSlcs > 3     % only allow half-scan if slices 4 or more
    prompt = cat(2,prompt, {'Half scan factor in FH (0.51 to 1):'});
    definput = cat(2,definput,{'1'});
end
answer = inputdlg(prompt,dlgtitle,dims,definput);

%--Golden angle increment--%
GA = str2double(answer{1});
dangleinc = GA*pi/180;

isRotating = str2double(answer{2});
Zh = 1;
if origSlcs > 3
    Zh  = str2double(answer{3});                    % Half scan in Z
end
nSlcs = ceil(origSlcs*Zh);                          % update slices

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
if nSlcs == 1
    kz = 0;
end
kz = repmat(kz,[1 size(KYY,2)/nSlcs]);
KZZ = repmat(kz,[length(kx) 1]);
