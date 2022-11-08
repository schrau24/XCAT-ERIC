function [KXX, KYY, KZZ, interleaves, nArms, nSlcs] = calculateSpiralTraj(app, nArms, origSlcs)

% other parameters are determined from popup inputdlg
prompt = {'Golden angle (degrees):','aligned (0) or rotating (1):',...
    'Half scan factor in FH (0.51 to 1):'};
dlgtitle = 'Stack-of-spirals parameters';
dims = [1 50];
definput = {'222.492','0','1'};
answer = inputdlg(prompt,dlgtitle,dims,definput);

%%%%% add calculation of one spiral here %%%%%
FOV = 374/10;                       % hard-coded, in mm
baseRes = app.matrixRL;
BW = 456.997;                       % hard-coded from Philips 3T Achieva gradient mode 'maximum'
dt = 0.5/(BW*baseRes);              % in s
Smax                = 16.670*1000;  % hard-coded from Philips 3T Achieva gradient mode 'maximum'
Gmax                = 1.0;          % hard-coded from Philips 3T Achieva gradient mode 'maximum'
krmax = 0.5/FOV*baseRes;

% iterative loop to determine number of interleaves based on TR
% here we assume spiral duration is limited to ~0.5*TR
targetSpiralDur = 0.5*app.TR_sim.Value/1000;     % in s

% initial guess using interleaves=8
interleaves = 5;
[k0,~] = vds(Smax,Gmax,dt,interleaves,FOV,krmax);
calcSpiralDur = length(k0)*dt;
it = 0;
clc
disp('Calculating best spiral interleave number')
fprintf('iteration=%i, interleaves=%i, delta duration=%0.1f ms\n', it, interleaves,...
    1000*(targetSpiralDur-calcSpiralDur));
while abs(targetSpiralDur-calcSpiralDur) > 0.0001
    it = it+1;
    interleaves = interleaves+1;
    [k0,~] = vds(Smax,Gmax,dt,interleaves,FOV,krmax);
    calcSpiralDur = length(k0)*dt;
    fprintf('iteration=%i, interleaves=%i, delta duration=%0.1f ms\n', it, interleaves,...
        1000*(targetSpiralDur-calcSpiralDur));
    
    % if sign of difference suddenly switches (>0), then we can just accept
    % the last value
    if sign(targetSpiralDur-calcSpiralDur)==1
        break
    end
end

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
KXX = zeros(length(k0),nArms);
KYY = zeros(length(k0),nArms);

count = 0;
for arm=1:nArms
    rot = count*dangleinc;
    KXX(:,arm) = (real(k0)*sin(rot)+imag(k0)*cos(rot));
    KYY(:,arm) = (-real(k0)*cos(rot)+imag(k0)*sin(rot));
    if isRotating
        count = count+1;
    else
        if mod(arm,nSlcs) == 0
            count = count+1;
        end
    end
end

% normalize KXX and KYY
KXX = 0.5*KXX/max(max([KXX(:) KYY(:)]));
KYY = 0.5*KYY/max(max([KXX(:) KYY(:)]));

% make z locations too, using reversed central ordering
steps = origSlcs-1:-1:1;
steps(2:2:end) = -steps(2:2:end);
loc_order = ones(nSlcs,1);
for s = 2:length(steps)+1
    loc_order(s) = loc_order(s-1) + steps(s-1);
end
kz = flip(loc_order,1)';
% slices to remove for half-scan
kz(kz<(origSlcs+1-nSlcs)) = [];
kz = repmat(kz,[1 size(KYY,2)/nSlcs]);

KZZ = repmat(kz,[length(k0) 1]);
