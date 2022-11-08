function [ky, kz, DCETimeResolution] = calculateCartesian(matrix, length_scan, TR, DCETimeResolution)
% Calculate the Cartesian sampling pattern
% approach.
% Inputs:
% matrix:               KY X KZ sampled Cartesian matrix
% length_scan:          total scan time, in seconds
% TR:                   in ms
% DCETimeResolution:    reconstructed frame length, in seconds
% Outputs:
% ky, kz:               sampled coordinate lists

% author: Eric Schrauben, e.m.schrauben@amsterdamumc.nl

Y   = matrix(1);                            % Y Matrix size
Z   = matrix(2);                            % Z Matrix size

% other parameters are determined from popup inputdlg
prompt = {'Elliptical scanning? (0:no, 1:yes):',...
    'Half scan factor in AP (0.51 to 1):', 'Half scan factor in FH (0.51 to 1):'};
dlgtitle = 'Cartesian parameters';
dims = [1 50];
definput = {'1','0.75','0.75'};
answer = inputdlg(prompt,dlgtitle,dims,definput);

isElliptical = str2double(answer{1});
Yh  = str2double(answer{2});                  % Half scan in Y
Zh  = str2double(answer{3});                  % Half scan in Z

% Define the positions of all phase encodes
cp = [round(.5+Y/2) round(.5+Z/2)];
kpos = 1E3 * ones(Y,Z);
for y = 1 : Y * Yh
    for z = 1 : Z * Zh
        kpos(y,z) = (y - cp(1)) / Y + i * (z - cp(2)) / Z;
    end
end

if isElliptical
    % Delete square corners from kpos
    for y = 1 : Y
        for z = 1 : Z
            if real(Y*kpos(y,z))^2/(Y/2)^2 + imag(Z*kpos(y,z))^2/(Z/2)^2 >= 1
                kpos(y,z) = 1E3;
            end
        end
    end
end

kpos = abs(kpos);
kpos(kpos<1000) = 1; kpos(kpos==1000) = 0;
kpos_vec = find(kpos(:));
% re-order kpos_vec to do centric sampling
cp_ind = sub2ind(size(kpos),cp(1),cp(2));
[~,ind] = sort(abs(kpos_vec-cp_ind));
kpos_vec = kpos_vec(ind);

% Calculate time to fill one k-space, if this is larger than
% DCETimeResolution, need to update app
timeForOneKspace = sum(kpos(:))*(TR/1000);
if timeForOneKspace > DCETimeResolution
    DCETimeResolution = timeForOneKspace;
    addedSamples = [];
else    % if less, then we can resample the central region?

    leftOverTRs = floor((DCETimeResolution - timeForOneKspace)/(TR/1000));
    
    % find the leftOverTRs of sampled values that are closest to cp
    samplesForward = cp_ind+1:1:cp_ind+2*leftOverTRs; samplesForward(samplesForward>(Y*Z)) = [];
    samplesBackward = cp_ind-1:-1:cp_ind-2*leftOverTRs; samplesBackward(samplesBackward<0) = [];
    ii = 1;
    addedSamples = [];
    while length(addedSamples) < leftOverTRs
        if kpos(samplesForward(ii))==1
            addedSamples = cat(1,addedSamples,samplesForward(ii));
            
            if length(addedSamples) == leftOverTRs
                break;
            end
        end
        
        if kpos(samplesBackward(ii))==1
            addedSamples = cat(1,addedSamples,samplesBackward(ii));
        end
        ii = ii+1;
    end
end

nFrames = floor(length_scan/DCETimeResolution);

outputKpos = repmat(cat(1,kpos_vec,addedSamples),[1 nFrames]);

% convert to PROUD sampling list
[ky, kz] = ind2sub([Y Z], outputKpos(:));