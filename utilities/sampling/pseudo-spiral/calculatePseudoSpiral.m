function [ky, kz, Sl] = calculatePseudoSpiral(matrix, length_scan, TR, DCETimeResolution)
% Calculate the pseudo-spiral sampling pattern using a dual-density caspr
% approach. 
% Inputs:
% matrix:               KY X KZ sampled Cartesian matrix
% length_scan:          total scan time, in seconds
% TR:                   in ms
% DCETimeResolution:    reconstructed frame length, in seconds
% Outputs:
% ky, kz:               sampled coordinate lists
% Sl:                   the length of one spiral interleaf

% author: Eric Schrauben, e.m.schrauben@amsterdamumc.nl

Y   = matrix(1);                            % Y Matrix size
Z   = matrix(2);                            % Z Matrix size

% other parameters are determined from popup inputdlg
prompt = {sprintf('Spiral interleaf length (%i to %i):',round(0.2*Y),round(0.5*Y)),...
    sprintf('Points per interleaf for central density (%i to %i):', round(0.1*Y), round(0.25*Y)),...
    'Half scan factor in AP (0.51 to 1):', 'Half scan factor in FH (0.51 to 1):'};
dlgtitle = 'Pseudo-spiral parameters';
dims = [1 50];
definput = {sprintf('%i',round(0.33*Y)),sprintf('%i',round(0.16667*Y)),'0.75','0.75'};
answer = inputdlg(prompt,dlgtitle,dims,definput);

Sl = str2double(answer{1});                 % Spiral interleaf length, hard-coded
Sl_time = TR/1000 * Sl;                     % Spiral interleaf time, in seconds

% calculate number of spiral interleaves based on scanTime
Ns = ceil(length_scan / Sl_time);

% Dual density region
% Number of points on the interleaf for central density part
dd_Sl = str2double(answer{2});             
% Number of shots required to fill the central density part
% this can also be used to determine a 'frame' length for DCE
% (e.g. how many shots fit in one frame)
dd_Ns = ceil(DCETimeResolution / (TR/1000 * Sl));
nFrames = floor(Ns/dd_Ns);

Ys    = 1;                  % Sense in Y, for now don't play with these
Zs    = 1;                  % Sense in Z
Yh  = str2double(answer{3});                  % Half scan in Y
Zh  = str2double(answer{4});                  % Half scan in Z

% accelerated Dual density CASPR
% this outputs the samples on a grid for all arms (Sl X Ns)
samp_scheme_ddcaspr = generate_ddcaspr_samp_scheme(Y,Z,Yh,Zh,Ys,Zs,Sl,dd_Sl,Ns,dd_Ns);

% loop over arms and create a separate sampling grid for each DCE time
% frame (Y X Z X nFrames)
Cart_grid = zeros(Y,Z,nFrames);
for n = 1 : nFrames
    % grab arms
    if n < nFrames
        idx = (n - 1)*dd_Ns + 1 : n*dd_Ns;
    else
        idx = (n-1)*dd_Ns+1 : Ns;
    end
    
    % convert to mask
    mask = zeros(Y,Z);
    tmp = samp_scheme_ddcaspr(:,idx);
    mask(tmp(:)) = 1;
    Cart_grid(:,:, n) = mask;
end

% convert to PROUD sampling list
[ky, kz] = ind2sub([Y Z], samp_scheme_ddcaspr(:));