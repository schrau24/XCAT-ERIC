function [samp_scheme,Ns_calib] = generate_ddcaspr_samp_scheme(Y,Z,Yh,Zh,Ys,Zs,Sl,dd_Sl,Ns,dd_Ns)
%Generate the dual-density cartesian spiral sampling scheme 
% Y      = matrix size Y
% Z      = matrix size Z
% Yh     = half scan factor Y [0.6-1]
% Zh     = half scan factor Z [0.6-1]
% Ys     = sense factor Y, must be an integer
% Zs     = sense factor Z, must be an integer
% Sl     = Shot length, i.e. number of phase encodes on a spiral interleaf
% dd_Sl  = Number of phase encodes per shot to use for central density
% Ns     = Total number of spiral interleaves, this dictates scan time
% dd_Nsc = Number of shots to fully sample central density

%% Check input
samp_scheme = [];
if mod(Y,1) > 0 || mod(Z,1) > 0 || mod(Ys,1) > 0 || mod(Zs,1) > 0 || mod(Ns,1) > 0 || mod(Sl,1) > 0 || mod(dd_Ns,1) > 0 || mod(dd_Sl,1) > 0
    display('Wrong input, most should be integers values except half scan.')
    return;
end

%% Run algorithm
% Create radius map for angle and magnitude
[rad_map,ang_map]     = caspr_sampling_radius_angle_map(Y,Z,Yh,Zh,Ys,Zs);
if (1+(dd_Sl-1)*dd_Ns) < numel(find(rad_map<1E03)) / 2
    Ns_max = ceil((numel(find(rad_map<1E03))-(1+(dd_Sl-1)*dd_Ns)) / (Sl-dd_Sl)); % Calculate maximum number of spirals for one time frame
else
    Ns_max = ceil((1+(dd_Sl-1)*dd_Ns) / (Sl-dd_Sl));
end

% Divide k-space into concentric rings
rings             = caspr_sampling_concentric_rings(rad_map,Ns_max,Sl,dd_Ns,dd_Sl);
ringsc            = caspr_sampling_concentric_rings(rad_map,dd_Ns,dd_Sl-1);

% Divide rings into angular indices
angular_segments  = caspr_sampling_angular_segments(ang_map,rings);
angular_segmentsc = caspr_sampling_angular_segments(ang_map,ringsc);

% Create caspr based on rings and angular segments
caspr             = caspr_sampling_create_spirals(rings,angular_segments);
dd_caspr          = caspr_sampling_create_spirals(ringsc,angular_segmentsc);
dd_caspr          = cat(1,repmat(sub2ind([Y Z],round(0.5+Y/2),round(0.5+Z/2)),[1 dd_Ns]),dd_caspr);

% Shuffle caspr to golden angle 
caspr             = caspr_sampling_shuffle_ga(caspr,ang_map);
dd_caspr          = caspr_sampling_shuffle_ga(dd_caspr,ang_map);

% Match calibration and non calibration with spiral angles and add k0 line 
caspr             = caspr_sampling_align_spiral_angles_calibration(caspr,dd_caspr,ang_map);
while size(caspr,1) < Sl
    caspr(end+1,:) = caspr(end,:);
end

% If sense is enabled first scan full calib region
Ns_calib = 0;
if Ys > 1 || Zs > 1
    Ns_calib  = ceil(dd_Sl * dd_Ns / Sl);
    calib_map = caspr_sampling_radius_angle_map(Y,Z,Yh,Zh,1,1);
    [~,idx]   = sort(calib_map(:));
    idx_calib = idx(1:Ns_calib*Sl);
    caspr     = cat(2,reshape(idx_calib,[Sl,Ns_calib]),caspr); 
end

% Add till all shots are filled
while size(caspr,2) < Ns
    caspr = cat(2,caspr,caspr(:,Ns_calib+1:end));
end
samp_scheme = caspr(:,1:Ns);
% END
end