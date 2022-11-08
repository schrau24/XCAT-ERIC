function rings = caspr_sampling_concentric_rings(rad_map,Ns,Sl,Nsc,Slc)
% Divide the phase encodes into concentric rings

if nargin < 5
    Slc = 0;
    Nsc = 0;    
else
    Sl = Sl - Slc;
end

% Case Sl == 0, only sample k0
if Sl == 0
    rings = zeros(size(rad_map)); 
    rings(round(0.5+size(rad_map,1)/2),round(0.5+size(rad_map,2)/2)) = 1;
    return;
end

% Sort radii
[rad_sorted,rad_idx]                 = sort(rad_map(:));
rad_idx(find(rad_sorted >= 1E03))    = [];
rad_sorted(find(rad_sorted >= 1E03)) = [];

% Calibration or non-calibration zone
if Slc > 0 && Nsc > 0 % non-calibration
    rad_idx(1:1+Nsc*(Slc-1))               = [];
    rad_sorted(1:1+Nsc*(Slc-1))            = [];
else   % calibration
    rad_idx(2+Ns*Sl:end)             = [];
    rad_sorted(2+Ns*Sl:end)          = [];
    rad_idx(1)                       = [];
    rad_sorted(1)                    = [];
end

rings = zeros(size(rad_map)); 
Npe   = numel(rad_idx);
dPE   = ceil(Npe / Sl);
for r = 1 : Sl
    idx = round(1 + (r-1) * dPE: r * dPE);
    idx(idx>numel(rad_idx)) = [];
    rings(rad_idx(idx)) = r;
end

% END
end