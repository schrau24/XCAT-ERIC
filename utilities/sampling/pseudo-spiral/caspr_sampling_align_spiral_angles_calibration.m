function caspr_joint = caspr_sampling_align_spiral_angles_calibration(caspr,casprc,ang_map,varargin)
%%Align the spiral angles between the calibration and high res regions

if size(casprc,2) == 0
	caspr_joint = caspr;
    return;
end

% Some dimensions
Sl  = size(caspr,1) + size(casprc,1);
Slc = size(casprc,1);
Ns  = size(caspr,2);
Nsc = size(casprc,2);

% Get angles of both regions
ang_caspr         = ang_map(caspr(2,:));
if size(casprc,1) > 1
    ang_casprc        = ang_map(casprc(2,:));
else
    ang_casprc        = 0;
end

% Reformat calibration region to golden angle
if nargin > 3
    Ns_tot = varargin{1};
else
    if Ns-mod(Ns,Nsc)==0
        Ns_tot = Ns;
    else
        Ns_tot = Ns-mod(Ns,Nsc);
        Ns = Ns_tot;
    end
end

caspr_joint = zeros(Sl,Ns);
for n = 1 : Ns_tot
    idx_calib = 1+mod(n-1,Nsc);
    %[~,pos] = min(abs(ang_casprc(idx_calib) - ang_caspr));
    %ang_caspr(pos) = 1E06;
    caspr_joint(1:Slc,n) = casprc(:,idx_calib);
    caspr_joint(Slc+1:end,n) = caspr(:,1+mod(n-1,size(caspr,2)));
end


% END
end