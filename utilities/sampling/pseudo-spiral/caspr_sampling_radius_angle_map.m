function [rad_map,ang_map] = caspr_sampling_radius_angle_map(X,Y,Xhalf,Yhalf,Xsense,Ysense)

% Check input for sense
if nargin < 6
    Xsense = 1;
    Ysense = 1;
end

% Check input for half scan
if nargin < 4
    Xhalf = 1;
    Yhalf = 1;
end
       
% Define the positions of all phase encodes
cp = [round(.5+X/2) round(.5+Y/2)];
kpos = 1E3 * ones(X,Y);
for x = 1 : X * Xhalf
	for y = 1 : Y * Yhalf
		kpos(x,y) = (x - cp(1)) / X + i * (y - cp(2)) / Y;
	end
end

% Delete square corners from kpos
for x = 1 : X
    for y = 1 : Y
        if real(X*kpos(x,y))^2/(X/2)^2 + imag(Y*kpos(x,y))^2/(Y/2)^2 >= 1
            kpos(x,y) = 1E3;
        end
    end
end

rad_map = abs(kpos);
ang_map = angle(kpos);

% Delete sense parts from entries
tmp_x = round([flip(round(X/2+0.5):-Xsense:1) round(X/2+0.5) + Xsense:Xsense:X]);
tmp_y = round([flip(round(Y/2+0.5):-Ysense:1) round(Y/2+0.5) + Ysense:Ysense:Y]);
mask_x = zeros(X,1);
mask_x(tmp_x) = 1;
mask_y = zeros(Y,1);
mask_y(tmp_y) = 1;
mask = mask_x*mask_y';

% Apply caipi
if Xsense > 1 && Ysense > 1
    if Xsense < Ysense
        %for n = 1 : 2*Ysense : Y
        for n = round([cp(2)-Ysense:-2*Ysense:1 cp(2)+Ysense:2*Ysense:Y])          
            shift = round((1./(Xsense-1)+1) ./ 2);
            mask(:,n) = mask([shift+1:end 1:shift],n);
            n
        end
    else
        %for n = 1 : 2*Xsense : X
        for n = round([cp(1)-Xsense:-2*Xsense:1 cp(1)+Xsense:2*Xsense:X])
            shift = round((1./(Ysense-1)+1) ./ 2);
            mask(n,:) = mask(n,[shift+1:end 1:shift]);
        end
    end
end

idx = find(mask==0);
rad_map(idx) = 1E3;
ang_map(idx) = 0;

% END
end