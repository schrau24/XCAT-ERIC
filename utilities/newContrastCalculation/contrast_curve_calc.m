function [enhanced_tissue,C,aif] = contrast_curve_calc(contrast_time,timing)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

Hct = 0.4;
aif = struct('ab',2.84/(1-Hct), 'mb', 22.8, 'ae', 1.36, 'me', 0.171, 't0',contrast_time/60, 'd', 0);
[cp,~,~] = Cosine4AIF(timing/60,aif);
% enhanced_tissue = [13 23 24 25 26 30 22 73];
enhanced_tissue = [13 30 22 73];

ve(13) = 0.29; % 13 liver
ve(23) = 0.8; % 23 right kydney cortex
ve(24) = 0.9; % 23 right kydney medulla
ve(25) = 0.8; % 23 left kydney cortex
ve(26) = 0.9; % 23 left kydney medulla
ve(30) = 0.28; % 30 spleen
ve(22) = 0.21; % 22 Pancreas (from literature)
ve(73) = 0.37; % 73 Pancreatic tumor (advanced state, from literature) 

vp(13) = 0.01; % 13 liver
vp(23) = 0.2; % 23 right kydney cortex
vp(24) = 0.1; % 23 right kydney medulla
vp(25) = 0.2; % 23 left kydney cortex
vp(26) = 0.1; % 23 left kydney medulla
vp(30) = 0.06; % 30 spleen
vp(22) = 0.04; % 22 Pancreas (from literature)
vp(73) = 0.02; % 73 Pancreatic tumor (advanced state, from literature) 

ke(13) = 1.09; % 13 liver
ke(23) = 0.7; % 23 right kydney cortex
ke(24) = 0.7; % 23 right kydney medulla
ke(25) = 0.7; % 23 left kydney cortex
ke(26) = 0.7; % 23 left kydney medulla
ke(30) = 1.66; % 30 spleen
ke(22) = 1.41; % 22 Pancreas (from literature)
ke(73) = 0.42; % 73 Pancreatic tumor (advanced state, from literature) 

dt(13) = 15; % 13 liver
dt(23) = 15; % 23 right kydney cortex
dt(24) = 15; % 23 right kydney medulla
dt(25) = 15; % 23 left kydney cortex
dt(26) = 15; % 23 left kydney medulla
dt(30) = 15; % 30 spleen
dt(22) = 15; % 22 Pancreas (from literature)
dt(73) = 15; % 73 Pancreatic tumor (advanced state, from literature) 

C = zeros([length(enhanced_tissue) length(timing)]);

for i = 1:length(enhanced_tissue)
    t = enhanced_tissue(i);
    pk = struct('ke',ke(t), 'dt',dt(t)/60, 've',ve(t),'vp',vp(t),'t0',contrast_time/60);
    C(i,:) = Cosine4AIF_ExtKety(timing/60, aif, pk);
end

i = i+1;
enhanced_tissue = [enhanced_tissue 36];
C(i,:) = cp;

end