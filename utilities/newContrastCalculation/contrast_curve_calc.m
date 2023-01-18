function [enhanced_tissue,C,aif] = contrast_curve_calc(contrast_time,timing)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

Hct = 0.4;
aif = struct('ab',2.84/(1-Hct), 'mb', 22.8, 'ae', 1.36, 'me', 0.171, 't0',contrast_time/60, 'd', 0);
[cp,~,~] = Cosine4AIF(timing/60,aif);
%  set enhanced tissue
% enhanced_tissue = [13 22 23 24 25 26 30 73];    
enhanced_tissue = [1 2 3 4 5 6 7 8 13 17 18 20 22 23 24 25 26 30 36 37 40 ...
    41 42 43 50 73];  

ve(1)  = 0.07;              % 1 Myocardium LV
ve(2)  = 0.07;              % 2 myocardium RV
ve(3)  = 0.07;              % 3 myocardium la
ve(4)  = 0.07;              % 4 myocardium ra
ve(5)  = 0.00;              % 5 Blood LV
ve(6)  = 0.00;              % 6 Blood RV
ve(7)  = 0.00;              % 7 Blood la
ve(8)  = 0.00;              % 8 Blood ra
ve(13) = 0.29;              % 13 liver
ve(17) = 0.45;              % 17 esophagus
ve(18) = 0.45;              % 18 esophagus cont
ve(20) = 0.27;              % 20 stomach wall
ve(22) = 0.21;              % 22 Pancreas (from literature)
ve(23) = 0.8;               % 23 right kydney cortex
ve(24) = 0.9;               % 23 right kydney medulla
ve(25) = 0.8;               % 23 left kydney cortex
ve(26) = 0.9;               % 23 left kydney medulla
ve(30) = 0.28;              % 30 spleen
ve(36) = 0.0;               % 36 artery
ve(37) = 0.0;               % 37 vein
ve(40) = 0.19;              % 40 asc lower intestine
ve(41) = 0.19;              % 41 trans lower intestine
ve(42) = 0.19;              % 42 desc lower intestine
ve(43) = 0.19;              % 43 small intestine
ve(50) = 0.07;              % 50 pericardium
ve(73) = 0.37;              % 73 Pancreatic tumor (advanced state, from literature) 

vp(1)  = 0.01;              % 1 Myocardium LV
vp(2)  = 0.01;              % 2 myocardium RV
vp(3)  = 0.01;              % 3 myocardium la
vp(4)  = 0.01;              % 4 myocardium ra
vp(5)  = 1.0 - 0.4;         % 5 Blood LV
vp(6)  = 1.0 - 0.4;         % 6 Blood RV
vp(7)  = 1.0 - 0.4;         % 7 Blood la
vp(8)  = 1.0 - 0.4;         % 8 Blood ra
vp(13) = 0.01;              % 13 liver
vp(17) = 0.01;              % 17 esophagus
vp(18) = 0.01;              % 18 esophagus cont
vp(20) = 0.01;              % 20 stomach wall
vp(22) = 0.04;              % 22 Pancreas (from literature)
vp(23) = 0.2;               % 23 right kydney cortex
vp(24) = 0.1;               % 23 right kydney medulla
vp(25) = 0.2;               % 23 left kydney cortex
vp(26) = 0.1;               % 23 left kydney medulla
vp(30) = 0.06;              % 30 spleen
vp(36) = 1.0 - 0.4;         % 36 artery
vp(37) = 1.0 - 0.4;         % 37 vein
vp(40) = 0.01;              % 40 asc lower intestine
vp(41) = 0.01;              % 41 trans lower intestine
vp(42) = 0.01;              % 42 desc lower intestine
vp(43) = 0.01;              % 43 small intestine
vp(50) = 0.01;              % 50 pericardium
vp(73) = 0.02;              % 73 Pancreatic tumor (advanced state, from literature) 

ke(1)  = 0.46;              % 1 Myocardium LV
ke(2)  = 0.46;              % 2 myocardium RV
ke(3)  = 0.46;              % 3 myocardium la
ke(4)  = 0.46;              % 4 myocardium ra
ke(5)  = 0.0;               % 5 Blood LV
ke(6)  = 0.0;               % 6 Blood RV
ke(7)  = 0.0;               % 7 Blood la
ke(8)  = 0.0;               % 8 Blood rake(13) = 1.09;              
ke(13) = 1.09;              % 13 liver
ke(17) = 1.29;              % 17 esophagus
ke(18) = 1.29;              % 18 esophagus cont
ke(20) = 1.09;              % 20 stomach wall
ke(22) = 1.41;              % 22 Pancreas (from literature)
ke(23) = 0.7;               % 23 right kydney cortex
ke(24) = 0.7;               % 24 right kydney medulla
ke(25) = 0.7;               % 25 left kydney cortex
ke(26) = 0.7;               % 26 left kydney medulla
ke(30) = 1.66;              % 30 spleen
ke(36) = 0.0;               % 36 artery
ke(37) = 0.0;               % 37 vein
ke(40) = 0.84;              % 40 asc lower intestine
ke(41) = 0.84;              % 41 trans lower intestine
ke(42) = 0.84;              % 42 desc lower intestine
ke(43) = 0.84;              % 43 small intestine
ke(50) = 0.46;              % 50 pericardium
ke(73) = 0.42;              % 73 Pancreatic tumor (advanced state, from literature) 

dt(1)  = 6;                 % 1 Myocardium LV
dt(2)  = 6;                 % 2 myocardium RV
dt(3)  = 6;                 % 3 myocardium la
dt(4)  = 6;                 % 4 myocardium ra
dt(5)  = 0.0;               % 5 Blood LV
dt(6)  = 5.0;               % 6 Blood RV
dt(7)  = 0.0;               % 7 Blood la
dt(8)  = 5.0;               % 8 Blood ra
dt(13) = 15.6;              % 13 liver
dt(17) = 8;                 % 17 esophagus
dt(18) = 8;                 % 18 esophagus cont
dt(20) = 10.6;              % 20 stomach wall
dt(22) = 15;                % 22 Pancreas (from literature)
dt(23) = 8;                 % 23 right kydney cortex
dt(24) = 8;                 % 24 right kydney medulla
dt(25) = 8;                 % 25 left kydney cortex
dt(26) = 8;                 % 26 left kydney medulla
dt(30) = 9;                 % 30 spleen
dt(36) = 7.0;               % 36 artery
dt(37) = 15.0;              % 37 vein
dt(40) = 12.7;              % 40 asc lower intestine
dt(41) = 12.7;              % 41 trans lower intestine
dt(42) = 10.7;              % 42 desc lower intestine
dt(43) = 10.7;              % 43 small intestine
dt(50) = 9;                 % 50 pericardium
dt(73) = 15;                % 73 Pancreatic tumor (advanced state, from literature) 

C = zeros([length(enhanced_tissue) length(timing)]);

for i = 1:length(enhanced_tissue)
    t = enhanced_tissue(i);
    pk = struct('ke',ke(t), 'dt',dt(t)/60, 've',ve(t),'vp',vp(t),'t0',contrast_time/60);
    C(i,:) = Cosine4AIF_ExtKety(timing/60, aif, pk);
end

% i = i+1;
% enhanced_tissue = [enhanced_tissue 36];
% C(i,:) = cp;

end