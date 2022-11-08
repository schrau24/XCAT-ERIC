function MR = XCAT_to_MR_DCE(XCAT,FA,TR,TE,Contrast,b0,fatFlag,dyn,enhanced_tissues,C,relaxivity)

% function MR = XCAT_to_MR(XCAT,FA,TR,TE,Contrast, b0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script converts XCAT tissue values to MR contrast based on the SSFP signal equation.
% Christopher W. Roy 2018-12-04
% fetal.xcmr@gmail.com

% T1 and T2 values are stored in a table with the format:
% [T1 @ 1.5 T, T2 @ 1.5 T, T1 @ 3.0 T, T2 @ 3.0 T]
% Note the current implimentation is for 1.5T

% Relaxation values are based off of:
% Stanisz GJ, Odrobina EE, Pun J, Escaravage M, Graham SJ, Bronskill MJ, Henkelman RM. T1, T2 relaxation and magnetization transfer in tissue at 3T. Magnetic resonance in medicine. 2005;54:507�12.
% Portnoy S, Osmond M, Zhu MY, Seed M, Sled JG, Macgowan CK. Relaxation properties of human umbilical cord blood at 1.5 Tesla. Magnetic Resonance in Medicine. 2016;00:1�13.
% https://www.itis.ethz.ch/virtual-population/tissue-properties/XCATbase/relaxation-times/

%Tissue legend:
% 1 Myocardium LV
% 2 myocardium RV
% 3 myocardium la
% 4 myocardium ra
% 5 Blood LV
% 6 Blood RV
% 7 Blood LA
% 8 Blood RA
% 9 body
% 10 muscle
% 11 Brain
% 12 Sinus
% 13 Liver
% 14 gall bladder
% 15 Right Lung
% 16 Left Lung
% 17 esophagus
% 18 esophagus cont
% 19 laryngopharynx
% 20 st wall
% 21 Stomach Contents
% 22 pancreas
% 23 Right kydney cortex
% 24 right kidney medulla
% 25 Left kidney cortex
% 26 left kidney medulla
% 27 adrenal
% 28 Right Renal Pelvis
% 29 Left Renal Pelvis
% 30 spleen
% 31 Ribs
% 32 Cortical Bone
% 33 Spine
% 34 spinal cord
% 35 Bone Marrow
% 36 Artery
% 37 Vein
% 38 bladder
% 39 prostate
% 40 asc lower intestine
% 41 trans lower intestine
% 42 desc lower intestine
% 43 small intestine
% 44 rectum
% 45 seminal vescile
% 46 vas deference
% 47 testicles
% 48 epididymus
% 49 ejac duct
% 50 pericardium
% 51 Cartilage
% 52 Intestine Cavity
% 53 ureter
% 54 urethra
% 55 Lymph
% 56 lymph abnormal
% 57 trach bronch
% 58 Airway
% 59 uterus
% 60 vagina
% 61 right ovary
% 62 left ovary
% 63 FAllopian tubes
% 64 Parathyroid
% 65 Thyroid
% 66 Thymus
% 67 salivary
% 68 Pituitary
% 69 Eye
% 70 eye lens
% 71 lesion
% 72 Fat
% 73 Pancreas tumor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Tissue(1,:)=[1030;40;1471;47];
Tissue(2,:)=[1030;40;1471;47];
Tissue(3,:)=[1030;40;1471;47];
Tissue(4,:)=[1030;40;1471;47];
Tissue(5,:)=[1441;290;1932;275];
Tissue(6,:)=[1441;290;1932;275];
Tissue(7,:)=[1441;290;1932;275];
Tissue(8,:)=[1441;290;1932;275];
Tissue(9,:)=[1008;44;1412;50.00];
Tissue(10,:)=[981.5;36;1232.9;37.20];
Tissue(11,:)=[884;72;1084;69];
Tissue(12,:)=[0;0;0;0];
Tissue(13,:)=[576,46,812,42];
Tissue(14,:)=[576;46;812;42];
Tissue(15,:)=[0;0;0;0];
Tissue(16,:)=[0;0;0;0];
Tissue(17,:)=[576;46;812;42];
Tissue(18,:)=[576;46;812;42];
Tissue(19,:)=[1045.5;37.3;1201;44];
Tissue(20,:)=[981.5;36;1232.9;37.20];
% Tissue(21,:)=[981.5;36;1232.9;37.20];
Tissue(21,:)=[0;0;0;0];
Tissue(22,:)=[584;46;725;43];
Tissue(23,:)=[828;71;1168;66];
Tissue(24,:)=[1412;85;1545;81];
Tissue(25,:)=[828;71;1168;66];
Tissue(26,:)=[1412;85;1545;81];
Tissue(27,:)=[576;46;812;42];
Tissue(28,:)=[200;0.5;302;0.25];
Tissue(29,:)=[200;0.5;302;0.25];
Tissue(30,:)=[1057;79;1328;61];
Tissue(31,:)=[200;0.5;302;0.25];
Tissue(32,:)=[200;0.5;302;0.25];
Tissue(33,:)=[200;0.5;302;0.25];
Tissue(34,:)=[745;74;993;78];
Tissue(35,:)=[549;49;586;49];
Tissue(36,:)=[1585;254;1664;147];
Tissue(37,:)=[1582;181;1584;66];
Tissue(38,:)=[576;46;812;42];
Tissue(39,:)=[1317;88;1597;74];
Tissue(40,:)=[576;46;812;42];
Tissue(41,:)=[576;46;812;42];
Tissue(42,:)=[576;46;812;42];
Tissue(43,:)=[576;46;812;42];
Tissue(44,:)=[576;46;812;42];
Tissue(45,:)=[1317;88;1597;74];
Tissue(46,:)=[576;46;812;42];
Tissue(47,:)=[576;46;812;42];
Tissue(48,:)=[576;46;812;42];
Tissue(49,:)=[576;46;812;42];
Tissue(50,:)=[981.5;36;1232.9;37.20];
Tissue(51,:)=[1024;30;1168;27];
% Tissue(52,:)=[;;;];
Tissue(53,:)=[1434.5;164;1498.3;164];
Tissue(54,:)=[1434.5;164;1498.3;164];
Tissue(55,:)=[5053;468;5053;468];
% Tissue(56,:)=[;;;];
Tissue(57,:)=[1045.5;37.3;1201;44];
Tissue(58,:)=[0;0;0;0];
Tissue(59,:)=[1309;117;1514;79];
Tissue(60,:)=[1135;58;1616;83];
Tissue(61,:)=[576;46;812;42];
Tissue(62,:)=[576;46;812;42];
Tissue(63,:)=[576;46;812;42];
Tissue(64,:)=[1317;88;1597;74];
Tissue(65,:)=[1317;88;1597;74];
Tissue(66,:)=[1317;88;1597;74];
Tissue(67,:)=[1317;88;1597;74];
Tissue(68,:)=[1317;88;1597;74];
Tissue(69,:)=[5053;468;5053;468];
Tissue(70,:)=[1138;26;1138;26];
% Tissue(71,:)=[;;;];
Tissue(72,:)=[312.4;117.4;377;97.5];
Tissue(73,:)=[0;0;1676;64]; % pancreatic tumor, only for 3T

fatTissues = [34,35,72];    % spinal cord, bone marrow, fat


MR=zeros(size(XCAT),'single');

for iTissue=1:length(Tissue)
    
    % if fat 'suppression' on, and tissue is fatty, skip MR contrast
    % calculation (leaves as zero)
    if fatFlag && numel(find(iTissue == fatTissues))==1
    else
        
        if b0 == 1.5
            T1=Tissue(iTissue,1);
            T2=Tissue(iTissue,2);
        else    % b0 = 3.0
            T1 = Tissue(iTissue,3);
            T2 = Tissue(iTissue,4);
        end
        
        if dyn && ismember(iTissue, enhanced_tissues)
            R1 = 1e3/T1;
            R1_new = R1 + relaxivity*C(enhanced_tissues==iTissue);
            T1 = 1e3/R1_new;
        end

        if T1>0 || T2 >0    % if values are 0, then do nothing
            if strcmp(Contrast,'GRE')           % GRE contrast
                %             MR=MR+(XCAT==(iTissue)).*sind(FA)*(1-exp(-TR/T1))*(exp(-TE/T2))/(1-(cosd(FA)*exp(-TR/T1)));
                % no T2(*) component with short TE, set bone to zero
                if isempty(intersect(iTissue,[31, 32, 33]))
                    MR=MR+(XCAT==(iTissue)).*sind(FA)*(1-exp(-TR/T1))/(1-(cosd(FA)*exp(-TR/T1)));
                end
            elseif strcmp(Contrast,'bSSFP')     % bSSFP contrast
                MR=MR+(XCAT==(iTissue)).*sind(FA)*(1-exp(-TR/T1))*exp(-TE/T2)/(1-(exp(-TR/T1)-exp(-TR/T2))*cosd(FA)-exp(-TR/T1)*exp(-TR/T2));
                %MR=MR+(XCAT==(iTissue)).*sind(FA)*(1-exp(-TR/T1))/(1-exp(-TR/T1)*cosd(FA));
            elseif strcmp(Contrast,'SE')        % spin echo, assume single-echo SE, with FA1=90,FA2=180
                MR=MR+(XCAT==(iTissue)).*(1-2*exp(-(TR-TE/2)/T1)+exp(-TR/T1))*exp(-TE/T2);
            else
                errordlg('Unknown MR contrast, use only GRE, bSSFP, or SE. Returning zeros')
            end
        end
    end
    
end



end

