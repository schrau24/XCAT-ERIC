% Tissue curve with Extended Kety leakage model using AIF with cosine 
% first-pass and exponential washout impulse response.  If pk.vp not
% present then assumes Kety model with no plasma component
% 
% [ct,ce,cp] = Cosine4AIF_ExtKety(t,aif,pk,param,speedFlag)
%
% t = time array
% aif.ab = bolus amplitude
%    .mb = bolus duration
%    .ae = washout amplitude
%    .me = washout rate constant
%    .t0 = onset time
%    .d  = duration of rectangular injection (optional, default = 0)
%
% tissue parameters - will accept any pair of {Kt,ve,ke}
% pk.Kt  = Ktrans
% pk.ve  = ve
% pk.vp  = vp   (default = 0)
% pk.ke  = kep
% pk.dt  = additional onset time offset
%
% param = 'PMB' or 'MRIW' to indicate parameterisation (default = 'PMB')
%
% speedFlag = 'slow' (default) or 'fast' - indicates how CosineBolusExp is called
%
% outputs are tissue curve and it's components.

function [ct,ce,cp] = Cosine4AIF_ExtKety(t,aif,pk,param,speedFlag)

% set defaults if necessary
if nargin==3, param = 'PMB'; end
if ~isfield(aif,'d'), aif.d = 0; end

if nargin<5
    speedFlag = 'slow';
end

% make sure ke is a field of the structure
if ~isfield(pk,'ke'), pk.ke = pk.Kt/pk.ve; end


% change parameterisation if indicated
if strcmp(param,'MRIW')
    tb = 2*pi/aif.mb;
    aif.ab = aif.ab/tb;
    sc = tb*SpecialCosineExp(tb*aif.me,tb*aif.mb);
    aif.ae = aif.ae/aif.ab/sc;
elseif strcmp(param,'PMB')
    %do nothing
else
    error('unexpected parameterisation indicator!')
end

% calculate curve components with and without convolution with injection
% rectangle
t = t - aif.t0 - pk.dt; % offset time array
if aif.d>0
    cpBolus = aif.ab*ConvBolusRect(t,aif.mb,aif.d)/aif.d;
    cpWashout = aif.ab*aif.ae*ConvBolusExpRect(t,aif.mb,aif.me,aif.d)/aif.d;
    ceBolus = pk.ke*aif.ab*ConvBolusExpRect(t,aif.mb,pk.ke,aif.d)/aif.d;
    ceWashout = pk.ke*aif.ab*aif.ae*ConvBolusExpExpRect(t,aif.mb,aif.me,pk.ke,aif.d)/aif.d;
else
    cpBolus = aif.ab*CosineBolus(t,aif.mb);
    cpWashout = aif.ab*aif.ae*ConvBolusExp(t,aif.mb,aif.me,speedFlag);
    ceBolus = pk.ke*aif.ab*ConvBolusExp(t,aif.mb,pk.ke,speedFlag);
    ceWashout = pk.ke*aif.ab*aif.ae*ConvBolusExpExp(t,aif.mb,aif.me,pk.ke,speedFlag);
end    

% AIF curve
cp = cpBolus + cpWashout;

% EES curve
ce = ceBolus + ceWashout;

% if neither ve or vp are specified then don't calculate ct, just return
if ~isfield(pk,'ve') && ~isfield(pk,'vp')
    ct = zeros(size(t));
    return
end

% set default for vp to give Kety model
if ~isfield(pk,'ve'), pk.ve = pk.Kt/pk.ke; end
if ~isfield(pk,'vp'), pk.vp = 0; end

% tissue curve
ct = zeros(size(t));
ct(t>0) = pk.vp*cp(t>0) + pk.ve*ce(t>0);