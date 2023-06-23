% AIF with cosine first-pass and exponential washout impulse response.
% 
% [cp,cpBolus,cpWashout] = Cosine4AIF(t,aif)
%
% t = time array
% aif.ab = bolus amplitude
%    .mb = bolus duration
%    .ae = washout amplitude
%    .me = washout rate constant
%    .t0 = onset time
%    .d  = duration of rectangular injection (optional, default = 0)
% param = 'PMB' or 'MRIW' to indicate parameterisation (default = 'PMB')
%
% outputs are AIF curve and it's components.

function [cp,cpBolus,cpWashout] = Cosine4AIF(t,aif,param)

% set defaults if necessary
if nargin==2, param = 'PMB'; end
if ~isfield(aif,'d'), aif.d = 0; end

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
if aif.d>0
    cpBolus = aif.ab*ConvBolusRect(t-aif.t0,aif.mb,aif.d)/aif.d;
    cpWashout = aif.ab*aif.ae*ConvBolusExpRect(t-aif.t0,aif.mb,aif.me,aif.d)/aif.d;
else
    cpBolus = aif.ab*CosineBolus(t-aif.t0,aif.mb);
    cpWashout = aif.ab*aif.ae*ConvBolusExp(t-aif.t0,aif.mb,aif.me);
end    

% AIF curve
cp = cpBolus + cpWashout;