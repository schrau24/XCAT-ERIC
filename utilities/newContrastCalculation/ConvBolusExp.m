function  [y,I] = ConvBolusExp(t,m,k,speedFlag)

if nargin<4
    speedFlag = [];
else
    if strcmp(speedFlag,'slow')
        speedFlag = [];
    elseif strcmp(speedFlag,'fast')
        speedFlag = 'no series';
    end
end


tB = 2*pi/m;

I1 = t>0 & t<tB;
I2 = t>=tB;

y = zeros(size(t));
if nargout==2
    I = zeros(size(t));
    [term1,~,~,I(I1)] = SpecialCosineExp(k*t(I1), m*t(I1));
    y(I1) = t(I1).*term1;
    [term2,~,~,I(I2)]= SpecialCosineExp(k*tB, m*tB);
    I(I2) = -I(I2);
    y(I2) = tB*term2*exp(-k*(t(I2)-tB));
else
    y(I1) = t(I1).*SpecialCosineExp(k*t(I1), m*t(I1),speedFlag);
    y(I2) = tB*SpecialCosineExp(k*tB, m*tB,speedFlag)*exp(-k*(t(I2)-tB));
end