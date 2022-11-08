function  y = ConvBolusGamma(t,m,k)

tB = 2*pi/m;

y = zeros(size(t));
I1 = t>0 & t<tB;
I2 = t>= tB;

ce = SpecialCosineExp(k*tB, m*tB);
cg = SpecialCosineGamma(k*tB, m*tB);

y(I1) = t(I1).^2.*SpecialCosineGamma(k*t(I1), m*t(I1));
y(I2) = tB*((t(I2)-tB)*ce + tB*cg).*exp(-k*(t(I2)-tB));
