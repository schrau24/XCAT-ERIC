function y = CosineBolus(t,m)

z = m*t;
I1 = z<0.2 & z>0;
I2 = z>=0.2 & z<(2*pi);

y = zeros(size(t));

if any(I1)
    z2 = z(I1).^2;
    y(I1) = z2/(1*2).*(1-z2/(3*4).*(1-z2/(5*6).*(1-z2/(7*8).*(1-z2/(9*10)))));
end

y(I2) = 1 - cos(z(I2));
