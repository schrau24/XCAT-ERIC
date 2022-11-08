function [y,nI] = ConvBolusExpExp(t,m,k1,k2,speedFlag)

if nargin<5
    speedFlag = 'slow';
end

tol = 1e-4;
% kT = 2*pi*abs(k2-k1)/tol;
% if (m<kT) && (kT<0.5*(k1+k2))
%     tT = 2*pi/m;
% else
    tT = tol/abs(k2-k1);
% end

Ig = t>0 & t<tT;
Ie = t>=tT;

y = zeros(size(t));

y(Ig) = ConvBolusGamma(t(Ig),m,0.5*(k1+k2));
y1 = ConvBolusExp(t(Ie),m,k1,speedFlag);
y2 = ConvBolusExp(t(Ie),m,k2,speedFlag);
y(Ie) = (y1-y2)/(k2-k1);

nI = [nnz(Ie) nnz(Ig)];
