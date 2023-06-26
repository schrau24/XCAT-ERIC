% adapted from Python code:
% https://github.com/oliverchampion/DCENET
function f = SpecialCosineGamma(x, y)
    x2 = x.^2;
    y2 = y.^2;
    expTerm = (3 + (y2 ./ x2) .* (1 - exp(-x))) - ((y2 + x2) ./ x) .* exp(-x);
    trigTerm = ((x2 - y2) .* (1 - cos(y))) - (2 * x .* y) .* sin(y);
    f = (trigTerm + (y2 .* expTerm)) ./ (y2 + x2).^2;
    f(isnan(f)) = 0;
end