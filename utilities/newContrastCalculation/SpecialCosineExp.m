% adapted from Python code:
% https://github.com/oliverchampion/DCENET
function f = SpecialCosineExp(x, y)
    expTerm = (1 - exp(-x)) ./ x;

    trigTerm = (x .* (1 - cos(y))) - (y .* sin(y));
    f = (trigTerm + (y.^2 .* expTerm)) ./ (x.^2 + y.^2);
    f(isnan(f)) = 0;
end