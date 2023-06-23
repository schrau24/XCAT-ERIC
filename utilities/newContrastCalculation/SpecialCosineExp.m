% This function is designed to return an accurate value for a wide range of
% (positive) inputs, which is ensured by using various Taylor series
% approximations.  The third input is a string that can be used to
% validate the function by forcing it to use only one of the approximate
% forms.  The function is designed so that the various Taylor series
% approximations are only coded in one location, which ensures that when
% the validation is used there is a guarantee that the same Taylor series
% code is being used.

function [f,ex,ey,I] = SpecialCosineExp(x,y,scheme)

% thresholds for taylor series approximations
ex = 0.6; %10^(-0.9);
ey = 0.45; %10^(-0.9);

% if scheme is not input then calculate function using 4 sub-forms
if nargin==2 || isempty(scheme)
    
    % indicators for the 4 sub-forms of the function
    I0  = x>ex  & y>ey;
    Ix  = x<=ex & y>ey;
    Iy  = x>ex  & y<=ey;
    Ixy = x<=ex & y<=ey;
    
    % if scheme is input then set indicators to make it return the sub-form
    % requested in scheme
elseif nargin==3 && ~isempty(scheme)
    switch scheme
        case 'no series'
            I0  = true(size(x));
            Ix  = false(size(y));
            Iy  = false(size(x));
            Ixy = false(size(x));
        case 'x series'
            I0  = false(size(x));
            Ix  = true(size(y));
            Iy  = false(size(x));
            Ixy = false(size(x));
        case 'y series'
            I0  = false(size(x));
            Ix  = false(size(y));
            Iy  = true(size(x));
            Ixy = false(size(x));
        case 'xy series'
            I0  = false(size(x));
            Ix  = false(size(y));
            Iy  = false(size(x));
            Ixy = true(size(x));
        otherwise
            error('Unknown computation scheme')
    end
end

[expTerm,trigTerm,f] = deal(zeros(size(x)));

% calculate the direct and Taylor series expansions for the exp term
I0y = I0 | Iy;
%Ixxy = Ix | Ixy;
expTerm(I0y) = (1 - exp(-x(I0y)))./x(I0y);
xH = x(Ix);
if any(Ix(:))
    expTerm(Ix) =   1 - xH/2.*( ...
        1 - xH/3.*( ...
        1 - xH/4.*( ...
        1 - xH/5.*( ...
        1 - xH/6.*( ...
        1 - xH/7.*( ...
        1 - xH/8.*( ...
        1 - xH/9.*( ...
        1 - xH/10.*( ...
        1 - xH/11.*( ...
        1 - xH/12.*( ...
        1 - xH/13.*( ...
        1 - xH/14.*( ...
        1 - xH/15)))))))))))));
end

% calculate the direct and Taylor series expansions for the cos term
I0x = I0 | Ix;
%Iyxy = Iy | Ixy;
trigTerm(I0x) = x(I0x).*(1 - cos(y(I0x))) - y(I0x).*sin(y(I0x));
yH2 = y(Iy).^2; 
xH = x(Iy);
if any(Iy(:))
    trigTerm(Iy) =   yH2/(1*2).*(xH - 2 - ...
        yH2/(3*4).*(xH - 4 - ...
        yH2/(5*6).*(xH - 6 - ...
        yH2/(7*8).*(xH - 8 - ...
        yH2/(9*10).*(xH - 10 - ...
        yH2/(11*12).*(xH - 12 - ...
        yH2/(13*14).*(xH - 14)))))));
end
I0xy = I0 | Ix | Iy;
f(I0xy) = (trigTerm(I0xy) + y(I0xy).^2.*expTerm(I0xy))./(x(I0xy).^2+y(I0xy).^2);

xH = x(Ixy);
yH2 = y(Ixy).^2;
if any(Ixy(:))
    f(Ixy) = yH2/(2*3).*(1 - xH/4.*(1-xH/5.*(1-xH/6.*(1-xH/7.*(1-xH/8.*(1-xH/9.*(1-xH/10.*(1-xH/11.*(1-xH/12.*(1-xH/13.*(1-xH/14.*(1-xH/15))))))))))) - ...
        yH2/(4*5).*(1 - xH/6.*(1-xH/7.*(1-xH/8.*(1-xH/9.*(1-xH/10.*(1-xH/11.*(1-xH/12.*(1-xH/13.*(1-xH/14.*(1-xH/15.*(1-xH/16.*(1-xH/17))))))))))) - ...
        yH2/(6*7).*(1 - xH/8.*(1-xH/9.*(1-xH/10.*(1-xH/11.*(1-xH/12.*(1-xH/13.*(1-xH/14.*(1-xH/15.*(1-xH/16.*(1-xH/17.*(1-xH/18.*(1-xH/19))))))))))) - ...
        yH2/(8*9).*(1 - xH/10.*(1-xH/11.*(1-xH/12.*(1-xH/13.*(1-xH/14.*(1-xH/15.*(1-xH/16.*(1-xH/17.*(1-xH/18.*(1-xH/19.*(1-xH/20.*(1-xH/21))))))))))) - ...
        yH2/(10*11).*(1 - xH/12.*(1-xH/13.*(1-xH/14.*(1-xH/15.*(1-xH/16.*(1-xH/17.*(1-xH/18.*(1-xH/19.*(1-xH/20.*(1-xH/21.*(1-xH/22.*(1-xH/23))))))))))) - ...
        yH2/(12*13).*(1 - xH/14.*(1-xH/15.*(1-xH/16.*(1-xH/17.*(1-xH/18.*(1-xH/19.*(1-xH/20.*(1-xH/21.*(1-xH/22.*(1-xH/23.*(1-xH/24.*(1-xH/25)))))))))))  ...
        ))))));
end

if nargout==4
    % requested set indicator array to show where the various forms were
    % used to calculate the function
    I = zeros(size(x));
    I(I0) =  1;
    I(Ix) =  2;
    I(Iy) =  3;
    I(Ixy) = 4;
end