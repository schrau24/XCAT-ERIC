% This function is designed to return an accurate value for a wide range of
% (positive) inputs, which is ensured by using various Taylor series
% approximations.  The third input is a string that can be used to
% validate the function by forcing it to use only one of the approximate
% forms.  The function is designed so that the various Taylor series
% approximations are only coded in one location, which ensures that when
% the validation is used there is a guarantee that the same Taylor series
% code is being used.

function [f,ex,ey,I] = SpecialCosineGamma(x,y,scheme)

% thresholds for taylor series approximations
ex = 1;
ey = 0.8;


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
expTerm(I0y) = (3+y(I0y).^2./x(I0y).^2).*(1-exp(-x(I0y))) - (y(I0y).^2+x(I0y).^2)./x(I0y).*exp(-x(I0y));
yE2 = y(Ix).^2;
xE = x(Ix);
if any(Ix(:))
   expTerm(Ix) =             yE2/2 - 0 ...   % special case 
                    - xE/1.*(yE2/3 - 2 ...
                    - xE/2.*(yE2/4 - 1 ...
                    - xE/3.*(yE2/5     ...
                    - xE/4.*(yE2/6 + 1 ...
                    - xE/5.*(yE2/7 + 2 ...
                    - xE/6.*(yE2/8 + 3 ...
                    - xE/7.*(yE2/9 + 4 ...
                    - xE/8.*(yE2/10 + 5 ...
                    - xE/9.*(yE2/11 + 6 ...
                    - xE/10.*(yE2/12 + 7 ...
                    - xE/11.*(yE2/13 + 8 ...
                    - xE/12.*(yE2/14 + 9 ...
                    - xE/13.*(yE2/15 + 10 ...
                    - xE/14.*(yE2/16 + 11 ...
                    - xE/15.*(yE2/17 + 12 ...
                    - xE/16.*(yE2/18 + 13 ...
                    - xE/17.*(yE2/19 + 14 ...
                    - xE/18.*(yE2/20 + 15 ...
                    - xE/19.*(yE2/21 + 16 ...
                    - xE/20.*(yE2/22 + 17))))))))))))))))))));
end     

% calculate the direct and Taylor series expansions for the cos term
I0x = I0 | Ix;
trigTerm(I0x) = (x(I0x).^2-y(I0x).^2).*(1 - cos(y(I0x))) - 2*x(I0x).*y(I0x).*sin(y(I0x));
yT2 = y(Iy).^2;
xT2 = x(Iy).^2;
xT = 4*x(Iy); % note the extra factor of 4
if any(Iy(:))
trigTerm(Iy) = yT2/(1*2).*  (xT2 - 1*xT + 0      - ... 	% note the last term on this line is an exception to the pattern
			   yT2/(3*4).*  (xT2 - 2*xT + 2*2*3  - ...
			   yT2/(5*6).*  (xT2 - 3*xT + 2*3*5  - ...
			   yT2/(7*8).*  (xT2 - 4*xT + 2*4*7  - ...
			   yT2/(9*10).* (xT2 - 5*xT + 2*5*9  - ...
			   yT2/(11*12).*(xT2 - 6*xT + 2*6*11 - ...
			   yT2/(13*14).*(xT2 - 7*xT + 2*7*13 - ...
			   yT2/(15*16).*(xT2 - 8*xT + 2*8*15 - ...
			   yT2/(17*18).*(xT2 - 9*xT + 2*9*17 - ...
			   yT2/(19*20).*(xT2 - 10*xT + 2*10*19))))))))));
end


% calculate the output function for these cases
I0xy = I0 | Ix | Iy;
f(I0xy) = (trigTerm(I0xy) + y(I0xy).^2.*expTerm(I0xy))./(y(I0xy).^2+x(I0xy).^2).^2;

% calculate the output function for small x and y using two parameter
% Taylor series
yH2 = y(Ixy).^2;
xH = x(Ixy);

if any(Ixy(:))
    f(Ixy) = yH2/(2*3*4).*(1 - xH/5.*(2-xH/6.*(3-xH/7.*(4-xH/8.*(5-xH/9.*(6-xH/10.*(7-xH/11.*(8-xH/12.*(9-xH/13.*(10-xH/14.*(11-xH/15.*(12-xH/16.*(13-xH/17.*(14-xH/18.*(15-xH/19.*(16-xH/20))))))))))))))) - ...
        yH2/(5*6).*  (1 - xH/7.*(2-xH/8.*(3-xH/9.*(4-xH/10.*(5-xH/11.*(6-xH/12.*(7-xH/13.*(8-xH/14.*(9-xH/15.*(10-xH/16.*(11-xH/17.*(12-xH/18.*(13-xH/19.*(14-xH/20.*(15-xH/21.*(16-xH/22))))))))))))))) - ...
        yH2/(7*8).*  (1 - xH/9.*(2-xH/10.*(3-xH/11.*(4-xH/12.*(5-xH/13.*(6-xH/14.*(7-xH/15.*(8-xH/16.*(9-xH/17.*(10-xH/18.*(11-xH/19.*(12-xH/20.*(13-xH/21.*(14-xH/22.*(15-xH/23.*(16-xH/24))))))))))))))) - ...
        yH2/(9*10).* (1 - xH/11.*(2-xH/12.*(3-xH/13.*(4-xH/14.*(5-xH/15.*(6-xH/16.*(7-xH/17.*(8-xH/18.*(9-xH/19.*(10-xH/20.*(11-xH/21.*(12-xH/22.*(13-xH/23.*(14-xH/24.*(15-xH/25.*(16-xH/26))))))))))))))) - ...
        yH2/(11*12).*(1 - xH/13.*(2-xH/14.*(3-xH/15.*(4-xH/16.*(5-xH/17.*(6-xH/18.*(7-xH/19.*(8-xH/20.*(9-xH/21.*(10-xH/22.*(11-xH/23.*(12-xH/24.*(13-xH/25.*(14-xH/26.*(15-xH/27.*(16-xH/28))))))))))))))) - ...
        yH2/(13*14).*(1 - xH/15.*(2-xH/16.*(3-xH/17.*(4-xH/18.*(5-xH/19.*(6-xH/20.*(7-xH/21.*(8-xH/22.*(9-xH/23.*(10-xH/24.*(11-xH/25.*(12-xH/26.*(13-xH/27.*(14-xH/28.*(15-xH/29.*(16-xH/30))))))))))))))) - ...
        yH2/(15*16).*(1 - xH/17.*(2-xH/18.*(3-xH/19.*(4-xH/20.*(5-xH/21.*(6-xH/22.*(7-xH/23.*(8-xH/24.*(9-xH/25.*(10-xH/26.*(11-xH/27.*(12-xH/28.*(13-xH/29.*(14-xH/30.*(15-xH/31.*(16-xH/32)))))))))))))))  ...
        )))))));
end

if nargout==4
    % requested set indicator array to show where the various forms were
    % used to calculate the function
    I = zeros(size(x));
    I(I0) = 1;
    I(Ix) = 2;
    I(Iy) = 3;
    I(Ixy) = 4;
end