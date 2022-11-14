function window = generatewindows(wsize, dim)

% Generate 2 or 3 dimensional hanning window 
%
% Inputs:   
%   wsize:  window size
%   dim:    image dimensions (2D or 3D)
%

if dim == 2
    wsize_x = wsize(2);
    wsize_y = wsize(1);
    wxy = convn(hann(wsize_y), hann(wsize_x)', 'full');
    window = wxy;
elseif dim == 3
    
    wsize_x = wsize(2);
    wsize_y = wsize(1);
    wsize_z = wsize(3);
    
    wxy = convn(hann(wsize_y), hann(wsize_x)', 'full');
    wz = zeros(wsize_y,wsize_x,wsize_z);
    wz((wsize_y+1)/2,(wsize_x+1)/2,:) = hann(wsize_z)';
    window = convn(wxy, wz,'full');
    window = window(wsize_y-(wsize_y-1)/2:wsize_y+(wsize_y-1)/2,wsize_x-(wsize_x-1)/2:wsize_x+(wsize_x-1)/2,:);
else
    error('Error: dim is wrongly selected. dim should be either 2 or 3.');
end
