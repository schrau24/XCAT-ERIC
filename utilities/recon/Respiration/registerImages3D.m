function [T, RegisteredImage] = registerImages3D(MOVING,FIXED)
%registerImages  Register grayscale images using auto-generated code from Registration Estimator app.
%  [MOVINGREG] = registerImages(MOVING,FIXED) Register grayscale images
%  MOVING and FIXED using auto-generated code from the Registration
%  Estimator app. The values for all registration parameters were set
%  interactively in the app and result in the registered image stored in the
%  structure array MOVINGREG.

% Auto-generated by registrationEstimator app on 28-Nov-2018
%-----------------------------------------------------------

% Default spatial referencing objects
xLimits = [-size(FIXED,1)/2+0.5 size(FIXED,1)/2-0.5];
yLimits = [-size(FIXED,2)/2+0.5 size(FIXED,2)/2-0.5];
zLimits = [-size(FIXED,3)/2+0.5 size(FIXED,3)/2-0.5];
fixedRefObj = imref3d(size(FIXED),xLimits,yLimits,zLimits);
movingRefObj = imref3d(size(FIXED),xLimits,yLimits,zLimits);

% Intensity-based registration
[optimizer, metric] = imregconfig('monomodal');
optimizer.GradientMagnitudeTolerance = 1.00000e-04;
optimizer.MinimumStepLength = 1.00000e-05;
optimizer.MaximumStepLength = 6.25000e-02;
optimizer.MaximumIterations = 100;
optimizer.RelaxationFactor = 0.500000;

% Align centers
[xFixed,yFixed,zFixed] = meshgrid(1:size(FIXED,2),1:size(FIXED,1),1:size(FIXED,3));
[xMoving,yMoving,zMoving] = meshgrid(1:size(MOVING,2),1:size(MOVING,1),1:size(MOVING,3));
sumFixedIntensity = sum(FIXED(:));
sumMovingIntensity = sum(MOVING(:));
fixedXCOM = (fixedRefObj.PixelExtentInWorldX .* (sum(xFixed(:).*double(FIXED(:))) ./ sumFixedIntensity)) + fixedRefObj.XWorldLimits(1);
fixedYCOM = (fixedRefObj.PixelExtentInWorldY .* (sum(yFixed(:).*double(FIXED(:))) ./ sumFixedIntensity)) + fixedRefObj.YWorldLimits(1);
fixedZCOM = (fixedRefObj.PixelExtentInWorldZ .* (sum(zFixed(:).*double(FIXED(:))) ./ sumFixedIntensity)) + fixedRefObj.ZWorldLimits(1);
movingXCOM = (movingRefObj.PixelExtentInWorldX .* (sum(xMoving(:).*double(MOVING(:))) ./ sumMovingIntensity)) + movingRefObj.XWorldLimits(1);
movingYCOM = (movingRefObj.PixelExtentInWorldY .* (sum(yMoving(:).*double(MOVING(:))) ./ sumMovingIntensity)) + movingRefObj.YWorldLimits(1);
movingZCOM = (movingRefObj.PixelExtentInWorldZ .* (sum(zMoving(:).*double(MOVING(:))) ./ sumMovingIntensity)) + movingRefObj.ZWorldLimits(1);
translationX = fixedXCOM - movingXCOM;
translationY = fixedYCOM - movingYCOM;
translationZ = fixedZCOM - movingZCOM;

% Coarse alignment
initTform = affine3d();
initTform.T(4,1:3) = [translationX, translationY, translationZ];

% Normalize images
movingInit = mat2gray(MOVING);
fixedInit = mat2gray(FIXED);

% Apply transformation
tform = imregtform(movingInit,movingRefObj,fixedInit,fixedRefObj,'translation',optimizer,metric,'PyramidLevels',3,'InitialTransformation',initTform);
T = tform;
RegisteredImage = imwarp(MOVING, movingRefObj, tform, 'OutputView', fixedRefObj, 'SmoothEdges', true);

% Store spatial referencing object
MOVINGREG.SpatialRefObj = fixedRefObj;

end