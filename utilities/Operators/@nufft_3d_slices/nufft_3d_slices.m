function  res = nufft_3d_slices(k,b1)

% Multicoil MCNUFFT3D operator
% Based on the NUFFT toolbox from Jeff Fessler
% Input
% k: k-space trajectory
% b1: coil sensitivity maps

%Li Feng, NYU, 12/18/2017
%Eric Schrauben, 2022/08/22, updated to use nufft_3d 
for tt=1:size(k,4)
	kk=k(:,:,:,tt);
    res.st{tt} = nufft_3d(kk,size(b1,1:3),'radial',1,'gpu',0);
end
res.adjoint = 0;
res.imSize = size(b1,1:3);
res.dataSize = size(kk,2:3);
res.b1 = b1;
res = class(res,'nufft_3d_slices');