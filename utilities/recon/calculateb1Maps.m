function [b1, recon] = calculateb1Maps(app, kspace)
nSlcs = length(unique(squeeze(app.kz_samples(1,:))));
% build k matrix from all currSpokes
x = app.ky_samples;
y = app.kx_samples;
Traj = x+1i*y;
Traj = Traj(:,1:nSlcs:end);
DensityComp = sqrt(x.^2+y.^2);
DensityComp = DensityComp(:,1:nSlcs:end);
DensityComp = DensityComp/max(DensityComp(:));

for sl = 1:nSlcs
    kdata(:,:,sl,:) = permute(kspace(:,sl:nSlcs:end,:),[1 2 4 3]);
end

[~,kz_ord] = sort(app.kz_samples(1,1:nSlcs));
kdata = kdata(:,:,kz_ord,:);
% pad kdata by zeros for half-scan
kdata = padarray(kdata,[0 0 app.matrixFH-nSlcs 0], 'pre');

kspaceSorted_slices = ifft(kdata,[],3);
% kdata1 = ifft(kdata,[],3);
shiftdim=3;
shiftdata = permute(kspaceSorted_slices,[shiftdim 1:shiftdim-1 shiftdim+1:length(size(kspaceSorted_slices))]);
pixelshift = ceil(size(shiftdata,1)/2);
s_img_size = size(shiftdata);
shiftdata = [shiftdata(pixelshift+1:end,:);shiftdata(1:pixelshift,:)];
shiftdata = reshape(shiftdata,s_img_size);
kdata1 = ipermute(shiftdata,[shiftdim 1:shiftdim-1 shiftdim+1:length(size(shiftdata))]);

[nx,ntviews,nz,nc]=size(kdata1);
filter_factor=20;
filter=kaiser(nx,filter_factor);
% kdata1=double(kdata1.*repmat(filter,[1,ntviews,nz,nc]));
kdata1=double(kdata1.*repmat(sqrt(DensityComp),[1,1,nz,nc]));

Img_Dim=[app.matrixRL, app.matrixRL, 1];
param.E = MCNUFFT3D(double(Traj),double(DensityComp),ones(Img_Dim));
clear ref;
for ch=1:nc
    ref(:,:,:,ch)=param.E'*kdata1(:,:,:,ch);
end
ref=ref/max(abs(ref(:)));

[recon,b1] = adapt_array_3d(ref,size(kdata1,4),1);
b1 = b1/max(abs(b1(:)));