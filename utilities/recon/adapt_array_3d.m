function [recon,cmap]=adapt_array_3d(yn,rn,norm)

% adapted from adapt_array_2d

yn=permute(yn,[4,1,2,3]);
[nc,nx,ny,nz]=size(yn);
if nargin<3, norm=0; end
if nargin<2, rn=eye(nc);end

% find coil with maximum intensity for phase correction
[mm,maxcoil]=max(sum(sum(sum(permute(abs(yn),[4 3 2 1])))));

% based on image size
if nx > 162 
    bs = 8; st = 4;
    if nz/2 < st
        st = round(nz/2);
    end
else
    bs = 4; st = 2;
end

% make sure  nz/st is integer
if mod(nz,st)~=0
    st = 2;
end

bs1=bs;  %x-block size
bs2=bs;  %y-block size
bs3=bs;  %z-block size
% st=4;   %increase to set interpolation step size

wsmall=zeros(nc,round(nx./st),round(ny./st),round(nz./st));
cmapsmall=zeros(nc,round(nx./st),round(ny./st),round(nz./st));

if nz==1
    stz = 1;
else
    stz = st;
end

for z=stz:stz:nz
    for x=st:st:nx
        for y=st:st:ny
            
            
            %Collect block for calculation of blockwise values
            ymin1=max([y-bs1./2 1]);
            xmin1=max([x-bs2./2 1]);
            zmin1=max([z-bs3./2 1]);
            % Cropping edges
            ymax1=min([y+bs1./2 ny]);
            xmax1=min([x+bs2./2 nx]);
            zmax1=min([z+bs3./2 nz]);
            
            ly1=length(ymin1:ymax1);
            lx1=length(xmin1:xmax1);
            lz1=length(zmin1:zmax1);
            m1=reshape(yn(:,xmin1:xmax1,ymin1:ymax1,zmin1:zmax1),nc,lx1*ly1*lz1);
            
            
            m=m1*m1'; %signal covariance
            
            % eignevector with max eigenvalue for optimal combination
            [e,v]=eig(inv(rn)*m);
            
            v=diag(v);
            [mv,ind]=max(v);
            
            mf=e(:,ind);
            mf=mf/(mf'*inv(rn)*mf);
            normmf=e(:,ind);
            
            % Phase correction based on coil with max intensity
            mf=mf.*exp(-1i*angle(mf(maxcoil)));
            normmf=normmf.*exp(-1i*angle(normmf(maxcoil)));
            
            wsmall(:,x./st,y./st,z./stz)=mf;
            cmapsmall(:,x./st,y./st,z./stz)=normmf;
        end
    end
end

sdky=linspace(1,round(ny./st),ny);
sdkx=linspace(1,round(nx./st),nx);
sdkz = linspace(1,round(nz./st),nz);
[Xq,Yq,Zq] = meshgrid(sdky,sdkx,sdkz);

% Interpolation of weights upto the full resolution
% Done separately for magnitude and phase in order to avoid 0 magnitude
% pixels between +1 and -1 pixels.
wfull = zeros(nc,nx,ny,nz);
cmap = zeros(nc,nx,ny,nz);

for i=1:nc
    if nz==1
        wfull(i,:,:,:)= conj(interp2(squeeze(abs(wsmall(i,:,:,:))),Xq,Yq,'bilinear').*exp(1i.*interp2(squeeze(angle(wsmall(i,:,:,:))),Xq,Yq,'nearest')));
        cmap(i,:,:,:) = interp2(squeeze(abs(cmapsmall(i,:,:,:))),Xq,Yq,'bilinear').*exp(1i.*interp2(squeeze(angle(cmapsmall(i,:,:,:))),Xq,Yq,'nearest'));
    else
        wfull(i,:,:,:)= conj(interp3(squeeze(abs(wsmall(i,:,:,:))),Xq,Yq,Zq,'bilinear').*exp(1i.*interp3(squeeze(angle(wsmall(i,:,:,:))),Xq,Yq,Zq,'nearest')));
        cmap(i,:,:,:) = interp3(squeeze(abs(cmapsmall(i,:,:,:))),Xq,Yq,Zq,'bilinear').*exp(1i.*interp3(squeeze(angle(cmapsmall(i,:,:,:))),Xq,Yq,Zq,'nearest'));
    end
end
recon=squeeze(sum(wfull.*yn));   %Combine coil signals.
% normalization proposed in the abstract by Griswold et al.
if norm, recon=recon.*squeeze(sum(abs(cmap))).^2; end

cmap=permute(cmap,[2,3,4,1]);
