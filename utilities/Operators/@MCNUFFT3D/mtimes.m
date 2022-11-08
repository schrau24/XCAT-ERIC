function ress = mtimes(a,bb)

if a.adjoint,
    % Multicoil non-Cartesian k-space to Cartesian image domain
    % nufft for each coil and time point
    for tt=1:size(bb,5),
        for ch=1:size(bb,4),
            for zz=1:size(bb,3),
                b = squeeze(bb(:,:,zz,ch,tt)).*a.w(:,:,tt);
                res(:,:,zz,ch) = reshape(nufft_adj(b(:),a.st{tt})/sqrt(prod(a.imSize2)),a.imSize(1),a.imSize(2));
            end
        end
        ress(:,:,:,tt)=sum(res.*conj(a.b1),4)./sum(abs((a.b1)).^2,4); clear res
    end
else
    % Cartesian image to multicoil non-Cartesian k-space
    for tt=1:size(bb,4),
        for ch=1:size(a.b1,4),
            for zz=1:size(a.b1,3),
                res=bb(:,:,zz,tt).*a.b1(:,:,zz,ch); 
                ress(:,:,zz,ch,tt) = reshape(nufft(res,a.st{tt})/sqrt(prod(a.imSize2)),a.dataSize(1),a.dataSize(2)).*a.w(:,:,tt);
            end
        end
    end
end

