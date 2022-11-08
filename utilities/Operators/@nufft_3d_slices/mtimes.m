function ress = mtimes(a,bb)

if a.adjoint
    % Multicoil non-Cartesian k-space to Cartesian image domain
    % nufft for each coil and time point
    maxit = 1; % 0 or 1 for gridding, higher values for conjugate gradient
    D = 0; % Tikhonov penalty on ||x||
    for tt=1:size(bb,5)
        b = squeeze(bb(:,:,:,:,tt));
        %%reconstruction (inverse transform)
        partial = 0.5; % Tikhobov penalty on ||imag(x))||
        tmp = a.st{tt}.iNUFT(b,maxit,D); % plain reconstruction
        ress(:,:,:,tt) = sum(tmp.*conj(a.b1),4)./sum(abs(a.b1).^2,4);
    end
    % cheat here and set nans to zero
    ress(isnan(ress)) = 0;
else
    % Cartesian image to multicoil non-Cartesian k-space
    for tt=1:size(bb,4)
        for ch=1:size(a.b1,4)
            res=bb(:,:,:,tt).*a.b1(:,:,:,ch);
            kk = a.st{tt}.fNUFT(res);
            ress(1,:,:,ch,tt) = reshape(kk, a.dataSize(1), []);%nufft(res,a.st{tt})/sqrt(prod(a.imSize2)),a.dataSize(1),a.dataSize(2)).*a.w(:,:,tt);
        end
    end
end

