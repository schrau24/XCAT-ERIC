function [x] = CSL1NlCg_4DRadial_SG(x0,param)

% starting point
x=x0;

% line search parameters
maxlsiter = 6;
gradToll = 1e-8;
param.l1Smooth = 1e-15;	
alpha = 0.01;  
beta = 0.6;
t0 = 1 ; 
k = 0;
% clear test
% compute g0  = grad(f(x))
g0 = grad(x,param);
dx = -g0;

% iterations
while(1)

    % backtracking line-search
	f0 = objective(x,dx,0,param);
	t = t0;
    f1 = objective(x,dx,t,param);
	lsiter = 0;
	while (f1 > f0 - alpha*t*abs(g0(:)'*dx(:)))^2 & (lsiter<maxlsiter)
		lsiter = lsiter + 1;
		t = t * beta;
		f1 = objective(x,dx,t,param);
	end

	% control the number of line searches by adapting the initial step search
	if lsiter > 2, t0 = t0 * beta; end 
	if lsiter < 1, t0 = t0 / beta; end
%     test(:,:,k+1)=x(:,:,103,9);
	x = (x + t*dx);

    % print some numbers for debug purposes	
    fprintf('  Update %d: obj: %f, L-S: %d\n',k,f1,lsiter);
    k = k + 1;
    
    % add quick update image
    figure(200); clf;
    imagesc(abs(cat(2,x0(:,:,round(size(x0,3)/2),round(size(x0,4)/2))/max(abs(x0(:))),...
        x(:,:,round(size(x0,3)/2),round(size(x0,4)/2))/max(abs(x(:))))))
    colormap gray; axis off equal
    title(sprintf('  Update %d ',k))
    drawnow;
    
    % stopping criteria (to be improved)
	if (k > param.nite) || (norm(dx(:)) < gradToll), break; end

    %conjugate gradient calculation
	g1 = grad(x,param);
	bk = g1(:)'*g1(:)/(g0(:)'*g0(:)+eps);
	g0 = g1;
	dx =  - g1 + bk* dx;
	
end
return;

function res = objective(x,dx,t,param) %**********************************

% L2-norm part
w=(param.E*(x+t*dx)-param.y).*param.SG;
L2Obj=w(:)'*w(:);

% TV part along time
if param.TVWeight
    w = param.TV*(x+t*dx); 
    TVObj = sum((w(:).*conj(w(:))+param.l1Smooth).^(1/2));
else
    TVObj = 0;
end

% TV part along respiration
if param.TVWeightRes
    w = param.TVRes*(x+t*dx); 
    TVResObj = sum((w(:).*conj(w(:))+param.l1Smooth).^(1/2));
else
    TVResObj = 0;
end

%Wavelet part
if param.L1Weight
    L1Obj=0;
    for jj=1:size(x,4)
        w = wavedec3(x(:,:,:,jj)+t*dx(:,:,:,jj),4,'db4');
        for ii=1:29
            temp=w.dec{ii};
            L1Obj = L1Obj+sum((temp(:).*conj(temp(:))+param.l1Smooth).^(1/2));
        end
    end
else
    L1Obj = 0;
end

res=L2Obj+param.TVWeight*TVObj+param.TVWeightRes*TVResObj+param.L1Weight*L1Obj;

function g = grad(x,param)%***********************************************

% L2-norm part
L2Grad = 2.*(param.E'*((param.E*x-param.y).*param.SG));

% TV part along time
if param.TVWeight
    w = param.TV*x;
    TVGrad = param.TV'*(w.*(w.*conj(w)+param.l1Smooth).^(-0.5));
else
    TVGrad=0;
end

% TV part along respiration
if param.TVWeightRes
    w = param.TVRes*x;
    TVResGrad = param.TVRes'*(w.*(w.*conj(w)+param.l1Smooth).^(-0.5));
else
    TVResGrad=0;
end

%Wavelet part
if param.L1Weight
    for jj=1:size(x,4)
        w = wavedec3(x(:,:,:,jj),4,'db4'); 
        for ii=1:29
            w.dec{ii}=(w.dec{ii}.*(w.dec{ii}.*conj(w.dec{ii})+param.l1Smooth).^(-0.5));
        end
        L1Grad(:,:,:,jj) = waverec3(w);
    end
else
    L1Grad=0;
end

g=L2Grad+param.TVWeight*TVGrad+param.TVWeightRes*TVResGrad+param.L1Weight*L1Grad;
