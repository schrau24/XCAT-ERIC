function outputImg = doBART_Resp(kspace)

% to image space along the readout direction
tmpKSP = (ifft(kspace,[],1));
nPE = size(kspace,1);

for SNR = 1%:5
    clc;
    SNR
    % perform CS reconstruction
    % save local tmp data
    tmpsavename = tempname;
    mkdir(tmpsavename)
    for i_FE = 1:nPE
        tmp = tmpKSP(i_FE,:,:,:,:,:);
        tmpsavename_full = fullfile(tmpsavename,sprintf('slice_%03d',i_FE));
        writecfl(tmpsavename_full,tmp);
    end
    
    tmp_out = tmpKSP(:,:,:,1,:,:);
    
    % create parpool
    delete(gcp('nocreate'))
    parpool(8);
    
    % make a waitbar for the parfor loop
    WaitMessage = parfor_wait(nFrames, 'Waitbar', true);
    
    parfor i_FE = 1:nPE % start of parfor single slice loop
        display(sprintf('CS recon: slice %03d',i_FE));
        WaitMessage.Send;
        
        % load local tmp slice
        tmpsavename_full = fullfile(tmpsavename,sprintf('slice_%03d',i_FE));
        tmp = readcfl(tmpsavename_full);
        
        %         cmdsens = 'ecalib -m1 -I -r20';
        cmdsens = 'caldir 21';
        [L,sensemap] = bart_evalc(cmdsens, sum(sum(tmp,5),6)./sum(sum(tmp~=0+eps,5),6) );
        sensemap_all(i_FE,:,:,:) = sensemap;
        
        % BART: PICS
        %         cmdpics = 'pics -S -d 5';
        cmdpics = 'pics -R T:7:0:0.01 -R T:2048:0:0.01 -i 50 -S -d 5';
        dims_2_bart = [1 2 3 4 6 7 8 9 10 11 12 5];
        tmp = permute(tmp, dims_2_bart);
        % [BART MRI DIMS: READ_DIM,	PHS1_DIM,	0.934PHS2_DIM,	COIL_DIM,	MAPS_DIM,	TE_DIM,	COEFF_DIM,	COEFF2_DIM,	ITER_DIM,	CSHIFT_DIM,	TIME_DIM,	TIME2_DIM,	LEVEL_DIM,	SLICE_DIM,	AVG_DIM
        [L,tmp] = bart_evalc(cmdpics,tmp,sensemap);
        tmp = ipermute(tmp, dims_2_bart);
        %                     tmp = bsxfun(@times,create_checkerboard([1,size(tmp,2),size(tmp,3)]),tmp); % redo checkerboard like in mrecon.Data
        tmp_out(i_FE,:,:,:,:,:) = tmp;
        
    end % end of parfor single slice loop
    delete(gcp('nocreate'))
    % Destroy the WaitMessage
    WaitMessage.Destroy;
    
    % delete local tmp data
    for i_FE = 1:nPE
        if (exist(sprintf(fullfile(tmpsavename,'/slice_%03d.hdr'),i_FE),'file'))
            delete(sprintf(fullfile(tmpsavename,'/slice_%03d.hdr'),i_FE));
        end
        if (exist(sprintf(fullfile(tmpsavename,'/slice_%03d.cfl'),i_FE),'file'))
            delete(sprintf(fullfile(tmpsavename,'/slice_%03d.cfl'),i_FE));
        end
    end
    rmdir(tmpsavename);
    clear tmpKSP tmp sensemap T
end

% perform shift along y
shiftdim=2;
shiftdata = permute(tmp_out,[shiftdim 1:shiftdim-1 shiftdim+1:length(size(tmp_out))]);
pixelshift = ceil(size(tmp_out,2)/2);
s_img_size = size(shiftdata);
shiftdata = [shiftdata(pixelshift+1:end,:);shiftdata(1:pixelshift,:)];
shiftdata = reshape(shiftdata,s_img_size);
shiftdata = ipermute(shiftdata,[shiftdim 1:shiftdim-1 shiftdim+1:length(size(shiftdata))]);

shiftdim=3;
shiftdata = permute(shiftdata,[shiftdim 1:shiftdim-1 shiftdim+1:length(size(shiftdata))]);
pixelshift = ceil(size(shiftdata,1)/2);
s_img_size = size(shiftdata);
shiftdata = [shiftdata(pixelshift+1:end,:);shiftdata(1:pixelshift,:)];
shiftdata = reshape(shiftdata,s_img_size);
shiftdata = ipermute(shiftdata,[shiftdim 1:shiftdim-1 shiftdim+1:length(size(shiftdata))]);

outputImg = squeeze(abs(shiftdata));

% add in functionality to save 4D gif here...