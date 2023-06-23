function saveSimToGIF(app)

imSize = size(app.IMG_CP);

contrast = app.ImageContrastSlider.Value;

figure(1);clf
set(gcf, 'color' ,'w')
dat = squeeze(app.IMG_CP(app.currentX_GT,app.currentY_GT,...
    app.currentSl_GT,:)); dat = dat/max(dat(:));
plotTimes = app.timing;

filename = 'tmp.gif';

for frame = 1:size(app.IMG_CP,4)
    % axial
    subplot 221
    aa = squeeze(app.IMG_CP(:,:,app.currentSl_GT,frame));
    imshow(aa,[0 contrast]);
    images.roi.Crosshair(gca, 'Position',[app.currentY_GT app.currentX_GT], ...
        'LineWidth', 1, 'Color', 'r');
    title('axial')
    
    % sagittal
    subplot 222
    aa = squeeze(app.IMG_CP(:,app.currentY_GT,:,frame));
    imshow(flipud(rot90(aa,-1)),[0 contrast]);
    images.roi.Crosshair(gca, 'Position',[imSize(1)-app.currentX_GT imSize(3)-app.currentSl_GT], ...
        'LineWidth', 1, 'Color', 'r');
    daspect([app.footheadRes.Value/app.anteriorposteriorRes.Value imSize(2)/imSize(1) 1]);
    title('sagittal')
    
    % coronal
    subplot 223
    aa = squeeze(app.IMG_CP(app.currentX_GT,:,:,frame));
    imshow(flipud(rot90(aa,-1)),[0 contrast]);
    images.roi.Crosshair(gca,'Position',[imSize(2)-app.currentY_GT imSize(3)-app.currentSl_GT], ...
        'LineWidth', 1, 'Color', 'r');
    daspect([app.footheadRes.Value/app.leftrightRes.Value 1 1]);
    title('coronal')
    
    % time plot
    subplot 224
    plot(plotTimes,dat,'k','LineWidth',2);
    xlim([0 max(plotTimes)]);
    hold on;
    xline(plotTimes(frame),'r','LineWidth',1);
    hold off;
    box off
    xlabel(gca, 'time (s)');
    title('relative voxel signal evolution')
    
    %    % Take a "screenshot" of the figure fh
    ff = getframe(gcf);
    % Turn screenshot into image
    im = frame2im(ff);
    % Turn image into indexed image (the gif format needs this)
    [imind,cm] = rgb2ind(im,256);
              
    if frame == 1
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
    else
        imwrite(imind,cm, filename,'gif','WriteMode','append','DelayTime',0.01);
    end
end
