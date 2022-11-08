function viewTrajPseudoSpiral(k, nSlcs, matrix,app, fileInfo)

cla(app.TrajectoryPlot, 'reset');
hold(app.TrajectoryPlot,'on');
c = parula(8);  % 8 arms plotted
ylabel(app.TrajectoryPlot,'k_y');
xlabel(app.TrajectoryPlot,'k_z');
axis(app.TrajectoryPlot,'equal')
ylim(app.TrajectoryPlot,[-5 matrix(1)+5])
xlim(app.TrajectoryPlot,[-5 matrix(2)+5])
count = 0;
for i = 1:nSlcs*8
    ky = squeeze(k(i,1));
    kz = squeeze(k(i,2));
    
    if mod(i,nSlcs) == 1
        count = count+1;
    end
    scatter(app.TrajectoryPlot,kz,ky, '.',...
        'MarkerFaceColor',c(count,:),'MarkerEdgeColor',c(count,:));
    
    pause(0.05)
    
    if ~isempty(fileInfo)
        frame = getframe(gcf);
        im{i} = frame2im(frame);
    end
end

if ~isempty(fileInfo)
    for idx = 1:size(im,2)
        [A,map] = rgb2ind(im{idx},256,'nodither');
        if idx == 1
            imwrite(A,map,'gifTest.gif','gif','LoopCount',Inf,'DelayTime',.012);
        else
            imwrite(A,map,'gifTest.gif','gif','WriteMode','append','DelayTime',.012);
        end
    end
end