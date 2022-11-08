function viewTrajRSOS(k, nSlcs, app, fileInfo)

cla(app.TrajectoryPlot, 'reset');
hold(app.TrajectoryPlot,'on');
c = parula(5);  % 8 stacks plotted

view(app.TrajectoryPlot,[-45 17])
axis(app.TrajectoryPlot,'equal')
xlim(app.TrajectoryPlot,[-.5 .5])
ylim(app.TrajectoryPlot,[-.5 .5])
zlim(app.TrajectoryPlot,[-.5 .5])

% normalize kz here for plotting
k(:,:,3) = k(:,:,3)/max(max(k(:,:,3)))-0.5;

count = 0;
for i = 1:nSlcs*5
    kx = squeeze(k(:,i,1));
    ky = squeeze(k(:,i,2));
    kz = squeeze(k(:,i,3));
    
    if mod(i,nSlcs) == 1
        count = count+1;
    end
    hline = plot3(app.TrajectoryPlot,kx,ky,kz, ...
        'color',c(count,:),'LineWidth',1.5);
    xlabel(app.TrajectoryPlot,'k_x');
    ylabel(app.TrajectoryPlot,'k_y');
    zlabel(app.TrajectoryPlot,'k_z');
%     drawnow
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