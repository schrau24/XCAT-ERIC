function viewTrajRSOS(k, nSlcs, nStacksPerFrame, app)

cla(app.TrajectoryPlot, 'reset');
hold(app.TrajectoryPlot,'on');
c = parula(nStacksPerFrame);  % 8 stacks plotted

view(app.TrajectoryPlot,[-43 32])
axis(app.TrajectoryPlot,'equal')
xlim(app.TrajectoryPlot,[-.5 .5])
ylim(app.TrajectoryPlot,[-.5 .5])
zlim(app.TrajectoryPlot,[-.5 .5])

% normalize kz here for plotting
if nSlcs > 1
    k(:,:,3) = k(:,:,3)/max(max(k(:,:,3)))-0.5;
end

for i = 1:nStacksPerFrame
    ind = (i-1)*nSlcs+1:(i*nSlcs);
    
    kx = squeeze(k(:,ind,1));
    ky = squeeze(k(:,ind,2));
    kz = squeeze(k(:,ind,3));
    
    hline = plot3(app.TrajectoryPlot,kx,ky,kz, ...
        'color',c(i,:),'LineWidth',1.5);
    xlabel(app.TrajectoryPlot,'k_x');
    ylabel(app.TrajectoryPlot,'k_y');
    zlabel(app.TrajectoryPlot,'k_z');
    pause(0.2)
end
legend(app.TrajectoryPlot,'readout','location','southeast')