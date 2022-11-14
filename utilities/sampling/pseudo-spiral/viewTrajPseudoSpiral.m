function viewTrajPseudoSpiral(k, Sl, matrix, app)
nArms = floor(app.nROperFrame/Sl);
cla(app.TrajectoryPlot, 'reset');
hold(app.TrajectoryPlot,'on');
c = parula(nArms);
ylabel(app.TrajectoryPlot,'k_y');
xlabel(app.TrajectoryPlot,'k_z');
axis(app.TrajectoryPlot,'square')
ylim(app.TrajectoryPlot,[-1 matrix(1)+1])
xlim(app.TrajectoryPlot,[-1 matrix(2)+1])
for i = 1:nArms
    ind = (i-1)*Sl+1:(i*Sl);
    ky = squeeze(k(ind,1));
    kz = squeeze(k(ind,2));
    scatter(app.TrajectoryPlot,kz,ky, '.',...
        'MarkerFaceColor',c(i,:),'MarkerEdgeColor',c(i,:));
    
    pause(0.2)
end