function extr2 = RespSort(sg,nPhases)

%% Create range of breathing locations, sort into nPhases resp phases
sg = smooth(sg,100);
limit_low = min(sg);%median(sg) - std(sg);

% so we want resp phases with roughly:
readoutsPerResp = numel(find(sg >= limit_low))/nPhases;

step = 0.0002; %abs(mean(diff(respRange))/99);

% start with max_beamLoc and add readouts until we get to readoutsPerResp
% the last phase will have the least amount of readouts
extr2 = zeros(size(sg))+100;               
limits = zeros(2,nPhases);                  % the limits for plotting
for ms = 1:nPhases
    temp_ind = [];
    if ms == 1
        startbeamLoc = max(sg);
    end
    limits(1,ms) = startbeamLoc;
    if ms < nPhases
        while numel(temp_ind) < readoutsPerResp
            temp_ind = cat(1,temp_ind,find(sg <= startbeamLoc & ...
                sg > (startbeamLoc - step)));
            startbeamLoc = startbeamLoc - step;
        end
    else    % in the last resp phase simply find points that haven't been sorted yet
        temp_ind = find(extr2 == 100 & sg >= limit_low);
        startbeamLoc = limit_low;
    end
    extr2(temp_ind) = ms;
    limits(2,ms) = startbeamLoc;
end

%% optional limit plots

c = lines(nPhases);
if 1
    figure(1);  clf;
    subplot 211;
    plot(sg,'Color','k'); hold on;  
    x = [1 length(sg)];
    for i = 1:nPhases
        I = patch([x fliplr(x)],[limits(1,i) limits(1,i) fliplr([limits(2,i) limits(2,i)])],'k');
        I.FaceColor = c(i,:); I.FaceAlpha = 0.2;
    end
    
    xlim([1 length(sg)])
    xlabel('readout number'); ylabel('sg signal (a.u.)')
    
    subplot 212; hold on
    histogram(sg,200, 'FaceColor','k');
    for i = 1:nPhases
        area([limits(2,i), limits(1,i)], [max(ylim), max(ylim)],'FaceColor',...
            c(i,:),'FaceAlpha',0.2);
    end
    xlabel('sg signal (a.u.)'); ylabel('bin count')
end

set(findall(gcf,'-property','FontSize'),'FontSize',16)
set(findall(gcf,'-property','FontName'),'FontName','Microsoft Yahei')
set(gcf,'Position',[964 53 954 913])