function plot_timetrace(Mean,Std,type,input,sessionNum,RorL)
% type -1) cell : cell average across trials
%	2) trial: trial average across cells
%	3) session: session average across trials and cells
 % input - corresponding number associated with type
% RorL - 1 is R, 2 is L
 figure
    plot(squeeze(Mean),'k')
    hold on
    Mean = reshape(Mean,[1 76]); Std = reshape(Std,[1 76]);
    h =fill([1:76, fliplr(1:76)],[Mean+Std , fliplr(Mean-Std)],'k');
    h.FaceAlpha = 0.2; h.EdgeColor = 'none'; h.EdgeAlpha = 0;
    xL = get(gca,'XLim');
    yL = get(gca,'YLim');
    if RorL==1
        Dir = 'R';
    else
        Dir = 'L';
    end
    
    if isequal(type,'session') 
    axis([0 max(xL) [yL]]);
    else
    axis([0 max(xL) -0.002 max(yL)]);
    end

    if isequal(type,'trial')
    Label= {sprintf('cell %d',input)};
    title(sprintf('Average across trial in Session %d (Stimulus-%s)',sessionNum, Dir))
 
    elseif isequal(type,'cell')
    Label= {sprintf('trial %d',input)};
    title(sprintf('Population average in Session %d (Stimulus-%s)',sessionNum,Dir))
   
    elseif  isequal(type,'session')
    Label= {sprintf('session %d',sessionNum)};
    title(sprintf('Average in Session %d (Stimulus-%s)',sessionNum, Dir))

    end

    legend(Label)
    xlabel('Time frame')
    plot(xL,[0 0],'Color','k','LineStyle',':')  % trial onset/cue onset??
    plot([13 13],yL,'Color','k','LineStyle',':')        % trial onset/cue onset??
%    plot([14 14],yL,'Color','k','LineStyle',':')         
%    area([39-12 39+12],[max(yL) max(yL)],'FaceColor','y','FaceAlpha',0.2, ...
%       'EdgeColor','none','EdgeAlpha',0)

    plot([39 39],yL,'Color','k','LineStyle',':')        % cue offset / delay period onset
%    area([54-12 54+12],[max(yL) max(yL)],'FaceColor','b','FaceAlpha',0.2, ...
%       'EdgeColor','none','EdgeAlpha',0)
    plot([54 54],yL,'Color','k','LineStyle',':')        % trial end / turn onset
    ylabel('Spike Activity (a.u)')
    if isequal(type,'trial')
    saveas(gcf,sprintf('AverageAcrossTrials_cell%d_%s.png',input,Dir))
 
    elseif isequal(type,'cell')
    saveas(gcf,sprintf('AverageAcrossCells_trial%d_%s.png',input,Dir))

    elseif isequal(type,'session')

    saveas(gcf,sprintf('AverageInSession%d_%s.png',input,Dir))
    end
end
